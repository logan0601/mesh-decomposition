#include "solver.h"
#include <omp.h>
#include <sstream>


Solver::Solver(Model& mod_, int lev_, vector<int> fid_) : model(mod_)
{
    level = lev_;

    if(fid_.size() == 0) {
        for(int i = 0; i < model.fnum; i++)
            fid.push_back(i);
    }
    else fid = fid_;
    fnum = fid.size();
    old_lab = model.faces[fid[0]].label;

    dmat = vector<vector<double>>(fnum, vector<double>(fnum, 0.0));
    avg_dis = 0;
    for(int i = 0; i < fnum; i++) {
        for(int j = 0; j < fnum; j++) {
            dmat[i][j] = model.dmat[fid[i]][fid[j]];
            avg_dis += dmat[i][j];
        }
    }
    avg_dis /= fnum * (fnum - 1);

    rid = vector<int>();
    rnum = 0;
}


void Solver::initial_repre()
{
    vector<double> tot(fnum, 0);
    #pragma omp parallel for
    for(int i = 0; i < fnum; i++)
    {
        double sum = 0;
        #pragma omp parallel for reduction(+:sum)
        for(int j = 0; j < fnum; j++) sum += dmat[i][j];
        tot[i] = sum;
    }
    int id = distance(tot.begin(), min_element(tot.begin(), tot.end()));

    vector<int> rep;
    vector<double> dis(fnum, INF), dif(NREP, 0);
    dis[id] = -1;
    rep.push_back(id);
    for(int i = 0; i < NREP; i++)
    {
        #pragma omp parallel for
        for(int j = 0; j < fnum; j++)
            dis[j] = std::min(dis[j], dmat[j][id]);
        id = distance(dis.begin(), max_element(dis.begin(), dis.end()));
        rep.push_back(id);
        dif[i] = dis[id];
        dis[id] = -1;
    }

    for(int i = 0; i < NREP - 1; i++)
        dif[i] -= dif[i + 1];
    dif[NREP-1] = 0;

    int num = std::distance(dif.begin(), std::max_element(dif.begin(), dif.end()));
    rid = vector<int>();
    unique = vector<int>();
    for(int i = 0; i <= num+1; i++) {
        rid.push_back(rep[i]);
        unique.push_back(rep[i]);
    }
    rnum = num + 2;
}


void Solver::compute_prob_matrix()
{
    pmat = vector<vector<double>>(fnum, vector<double>(rnum, 0));
    vector<int> rep(fnum, 0);
    int num = unique.size();
    for(int i = 0; i < num; i++) {
        rep[unique[i]] = 1;
        pmat[unique[i]][i] = 1;
    }
    #pragma omp parallel for
    for(int i = 0; i < fnum; i++) {
        if(rep[i] == 1)
            continue;
        double sum = 0;
        for(int j = 0; j < num; j++) pmat[i][j] = 1.0 / dmat[i][unique[j]];
        for(int j = 0; j < num; j++) sum += pmat[i][j];
        for(int j = 0; j < num; j++) pmat[i][j] /= sum;
    }
}


void Solver::assign()
{
    int p[3], num = unique.size();
    double eps = (num <= 3) ? 0.04 : 0.02;
    for(int i = 0; i < fnum; i++) {
        auto s = pmat[i].begin(), e = pmat[i].end();
        p[0] = distance(s, max_element(s, e));
        p[1] = (p[0] == 0) ? 1 : distance(s, max_element(s, s + p[0]));
        p[2] = (p[0] == rnum - 1) ? rnum - 2 : distance(s, max_element(s + (p[0] + 1), e));
        p[1] = pmat[i][p[1]] < pmat[i][p[2]] ? p[2] : p[1];

        Face& f = model.faces[fid[i]];
        if(pmat[i][p[0]] - pmat[i][p[1]] >= eps) f.label = lab_idx + p[0];
        else f.label = fuz_idx + p[0] * rnum + p[1];
    }
}


void Solver::assign_fuzzy()
{
    for(int i = 0; i < rnum; i++) {
        for(int j = 0; j < rnum; j++) {
            if(j <= i)
                continue;
            vector<int> ftype(model.fnum, 0);
            for(int k : fid) {
                const Face& f = model.faces[k];
                if(f.label == fuz_idx + i * rnum + j || f.label == fuz_idx + j * rnum + i) {
                    ftype[k] = 3;
                    for(const Edge& e : f.edges) {
                        const Face& n = model.faces[e.fid];
                        if(n.label == lab_idx + i) ftype[e.fid] = 1;
                        if(n.label == lab_idx + j) ftype[e.fid] = 2;
                    }
                }
            }
            model.compute_flow(ftype);

            for(int k = 0; k < model.fnum; k++) {
                if(ftype[k] == 1) model.faces[k].label = lab_idx + i;
                if(ftype[k] == 2) model.faces[k].label = lab_idx + j;
            }
        }
    }  
}


vector<int> Solver::recompute_repre()
{
    assign();
    vector<vector<double>> dis(rnum, vector<double>(fnum, 0));
    vector<int> cnt(rnum, 0);
    #pragma omp parallel for
    for(int i = 0; i < fnum; i++) {
        int k = model.faces[fid[i]].label - lab_idx;
        if(k >= rnum)
            continue;
        cnt[k]++;
        #pragma omp parallel for
        for(int j = 0; j < fnum; j++) dis[k][j] += dmat[i][j];
    }
    for(int i = 0; i < rnum; i++)
        for(int j = 0; j < fnum; j++)
            dis[i][j] = (cnt[i] == 0) ? INF : dis[i][j] / cnt[i];
    #pragma omp parallel for
    for(int i = 0; i < fnum; i++) {
        double sum = 0.0;
        for(int j = 0; j < rnum; j++) pmat[i][j] = 1 / (dis[j][i] + EPS);
        for(int j = 0; j < rnum; j++) sum += pmat[i][j];
        for(int j = 0; j < rnum; j++) pmat[i][j] /= sum;
    }

    cmat = vector<vector<double>>(rnum, vector<double>(fnum, 0));
    for(int i = 0; i < rnum; i++) {
        #pragma omp parallel for
        for(int j = 0; j < fnum; j++) {
            double sum = 0;
            #pragma omp parallel for reduction(+:sum)
            for(int k = 0; k < fnum; k++)
                sum += pmat[k][i] * dmat[j][k];
            cmat[i][j] = sum;
        }
    }

    vector<int> rep(rnum, 0);
    for(int i = 0; i < rnum; i++)
        rep[i] = distance(cmat[i].begin(), min_element(cmat[i].begin(), cmat[i].end()));
    return rep;
}


bool Solver::check()
{
    if(fnum < 200)
        return true;

    if(level >= LEVEL)
        return true;

    if((avg_dis / model.avg_dis) < 0.3)
        return true;
    
    return false;
}


void Solver::solve()
{
    double start = omp_get_wtime();
    initial_repre();

    lab_idx = model.lnum;
    fuz_idx = model.lnum + rnum;
    printf("[INIT] Repre Number %d\n", rnum);
    vector<int> rep;
    for(int i = 0; i < NREP; i++) {
        compute_prob_matrix();
        rep = recompute_repre();
        bool change = false;
        for(int j = 0; j < rnum; j++)
            if(rep[j] != rid[j]) {
                change = true;
                break;
            }
        if(change) {
            rid = rep;
            unique = vector<int>();
            std::map<int, int> rmap;
            for(int id : rid) rmap[id] = 1;
            for(auto& p : rmap) unique.push_back(p.first);
        }
        else break;
    }
    recompute_repre();
    assign();
    assign_fuzzy();

    vector<int> cnt(rnum, 0);
    for(const Face& f : model.faces)
        if(lab_idx <= f.label && f.label < fuz_idx) cnt[f.label - lab_idx] ++;
    
    printf("{%d: %d} -> { ", old_lab, fnum);
    for(int i = 0; i < rnum; i++)
        printf("%d: %d, ", lab_idx + i, cnt[i]);
    printf("}\n");

    model.lnum += rnum;

    if(check()) return;
    vector<Solver> solvers;
    for(int i = 0; i < rnum; i++) {
        vector<int> id;
        for(int j = 0; j < model.fnum; j++)
            if(model.faces[j].label == lab_idx + i)
                id.push_back(j);
        if(id.size() > 0)
            solvers.push_back(Solver(model, level+1, id));
    }
    for(Solver& solver : solvers)
        solver.solve();

    double end = omp_get_wtime();
    printf("Using %.2lf seconds\n", end - start);
}
