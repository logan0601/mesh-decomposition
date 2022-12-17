#include "solver.h"
#include <algorithm>
#include <omp.h>


Solver::Solver(Model& model_) : model(model_)
{
    lnum = 0;
    lmat = vector<vector<int>>();
}


void Solver::solve()
{
    lnum = 0;
    labels = vector<int>( model.fnum, 0 );
    lmat = vector<vector<int>>( 1, vector<int>( model.fnum, 0 ) );
    for(int i = 0; i < model.fnum; i++)
        lmat[0][i] = i;

    double s = omp_get_wtime();
    kernel(0);
    double e = omp_get_wtime();
    printf("Using %.2lf seconds\n", e - s);

    model.lnum = lnum;
    model.labels = labels;
}


void Solver::kernel(int l)
{
    Pair_vd p = init_seed(l);
    if(check(l, p.first, p.second))
        return;
    assign(l, p.first);
}


Pair_vd Solver::init_seed(int l)
{
    vector<int>& id = lmat[l];
    int n = id.size();

    vector<double> distant(n, 0);
    #pragma omp parallel for
    for(int i = 0; i < n; i++)
    {
        double sum = 0;
        #pragma omp parallel for reduction(+:sum)
        for(int j = 0; j < n; j++)
            sum += model.dmat[id[i]][id[j]];
        distant[i] = sum;
    }
    int sid = std::distance(distant.begin(), std::min_element(distant.begin(), distant.end()));

    vector<int> vis( n, 0 ), seed;
    vector<double> dis( n, Utils::UPPER ), dif( Utils::NUM_SEED, 0 );

    vis[sid] = 1;
    seed.push_back(sid);
    for(int i = 1; i < Utils::NUM_SEED; i++)
    {
        #pragma omp parallel for
        for(int j = 0; j < n; j++)
            dis[j] = std::min(dis[j], model.dmat[id[j]][id[sid]]);
        dis[id[sid]] = -1;

        sid = std::distance(dis.begin(), std::max_element(dis.begin(), dis.end()));
        vis[sid] = 1;
        seed.push_back(sid);
        dif[i] = dis[sid];
    }
    for(int i = 1; i < Utils::NUM_SEED-1; i++)
        dif[i] -= dif[i+1];
    int snum = std::distance(dif.begin(), std::max_element(dif.begin(), dif.end()-1));

    Pair_vd p( vector<int>(), 0 );
    for(int i = 0; i <= snum; i++) p.first.push_back(seed[i]);
    for(int i = 0; i < n; i++) p.second += distant[i];
    p.second /= n * (n - 1);
    return p;
}


bool Solver::check(int l, vector<int>& seed, double avg_dis)
{
    vector<int>& id = lmat[l];
    int snum = seed.size();
    int fnum = id.size();

    if(fnum < 50)
        return true;

    if((avg_dis / model.avg_dis) < Utils::AVG_EPS)
        return true;
    
    vector<double> dis(snum, 0);
    for(int i = 0; i < snum; i++)
    {
        #pragma omp parallel for
        for(int j = 0; j < snum; j++)
            dis[j] = std::max(dis[j], model.dmat[id[seed[j]]][id[seed[i]]]);
    }
    double max_dis = *std::max_element(dis.begin(), dis.end());
    if(max_dis < Utils::DIS_EPS)
        return true;
    
    double angles[2] = { 0, acos(-1) };
    for(int i = 0; i < snum; i++)
    {
        for(const Edge& e : model.faces[i].nbrs)
        {
            angles[0] = std::max(angles[0], e.angle);
            angles[1] = std::min(angles[1], e.angle);
        }
    }
    if((angles[0] - angles[1]) < Utils::ANG_EPS)
        return true;
    
    return false;
}


void Solver::assign(int l, vector<int>& seed)
{
    vector<int>& id = lmat[l];
    int n = id.size();
    int snum = seed.size();

    vector<int> nlabel(n, -1);
    for(int i = 0; i < snum; i++)
        nlabel[seed[i]] = i;
    
    vector<vector<double>> pmat( n, vector<double>(snum, 0) );
    #pragma omp parallel for
    for(int i = 0; i < n; i++)
    {
        if(nlabel[i] == -1)
        {
            double tot = 0;
            for(int j = 0; j < snum; j++) pmat[i][j] = 1.0 / model.dmat[id[i]][id[seed[j]]];
            for(int j = 0; j < snum; j++) tot += pmat[i][j];
            for(int j = 0; j < snum; j++) pmat[i][j] /= tot;
        }
    }

    int p[3];
    std::map<Pair_ii, vector<int>> fuzzy;
    for(int i = 0; i < n; i++)
    {
        if(nlabel[i] != -1) continue;

        vector<double>::iterator s = pmat[i].begin();
        vector<double>::iterator e = pmat[i].end();
        p[0] = std::distance(s, std::max_element(s, e));
        if(p[0] == 0) p[1] = 1;
        else p[1] = std::distance(s, std::max_element(s, s+(p[0]-1)));
        if(p[0] == snum-1) p[2] = snum - 2;
        else p[2] = std::distance(s, std::max_element(s+p[0]+1, e));
        p[1] = pmat[i][p[1]] < pmat[i][p[2]] ? p[2] : p[1];

        if(pmat[i][p[0]] > 0.5 || (pmat[i][p[0]] - pmat[i][p[1]]) > Utils::PROB_EPS)
            nlabel[i] = p[0];
        else
        {
            Pair_ii t = sort( p[0], p[1] );
            if(fuzzy.find(t) == fuzzy.end()) fuzzy[t] = vector<int>();
            fuzzy[t].push_back(i);
        }
    }

    for(int i = 0; i < n; i++)
        if(nlabel[i] != -1)
            model.faces[id[i]] = lnum + nlabel[i];
    
    assign_fuzzy(l, seed, fuzzy);
}


void Solver::assign_fuzzy(int l, vector<int>& seed, std::map<Pair_ii, vector<int>>& fuzzy)
{
    vector<int>& id = lmat[l];
    int n = fuzzy.size();
    int snum = seed.size();

    
}
