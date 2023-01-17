#include "model.h"

#include <fstream>
#include <string>
#include <string.h>
#include <omp.h>


pii sort(int a, int b) { return (a < b) ? pii(a, b) : pii(b, a); }


Model::Model() { vnum = 0; fnum = 0; lnum = 0; }


void Model::parse(const char* path)
{
    std::ifstream in;
    in.open(path, std::fstream::in);

    std::string line;
    while(getline(in, line)) {
        if(line.find("element vertex") == 0) vnum = std::stoi(line.substr(15));
        if(line.find("element face") == 0)   fnum = std::stoi(line.substr(13));
        if(line.find("end_header") == 0) break;
    }

    double x[3];
    int p[4] = {0}, k[3];
    for(int i = 0; i < vnum; i++) {
        getline(in, line);
        for(int j = 1; j < 3; j++) p[j] = line.find(" ", p[j-1]) + 1;
        for(int j = 0; j < 3; j++) x[j] = std::stod(line.substr(p[j]));
        vecs.push_back(Vec(x[0], x[1], x[2]));
    }
    for(int i = 0; i < fnum; i++) {
        getline(in, line);
        for(int j = 1; j < 4; j++) p[j] = line.find(" ", p[j-1]) + 1;
        for(int j = 0; j < 3; j++) k[j] = std::stoi(line.substr(p[j+1]));
        faces.push_back(Face(k[0], k[1], k[2], vecs[k[0]], vecs[k[1]], vecs[k[2]]));
    }
    in.close();
}


void Model::build_edge(int i, int j, pii vid)
{
    Face& f1 = faces[i];
    Face& f2 = faces[j];

    double angle = acos(Vec::dot(f1.normal, f2.normal));
    bool is_convex = Vec::dot(f1.normal, f2.center-f1.center) < EPS;
    double eta = is_convex ? ETA : 1;
    double ang_dis = eta * (1 - Vec::dot(f1.normal, f2.normal));

    Vec ed = vecs[vid.second] - vecs[vid.first];
    Vec v1 = f1.center - vecs[vid.first];
    Vec v2 = f2.center - vecs[vid.first];
    double elen = ed.len();
    double len1 = v1.len();
    double len2 = v2.len();
    double angle1 = acos(Vec::dot(ed, v1) / elen / len1);
    double angle2 = acos(Vec::dot(ed, v2) / elen / len2);
    double geo_dis = v1.len2() + v2.len2() - 2 * len1 * len2 * cos(angle1 + angle2);

    f1.edges.push_back(Edge(vid.first, vid.second, j, angle, ang_dis, geo_dis));
    f2.edges.push_back(Edge(vid.first, vid.second, i, angle, ang_dis, geo_dis));
}


void Model::build_graph()
{
    int emap[3][2] = {{0,1}, {0,2}, {1,2}};
    std::map<pii, int> visit;
    for(int i = 0; i < fnum; i++) {
        const Face& f = faces[i];
        for(int j = 0; j < 3; j++) {
            pii p = sort(f.vid[emap[j][0]], f.vid[emap[j][1]]);
            if(visit.find(p) != visit.end()) {
                build_edge(visit[p], i, p);
            }
            else visit[p] = i;
        }
    }

    int count = 0;
    double ang_dis = 0, geo_dis = 0;
    for(const Face& f : faces) {
        count += f.edges.size();
        for(const Edge& e : f.edges) {
            ang_dis += e.ang_dis;
            geo_dis += e.geo_dis;
        }
    }

    avg_ang_dis = ang_dis / count;
    double avg_geo_dis = geo_dis / count;
    for(Face& f : faces) {
        for(Edge& e : f.edges)
            e.tot_dis = DELTA * e.ang_dis / avg_ang_dis + (1-DELTA) * e.geo_dis / avg_geo_dis;
    }
}


void Model::build_short_path()
{
    vector<int> visit;
    vector<double> distant;
    std::priority_queue<pdi, vector<pdi>, std::greater<pdi>> q;
    dmat = vector<vector<double>>(fnum, vector<double>(fnum, 0));

    #pragma omp parallel for private(visit) private(distant) private(q)
    for(int i = 0; i < fnum; i++) {
        visit = vector<int>(fnum, 0);
        distant = vector<double>(fnum, INF);
        q = std::priority_queue<pdi, vector<pdi>, std::greater<pdi>>();

        distant[i] = 0;
        q.push(pdi(0, i));
        while(!q.empty()) {
            int u = q.top().second;
            q.pop();

            if(visit[u]) continue;
            visit[u] = 1;
            for(const Edge& e : faces[u].edges) {
                int v = e.fid;
                double w = e.tot_dis;
                if(distant[u] + w < distant[v]) {
                    distant[v] = distant[u] + w;
                    q.push(pdi(distant[v], v));
                }
            }
        }
        #pragma omp parallel for
        for(int j = 0; j < fnum; j++) dmat[i][j] = distant[j];
    }
    avg_dis = 0;
    for(int i = 0; i < fnum; i++)
        for(int j = 0; j < fnum; j++) avg_dis += dmat[i][j];
    avg_dis /= fnum * (fnum - 1);
}


void Model::read(const char* path)
{
    double start = omp_get_wtime();
    parse(path);
    build_graph();
    build_short_path();
    double end = omp_get_wtime();
    printf("Using %.2lf seconds\n", end - start);
}


void Model::write(const char* path)
{
    std::fstream os;
    os.open(path, std::fstream::out);
    os << std::fixed;
    os << std::showpoint;
    os.precision(4);
    os << "ply\n";
    os << "format ascii 1.0\n";
    os << "element vertex " << vnum << "\n";
    os << "property float x\n";
    os << "property float y\n";
    os << "property float z\n";
    os << "element face " << fnum << "\n";
    os << "property list uchar int vertex_indices\n";
    os << "property uint8 red\n";
    os << "property uint8 green\n";
    os << "property uint8 blue\n";
    os << "end_header\n";
    for(const Vec& v : vecs)
        os << v.x << " " << v.y << " " << v.z << "\n";
    for(const Face& f : faces)
        os << "3 " << f.vid[0] << " " << f.vid[1] << " " << f.vid[2] << " " <<
            (60*(f.label % 4 + 1)) << " " << (80*((f.label + 1) % 3 + 1)) << " " << (50*((f.label + 2) % 5 + 1)) << "\n";
    os.close();
}


Flow::Flow(int s_, int t_, double cap_, double flow_) { s = s_; t = t_; cap = cap_; flow = flow_; }


void Model::compute_flow(vector<int>& ftype)
{
    vector<vector<Flow>> flows(fnum + 2, vector<Flow>());
    int start = fnum;
    int end = fnum + 1;

    for(int i = 0; i < fnum; i++) {
        if(ftype[i] == 0)
            continue;
        for(const Edge& e : faces[i].edges)
            if(ftype[e.fid] != 0)
                flows[i].push_back(Flow(i, e.fid, 1.0 / (1.0 + e.ang_dis / avg_ang_dis), 0));
        if(ftype[i] == 1)
            flows[start].push_back(Flow(start, i, INF, 0));
        if(ftype[i] == 2)
            flows[i].push_back(Flow(i, end, INF, 0));
    }

    vector<int> p, h;
    vector<double> a;
    std::queue<int> q;
    while(true) {
        p = vector<int>( fnum + 2, -1 );
        h = vector<int>( fnum + 2, 0 );
        a = vector<double>( fnum + 2, 0 );
        q = std::queue<int>();

        q.push(start);
        a[start] = INF;
        while(!q.empty()) {
            int u = q.front();
            q.pop();
            for(const Flow& e : flows[u]) {
                if(h[e.t] == 0 && e.cap > e.flow) {
                    p[e.t] = u;
                    h[e.t] = 1;
                    a[e.t] = std::min(a[u], e.cap - e.flow);
                    q.push(e.t);
                }
            }
            if(h[end] == 1) break;
        }
        if(h[end] == 0) {
            for(int i = 0; i < fnum; i++) if(h[i] == 1) ftype[i] = 1;
            for(int i = 0; i < fnum; i++) if(ftype[i] == 3) ftype[i] = 2;
            break;
        }
        int cur = end;
        while(cur != start) {
            for(Flow& e : flows[p[cur]]) if(e.t == cur) e.flow += a[end];
            for(Flow& e : flows[cur]) if(e.t == p[cur]) e.flow -= a[end];
            cur = p[cur];
        }
    }
}
