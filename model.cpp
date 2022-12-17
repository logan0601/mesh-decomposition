#include "model.h"

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <omp.h>
#include <queue>


Pair_ii sort(int a, int b)
{
    return (a < b) ? Pair_ii(a, b) : Pair_ii(b, a);
}


Model::Model()
{
    vnum = 0;
    fnum = 0;
    lnum = 0;
    dmat = vector<vector<double>>( Utils::MAX_FACES, vector<double>(Utils::MAX_FACES, 0) );
}


void Model::read(const char* path)
{
    double start_time = omp_get_wtime();

    parse(path);

    build_graph();

    build_short_path();

    double end_time = omp_get_wtime();
    printf("Using %.2lf seconds\n", end_time - start_time);
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

    for(int i = 0; i < fnum; i++)
    {
        const Face& f = faces[i];
        int r = 60 * (labels[i] % 4 + 1);
        int g = 80 * ((labels[i] + 1) % 3 + 1);
        int b = 50 * ((labels[i] + 2) % 5 + 1);
        os << "3 " << f.vid[0] << " " << f.vid[1] << " " << f.vid[2] << " ";
        os << r << " " << g << " " << b << "\n";
    }
    
    os.close();
}


void Model::parse(const char* path)
{
    std::ifstream in;
    in.open(path, std::fstream::in);

    std::string line;
    while(getline(in, line))
    {
        if(line.find("element vertex") == 0) vnum = std::stoi(line.substr(15));
        if(line.find("element face") == 0)   fnum = std::stoi(line.substr(13));
        if(line.find("end_header") == 0) break;
    }

    double x[3];
    int p[4] = {0}, k[3];
    for(int i = 0; i < vnum; i++)
    {
        getline(in, line);
        for(int j = 1; j < 3; j++) p[j] = line.find(" ", p[j-1]) + 1;
        for(int j = 0; j < 3; j++) x[j] = std::stod(line.substr(p[j]));

        vecs.push_back(Vec(x[0], x[1], x[2]));
    }

    for(int i = 0; i < fnum; i++)
    {
        getline(in, line);
        for(int j = 1; j < 4; j++) p[j] = line.find(" ", p[j-1]) + 1;
        for(int j = 0; j < 3; j++) k[j] = std::stoi(line.substr(p[j+1]));

        faces.push_back(Face(k[0], k[1], k[2], vecs[k[0]], vecs[k[1]], vecs[k[2]]));
    }

    in.close();
}


void Model::build_edge(int i, int j, Pair_ii vid)
{
    Face& f1 = faces[i];
    Face& f2 = faces[j];

    double angle = acos( Vec::dot( f1.normal, f2.normal ) );
    bool is_convex = Vec::dot( f1.normal, f2.center-f1.center ) < Utils::EPS;
    double eta = is_convex ? Utils::ETA : 1;
    double ang_dis = eta * ( 1 - Vec::dot( f1.normal, f2.normal ) );

    Vec ed = vecs[vid.second] - vecs[vid.first];
    Vec v1 = f1.center - vecs[vid.first];
    Vec v2 = f2.center - vecs[vid.first];
    double elen = ed.len();
    double len1 = v1.len();
    double len2 = v2.len();
    double angle1 = acos( Vec::dot( ed, v1 ) / elen / len1 );
    double angle2 = acos( Vec::dot( ed, v2 ) / elen / len2 );
    double geo_dis = sqrt( v1.len2() + v2.len2() - 2 * len1 * len2 * cos( angle1 + angle2 ) );

    f1.nbrs.push_back(Edge(vid.first, vid.second, j, angle, ang_dis, geo_dis));
    f2.nbrs.push_back(Edge(vid.first, vid.second, i, angle, ang_dis, geo_dis));
}


void Model::build_graph()
{
    int emap[3][2] = { {0,1}, {0,2}, {1,2} };
    std::map<Pair_ii, int> visited;
    for(int i = 0; i < fnum; i++)
    {
        const Face& f = faces[i];
        for(int j = 0; j < 3; j++)
        {
            Pair_ii p = sort(f.vid[emap[j][0]], f.vid[emap[j][1]]);
            if(visited.find(p) != visited.end())
            {
                build_edge(visited[p], i, p);
            }
            else visited[p] = i;
        }
    }

    int count = 0;
    double ang_dis = 0, geo_dis = 0;
    for(const Face& f : faces)
    {
        count += f.nbrs.size();
        for(const Edge& e : f.nbrs)
        {
            ang_dis += e.ang_dis;
            geo_dis += e.geo_dis;
        }
    }

    double avg_ang_dis = ang_dis / count;
    double avg_geo_dis = geo_dis / count;
    for(Face& f : faces)
    {
        for(Edge& e : f.nbrs)
            e.tot_dis = Utils::DELTA * e.ang_dis / avg_ang_dis
                + ( 1 - Utils::DELTA ) * e.geo_dis / avg_geo_dis;
    }
}


void Model::build_short_path()
{
    vector<int> visited;
    vector<double> distant;
    std::priority_queue<Pair_di, vector<Pair_di>, std::greater<Pair_di>> q;

    #pragma omp parallel for private(visited) private(distant) private(q)
    for(int i = 0; i < fnum; i++)
    {
        visited = vector<int>( fnum, 0 );
        distant = vector<double>( fnum, Utils::UPPER );
        q = std::priority_queue<Pair_di, vector<Pair_di>, std::greater<Pair_di>>();

        distant[i] = 0.0;
        q.push(Pair_di(0.0, i));
        while(!q.empty())
        {
            int u = q.top().second;
            q.pop();

            if(visited[u]) continue;
            visited[u] = 1;

            for(const Edge& e : faces[u].nbrs)
            {
                int v = e.fid;
                double w = e.tot_dis;
                if(distant[u] + w < distant[v])
                {
                    distant[v] = distant[u] + w;
                    q.push(Pair_di(distant[v], v));
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
