#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include <queue>
#include <map>
#include "common.h"

#define ETA 0.2
#define DELTA 0.2
#define INF 1e100
#define EPS 1e-12

typedef std::pair<int, int> pii;
typedef std::pair<double, int> pdi;

pii sort(int a, int b);


class Model
{
public:
    int vnum;                       // vertex number
    int fnum;                       // face number
    int lnum;                       // label number
    vector<Vec> vecs;               // vertex coordinate
    vector<Face> faces;             // faces
    vector<vector<double>> dmat;    // distance matrix

    double avg_dis;                 // Average Total Distance(for stopping check)
    double avg_ang_dis;             // Average Angle Distance

    Model();
    void read(const char* path);    // Read PLY file
    void write(const char* path);   // Write PLY file

    void compute_flow(vector<int>& ftype);  // assign fuzzy faces

private:
    void parse(const char* path);
    void build_edge(int i, int j, pii vid);     // build edge in dual graph
    void build_graph();                         // build dual graph
    void build_short_path();                    // build distance matrix

};


struct Flow
{
public:
    int s, t;
    double cap, flow;
    Flow(int s_, int t_, double cap_, double flow_);

};


#endif // MODEL_H
