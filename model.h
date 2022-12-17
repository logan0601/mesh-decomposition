#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include <queue>
#include <map>
#include "common.hpp"

typedef std::pair<int, int> Pair_ii;
typedef std::pair<double, int> Pair_di;

Pair_ii sort(int a, int b);


class Model
{
public:
    int vnum;                       // vertex number
    int fnum;                       // face number
    vector<Vec> vecs;               // vertex coordinate
    vector<Face> faces;             // faces
    vector<vector<double>> dmat;    // distance matrix

    double avg_dis;                 // Average Total Distance(for stopping check)

    Model();
    void read(const char* path);    // Read PLY file
    void write(const char* path);   // Write PLY file

private:
    void parse(const char* path);
    void build_edge(int i, int j, Pair_ii vid); // build edge in dual graph
    void build_graph();                         // build dual graph
    void build_short_path();                    // build distance matrix

};


#endif // MODEL_H
