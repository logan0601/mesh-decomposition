#ifndef SOLVER_H
#define SOLVER_H

#include <algorithm>
#include "model.h"

#define NREP 4
#define LEVEL 1

using std::distance;
using std::min_element;
using std::max_element;


class Solver
{
public:
    Model& model;                   // temporal model
    int fnum;                       // face number
    int rnum;                       // repre number
    int level;                      // iterate level
    int lab_idx, fuz_idx;           // index

    vector<int> fid;                // faces id
    vector<int> rid;                // repre id
    vector<int> unique;             // unique repre id
    vector<vector<double>> dmat;    // distance matrix
    vector<vector<double>> pmat;    // probability matrix
    vector<vector<double>> cmat;    // cost matrix

    double avg_dis;                 // average distance
    int old_lab;                    // old label for print

    Solver(Model& mod_, int lev_, vector<int> fid_=vector<int>());
    void solve();

private:
    void initial_repre();
    void compute_prob_matrix();
    void assign();
    void assign_fuzzy();
    vector<int> recompute_repre();
    bool check();
};


#endif // SOLVER_H
