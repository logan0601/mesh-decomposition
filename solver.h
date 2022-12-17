#ifndef SOLVER_H
#define SOLVER_H

#include "model.h"

typedef std::pair<vector<int>, double> Pair_vd;

class Solver
{
public:
    int lnum;                   // label number
    vector<int> labels;         // label for every face
    vector<vector<int>> lmat;   // label matrix

    Model& model;               // temporal model

    Solver(Model& model_);
    void solve();

private:
    void kernel(int l);
    Pair_vd init_seed(int l);
    bool check(int l, vector<int>& seed, double avg_dis);
    void assign(int l, vector<int>& seed);
    void assign_fuzzy(int l, vector<int>& seed, std::map<Pair_ii, vector<int>>& fuzzy);

};

#endif // SOLVER_H
