#include "model.h"
#include "solver.h"


int main()
{
    Model model;
    model.read("./data/dino.ply");

    Solver solver(model);
    solver.solve();

    model.write("./1.ply");
    return 0;
}
