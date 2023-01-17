#include "model.h"
#include "solver.h"


int main()
{
    Model model;
    model.read("./data/dino.ply");

    Solver solver(model, 0);
    solver.solve();

    model.write("./out/dino.ply");
    return 0;
}
