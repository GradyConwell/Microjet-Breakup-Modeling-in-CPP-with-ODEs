#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <heyoka/heyoka.hpp>

using namespace std;

// Simulation parameters, mostly pulled from Van houve et al. III.a
// "A 40% w/w glycerol-saline solution"
struct PhysicalParameters {
    const double rho = 1098 // density of water in kg / m^3
    const double gamma = 0.0679; // surface tension of water in N / m
    const double mu = 0.00365; // dynamic viscosity of water in Pa * s
    const double h0 = 18.5e-6; // initial jet radius in m
    const double v0 = 5.4; // initial jet velocity in m / s
    
};
struct SimulationParameters {
    const int Nodes = 120; // number of grid nodes
    const double Length = 500e-6; // domain length, default of 500µm
    const double dz = Length / (Nodes - 1); // grid spacing
};

int main() {
    // Create state variables for each nodes in the domain
    // let h be the jet radius at each node, and v be the axial velocity at each node
    vector<expression> h(SimulationParameters.Nodes), v(SimulationParameters.Nodes);
    // explicitly allocate space for each node
    h.reserve(SimulationParameters.Nodes);
    v.reserve(SimulationParameters.Nodes);

    // Initialize each node as its own expression variable
    for (int i = 0; i < SimulationParameters.Nodes; i++) {
        h[i] = make_vars("h_" + to_string(i));
        v[i] = make_vars("v_" + to_string(i));
    }

    // https://en.wikipedia.org/wiki/Stencil_(numerical_analysis)
    // Next step: create expressions for stencils IOT create set of solvable ODEs for heyoka

    return 0;
}