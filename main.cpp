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

    // stencils
    auto d1 = [&](vector<expression>& u, int i) {
        return (u[i + 1] - u[i - 1]) / (2.0 * SimulationParameters.dz);
    };
    // create second derivative stencil
    auto d2 = [&](vector<expression>& u, int i) {
        return (u[i + 1] - 2.0 * u[i] + u[i - 1]) / (SimulationParameters.dz * SimulationParameters.dz);
    };

    // calculate curvature at a node (C)
    // pass in curvature(expression, index)
    auto curvature = [&](vector<expression>& u, int i) {
        auto hp = d1(u, i); // h prime at i
        auto hpp = d2(u, i); // h double prime at i
        auto term1 = 1.0 + hp * hp; // (1 + (h')^2) term    
        return 1.0 / (u[i] * sqrt(term1)) - hpp / pow(term1, 1.5);
    };

    // calculate spatial derivative of curvature (C')
    auto dC = [&](vector<expression>& u, int i) {
        auto Cp = (curvature(u, i + 1) - curvature(u, i - 1)) / (2.0 * SimulationParameters.dz);
        return Cp;
    };
    
    // create curvature vector for all points i (excluding boundaries)
    vector<expression> C(SimulationParameters.Nodes);   
    C.reserve(SimulationParameters.Nodes);

    for (int i = 1; i < SimulationParameters.Nodes - 1; i++) {
        C[i] = curvature(h, i);
    }
    // boundaries are "closed" to not influence the interior flow (ie, not changing at edges)
    // first node is equal to the second node
    C[0] = C[1];
    // last node is equal to the second to last node
    C[SimulationParameters.Nodes - 1] = C[SimulationParameters.Nodes - 2];

    // Create ODEs for Heyoka
    vector<pair<expression,expression>> odes;

    // Boundary Node 0 (only boundary that I control)
    odes.push_back({h[0], 0.0}); // h' = 0 (radius not changing, equal to initial radius)
    odes.push_back({v[0], 0.0}); // v' = 0 (velocity not changing, equal to initial velocity)

    // Nodes
    for (int i = 0, i < N-1, i++){
        // set up terms of equations.
        // patterned in stencils above.
        auto hp = d1(h, i); // h prime
        auto vp = d1(v, i); // v prime
        auto vpp = d2(v, i);// v double prime
        auto Cp = dC(h, i); // C prime

        auto dvdt = -v[i] * vp - (gamma / rho) * Cp + 3 * mu / rho * (vpp + 2 * hp / h[i] * vp);
        auto dhdt =

        
        
    }

    
return 0;
}