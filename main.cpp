#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <algorithm>

using namespace std;

// Simulation parameters
struct PhysicalParameters {
    const double rho = 1098 // density of water in kg / m^3
    const double gamma = 0.0679; // surface tension of water in N / m
    const double mu = 0.00365; // dynamic viscosity of water in Pa * s
    const double h0 = 18.5e-6; // initial jet radius in m
    const double v0 = 5.4; // initial jet velocity in m / s
    
};
struct SimulationParameters {
    const int Nodes = 120; // number of grid nodes
    const double L = 500e-6 // domain length, default of 500µm
};