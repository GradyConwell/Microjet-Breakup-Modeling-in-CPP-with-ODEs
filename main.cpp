#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <ctime>
#include <numbers>
#include <heyoka/heyoka.hpp>

using namespace std;
using namespace heyoka;

// Simulation parameters, mostly pulled from Van houve et al. III.a
// "A 40% w/w glycerol-saline solution"
struct PhysicalParameters {
    // const double rho = 1098; // density of water in kg / m^3
    // const double gamma = 0.0679; // surface tension of water in N / m
    // const double mu = 0.00365; // dynamic viscosity of water in Pa * s

    // 40% glycerol saline solution parameters
    const double rho = 1003; // density of saline solution in kg / m^3
    const double gamma = 0.072; // surface tension of saline solution in N / m
    const double mu = 0.001; // dynamic viscosity of saline solution in Pa * s

    const double h0 = 1.25e-6; // initial jet radius in m
    const double v0 = 35; // initial jet velocity in m / s
    const double amp = 0.1; // amplitude of initial perturbation (very small)
    const double f = 33e3; // perturbation frequency Hz
} PhysicalParameters;
const auto [rho, gamma, mu, h0, v0, amp, f] = PhysicalParameters;

struct SimulationParameters {
    const int Nodes = 120; // number of grid nodes
    const double Length = 500e-6; // domain length, default of 500µm
    const double dz = Length / (Nodes - 1); // grid spacing
} SimulationParameters;
const int N = SimulationParameters.Nodes;
constexpr double PI = numbers::pi;

void initCSV(const string& filename) {
    // open file
    ofstream file(filename);
    // write header
    file << "t";
    for (int i = 0; i < N; i++) file << ",h_" << i << ",v_" << i;
    file << endl;
    // close file
    file.close();
}

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
    // 1OCD
    auto d1 = [&](vector<expression>& u, int i) {
        return (u[i + 1] - u[i - 1]) / (2.0 * SimulationParameters.dz);
    };
    // create second derivative stencil
    // 2OCD
    auto d2 = [&](vector<expression>& u, int i) {
        return (u[i + 1] - 2.0 * u[i] + u[i - 1]) / (SimulationParameters.dz * SimulationParameters.dz);
    };

    // smoothing coefficient (m^2/s), scaled to grid spacing and jet velocity
    const double smooth = 0.1 * SimulationParameters.dz * v0;

    // Build the stencils once per node and store them,
    // so curvature and the ODEs below reuse the same expressions.
    vector<expression> hP(SimulationParameters.Nodes), hPP(SimulationParameters.Nodes);
    vector<expression> vP(SimulationParameters.Nodes), vPP(SimulationParameters.Nodes);
    for (int i = 1; i < SimulationParameters.Nodes - 1; i++) {
        hP[i] = d1(h, i);
        hPP[i] = d2(h, i);
        vP[i] = d1(v, i);
        vPP[i] = d2(v, i);
    }

    // calculate curvature at a node (C)
    // pass in curvature(expression, index)
    auto curvature = [&](vector<expression>& u, int i) {
        auto hp = hP[i]; // h prime at i
        auto hpp = hPP[i]; // h double prime at i
        auto term1 = 1.0 + hp * hp; // (1 + (h')^2) term    
        return 1.0 / (u[i] * sqrt(term1)) - hpp / pow(term1, 1.5);
    };

    // calculate spatial derivative of curvature (C')
    // uses the precomputed C vector (below) so boundary values C[0], C[N-1] are respected
    vector<expression> C(SimulationParameters.Nodes);
    auto dC = [&](int i) {
        auto Cp = (C[i + 1] - C[i - 1]) / (2.0 * SimulationParameters.dz);
        return Cp;
    };
    
    // create curvature vector for all points i (excluding boundaries)

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
    odes.push_back({h[0], expression(0.0)}); // h' = 0 (radius not changing, equal to initial radius)
    odes.push_back({v[0], expression(0.0)}); // v' = 0 (velocity not changing, equal to initial velocity)

    // Nodes
    // interior nodes only: stencils need i-1 and i+1
    for (int i = 1; i < N-1; i++){
        // set up terms of equations.
        // patterned in stencils above.
        auto hp = hP[i]; // h prime
        auto vp = vP[i]; // v prime
        auto vpp = vPP[i];// v double prime
        auto Cp = dC(i); // C prime

        // Two equations from paper:
        auto dvdt = -v[i] * vp - (gamma / rho) * Cp + 3 * mu / rho * (vpp + 2.0 * hp / h[i] * vp);
        auto dhdt = -v[i] * hp - 0.500 * h[i] * vp;

        // smoothing. Just pave over the sawtooth behavior by running it over with the second derivative.
        dvdt += smooth * vpp;
        dhdt += smooth * d2(h, i);

        // store node in vectors
        odes.push_back({h[i], dhdt});
        odes.push_back({v[i], dvdt});
    }

    // Boundary Node N-1 (outflow)
    // values are carried out of the domain by the flow: upwind difference with the node before it
    odes.push_back({h[N-1], -v[N-1] * (h[N-1] - h[N-2]) / SimulationParameters.dz});
    odes.push_back({v[N-1], -v[N-1] * (v[N-1] - v[N-2]) / SimulationParameters.dz});

    // Create integrator object
    // taylor_adaptive<type>{expressions, initial values}
    vector<double> initValues(2 * SimulationParameters.Nodes);
    for (int i = 0; i < SimulationParameters.Nodes; i++) {
        double z = i * SimulationParameters.dz;
        // modulation of "nozzle radius" to perturb the jet
        double hInit = h0*(1 + amp*sin(2 * PI * z * f / v0));
        // maintain constant flowrate thru nozzle 
        double vInit = v0 * pow(h0, 2) / pow(hInit, 2);
        // Store in array
        initValues[2*i] = hInit;
        initValues[2*i + 1] = vInit;
    }
   
    // Build the integrator (slow)
    // type = taylor_adaptive
    // tol = 1e-9
    // compact_mode = true
    // opt_level = 1
    cout << "Building integrator" << endl;
    clock_t buildStart = clock();
    auto ta = taylor_adaptive<double>{odes, initValues, kw::tol = 1e-12, kw::compact_mode = true, kw::opt_level = 1u};
    clock_t buildEnd = clock();
    double buildSeconds = double(buildEnd - buildStart) / CLOCKS_PER_SEC;
    cout << "Build time: " << buildSeconds << " s, Taylor order: " << ta.get_order() << endl;

    // Integrate
    double initialTime = 0; // initial time
    double finalTime = 50 * 1e-6; // 50 microseconds, final time
    double sampleTime = 1 * 1e-6; // sample every microsecond

    // set initial time
    ta.set_time(initialTime);
    double t = ta.get_time();

    // get date
    time_t now = std::time(0);
    tm* ltm = localtime(&now);
    string date = to_string(ltm->tm_mon + 1) + "-" + to_string(ltm->tm_mday) + "-" + to_string(ltm->tm_year + 1900);

    // Integrate until final time
    string filename = "microjet" + date + ".csv";
    initCSV(filename);

    // open CSV once (state is ordered h_0, v_0, h_1, v_1, ...)
    ofstream file(filename, ios::app);
    file << setprecision(17);

    // write initial state
    file << t;
    for (double x : ta.get_state()) file << "," << x;
    file << "\n";

    // Step to each sample time and write that row right away
    // (instead of storing every adaptive step in memory)
    clock_t stepStart = clock();
    double nextSample = initialTime + sampleTime;
    while (t < finalTime) {
        ta.propagate_until(nextSample); // work until next checkpoint
        t = ta.get_time();

        // write row and flush
        file << t;
        for (double x : ta.get_state()) file << "," << x;
        file << "\n";
        file.flush();

        nextSample += sampleTime;
    }
    clock_t stepEnd = clock();
    double stepSeconds = double(stepEnd - stepStart) / CLOCKS_PER_SEC;
    cout << "Total integration time: " << stepSeconds << " s" << endl;

    // write the last frame in "long" format
    // columns: node, z, h, v   (z is implied by the node index: z = node * dz)
    ofstream lastFile("microjet_last_frame.csv");
    lastFile << setprecision(17);
    lastFile << "node,z,h,v" << endl;
    vector<double> lastState = ta.get_state();
    for (int i = 0; i < N; i++) {
        lastFile << i << "," << i * SimulationParameters.dz << "," << lastState[2*i] << "," << lastState[2*i + 1] << "\n";
    }

    file.close();

    // write the fifth frame in long format from the microjetDATE.csv file
    // read microjetDATE.csv, find the row at t = 5 * sampleTime (5.0e-6 s), and write it to a new CSV
    ifstream inFile(filename);
    ofstream fifthFile("microjet_5th_frame.csv");
    fifthFile << setprecision(17);
    fifthFile << "node,z,h,v\n";

    double targetTime = 5 * sampleTime; // 5 µs
    string line;
    // Skip header line
    if (getline(inFile, line)) {
        while (getline(inFile, line)) {
            stringstream ss(line);
            string cell;
            if (!getline(ss, cell, ',')) continue;
            
            double rowTime = stod(cell);
            if (abs(rowTime - targetTime) < 0.5 * sampleTime) {
                // Parse the state values: h_0, v_0, h_1, v_1, ...
                vector<double> stateValues;
                while (getline(ss, cell, ',')) {
                    stateValues.push_back(stod(cell));
                }
                for (int i = 0; i < N && (2 * i + 1) < static_cast<int>(stateValues.size()); i++) {
                    fifthFile << i << "," << i * SimulationParameters.dz << ","
                              << stateValues[2 * i] << "," << stateValues[2 * i + 1] << "\n";
                }
                break;
            }
        }
    }
    fifthFile.close();
    inFile.close();

    return 0;
}