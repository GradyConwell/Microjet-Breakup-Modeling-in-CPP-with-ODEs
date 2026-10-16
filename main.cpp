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

// Simulation parameters from Van Hoeve et al. (2010), Section III.A
// 40% w/w glycerol-saline solution (18.5 µm microjet)
struct PhysicalParameters {
    const double rho = 1098.0;   // density of 40% glycerol-saline in kg / m^3
    const double gamma = 0.0679; // surface tension in N / m
    const double mu = 0.00365;   // dynamic viscosity in Pa * s
    const double h0 = 18.5e-6;   // initial jet radius in m
    const double v0 = 5.4;       // initial jet velocity in m / s
    const double amp = 0.05;     // amplitude of initial perturbation
    const double f = 32849.33;   // optimum Rayleigh frequency (v0 / lambda_opt) in Hz
} PhysicalParameters;
// create an object to hold the physical parameters
const auto [rho, gamma, mu, h0, v0, amp, f] = PhysicalParameters;

struct SimulationParameters {
    const int Nodes = 175; // number of grid nodes
    const double Length = 700e-6; // domain length: 700 µm
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

vector<double> getInitialValues() {
    // create empty vector of initial values
    vector<double> initValues(2 * SimulationParameters.Nodes);
    // populate vector with initial conditions
    for (int i = 0; i < SimulationParameters.Nodes; i++) {
        double z = (i - 0.5) * SimulationParameters.dz;
        // modulation of "nozzle radius" to perturb the jet
        double hInit = h0 * (1 + amp * sin(2 * PI * z * f / v0));
        // maintain constant flowrate thru nozzle 
        double vInit = v0 * pow(h0, 2) / pow(hInit, 2);
        // Store in array
        initValues[2 * i] = hInit;
        initValues[2 * i + 1] = vInit;
    }
    return initValues;
}

// A detached child droplet, approximated as a point with volume, velocity and position.
struct Droplet {
    double volume;   // m^3
    double velocity; // m/s (axial)
    double position; // m (axial, centroid)
    double time;     // s (creation time)
};

// Steps 2-3 of pinch-off. Operates on the domain [0, lastNode].
// finds the neck (min radius) and the contiguous thin region around it (h <= some threshold)
// collects the neck + the fragment downstream of it (up to the next local minimum or the
// domain end) into one point-like child droplet (volume, volume-weighted velocity/position)
// closes the jet: lastNode moves to the node just upstream of the neck, and every node
// beyond it is overwritten with the new boundary values so it stays finite (it is frozen later)
// Returns false (state untouched) if the pinch-off can't be performed.
bool performPinchoff(vector<double>& state, int& lastNode, double threshold, double t, vector<Droplet>& droplets) {
    // pull in grid spacing from the global object SimulationParameters
    const double dz = SimulationParameters.dz;

    // Helper functions to return values at a given node (i)
    auto hAt = [&](int i) { return state[2 * i]; };
    auto vAt = [&](int i) { return state[2 * i + 1]; };

    // find neck = minimum radius in the interior nodes (exclude node 0 and lastNode)
    int k = 1;
    for (int i = 1; i < lastNode; i++) if (hAt(i) < hAt(k)) k = i;

    // find extent of the thin neck region (nodes included), and map into neckHi, neckLo.
    int neckLo = k, neckHi = k; 
    while (neckLo > 1 && hAt(neckLo - 1) <= threshold) neckLo--;
    while (neckHi < lastNode - 1 && hAt(neckHi + 1) <= threshold) neckHi++;

    // find the end of the fragment, and map into dropEnd.
    int dropEnd = lastNode;
    for (int i = neckHi + 1; i < lastNode; i++) {
        if (hAt(i) <= hAt(i - 1) && hAt(i) <= hAt(i + 1)) { dropEnd = i; break; }
    }

    // store the new last node that describes the jet. (1 less than the pinched regions lower bound)
    int newLast = neckLo - 1;

    // check if not enough jet left upstream (prevent negative last node)
    // shouldn't ever be invoked unless simulation was setup wrong.
    if (newLast < 2) return false;

    // redistribute neck + fragment volume into one droplet
    // create all my entries, loop through, and collect all the mass and momentum from each node.
    double volume = 0.0, momentum = 0.0, moment = 0.0;
    for (int i = neckLo; i <= dropEnd; i++) {
        double dV = PI * hAt(i) * hAt(i) * dz;
        volume += dV;
        momentum += dV * vAt(i);
        moment += dV * i * dz;
    }

    // add the droplet to the list
    droplets.push_back({volume, momentum / volume, moment / volume, t});

    // close the jet at the new boundary; freeze everything past it at finite values 
    // at the end of the live region, assign values from the closest node in the live region
    for (int i = newLast + 1; i < SimulationParameters.Nodes; i++) {
        state[2 * i] = hAt(newLast); // radius
        state[2 * i + 1] = vAt(newLast); // velocity
    }

    // update the last node in the jet domain to the new pinched off value.
    lastNode = newLast;
    return true;
}

// Builds the ODE system for the live domain [0, lastNode].
// Node 0 is held fixed, lastNode is an outflow boundary, nodes beyond lastNode are frozen.
vector<pair<expression, expression>> buildODEs(int lastNode) {
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
    for (int i = 1; i < lastNode; i++) {
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

    for (int i = 1; i < lastNode; i++) {
        C[i] = curvature(h, i);
    }
    // boundaries are "closed" to not influence the interior flow (ie, not changing at edges)
    // first node is equal to the second node
    C[0] = C[1];

    // last live node is equal to the one before it
    C[lastNode] = C[lastNode - 1];

    // Create ODEs for Heyoka
    vector<pair<expression,expression>> odes;

    // Boundary Node 0 (only boundary that I control)
    odes.push_back({h[0], expression(0.0)}); // h' = 0 (radius not changing, equal to initial radius)
    odes.push_back({v[0], expression(0.0)}); // v' = 0 (velocity not changing, equal to initial velocity)

    // Nodes
    // interior nodes only: stencils need i-1 and i+1
    for (int i = 1; i < lastNode; i++){
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

    // Boundary Node lastNode (outflow)
    // values are carried out of the domain by the flow: upwind difference with the node before it
    odes.push_back({h[lastNode], -v[lastNode] * (h[lastNode] - h[lastNode - 1]) / SimulationParameters.dz});
    odes.push_back({v[lastNode], -v[lastNode] * (v[lastNode] - v[lastNode - 1]) / SimulationParameters.dz});

    // Frozen nodes beyond the pinch-off point (removed from the live jet)
    for (int i = lastNode + 1; i < N; i++) {
        odes.push_back({h[i], expression(0.0)});
        odes.push_back({v[i], expression(0.0)});
    }

    return odes;
}

// Build the integrator (slow)
// Create integrator object
// taylor_adaptive<type>{expressions, initial values}
taylor_adaptive<double> buildIntegrator(const vector<pair<expression, expression>>& odes,
                                        const vector<double>& state, double t) {
    cout << "Building integrator" << endl;
    clock_t buildStart = clock();
    auto ta = taylor_adaptive<double>{odes, state, kw::tol = 1e-9, kw::compact_mode = true, kw::opt_level = 1u};
    ta.set_time(t);
    double buildSeconds = double(clock() - buildStart) / CLOCKS_PER_SEC;
    cout << "Build time: " << buildSeconds << " s, Taylor order: " << ta.get_order() << endl;
    return ta;
}

int main() {
    // domain is nodes [0, lastNode]; shrinks at each pinch-off
    int lastNode = N - 1;
    vector<Droplet> droplets;
    auto odes = buildODEs(lastNode);

    // get initial values for integrator
    vector<double> initValues = getInitialValues();

    auto ta = buildIntegrator(odes, initValues, 0.0);

    // Integrate
    double initialTime = 0; // initial time
    double finalTime = 70 * 1e-6; // 70 microseconds, final time
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
    file << setprecision(16);

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
        vector<double> currentState = ta.get_state();

        // Check minimum radius to detect state prior to breakup
        // only the live interior [1, lastNode - 1] is checked for the neck
        double min_h = currentState[2];
        bool has_nan = false;
        for (int i = 0; i < N; i++) {
            double hi = currentState[2 * i];
            if (std::isnan(hi)) { has_nan = true; break; }
            if (i >= 1 && i < lastNode && hi < min_h) min_h = hi;
        }

        // Check to see if the math has broken.
        if (has_nan) {
            cout << "Non-finite values at t = " << t * 1e6 << " µs, halting." << endl;
            break;
        }

        // write row and flush
        file << t;
        for (double x : currentState) file << "," << x;
        file << "\n";
        file.flush();

        // if jet has developed deep neck prior to breakup (pinch-off)
        if (min_h <= 5.0e-6) {
            cout << "Breakup: h_min = " << min_h * 1e6 << " µm at t = " << t * 1e6 << " µs." << endl;
            // perform pinchoff operations
            // 1. stop time integration (propagate_until has returned; nothing is stepping)
            // 2-3. redistribute neck volume into a child droplet (point: volume, velocity, position)
            // close the slender jet at the new boundary at the pinch-off point
            if (performPinchoff(currentState, lastNode, 5.0e-6, t, droplets)) {
                const Droplet& d = droplets.back();
                cout << "Droplet #" << droplets.size() << ": V = " << d.volume
                     << " m^3, v = " << d.velocity << " m/s, z = " << d.position * 1e6
                     << " µm. Jet now ends at node " << lastNode << endl;
                // 4. recalculate the integrator from the current state and time
                odes = buildODEs(lastNode);
                ta = buildIntegrator(odes, currentState, t);
                // 5. continue stepping (the loop continues with the new integrator)
            }
        }

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
    fifthFile << setprecision(16);
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