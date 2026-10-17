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
    const double f = v0 / (2.0 * sqrt(2.0) * numbers::pi * h0); // optimum Rayleigh frequency (v0 / lambda_opt) in Hz
} PhysicalParameters;
// create an object to hold the physical parameters
const auto [rho, gamma, mu, h0, v0, amp, f] = PhysicalParameters;

struct SimulationParameters {
    const int Nodes = 1e3; // number of grid nodes
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
    vector<double> initValues(2 * SimulationParameters.Nodes, 0.0);
    // populate vector with initial conditions: hemispherical volume at z=0, not a pre-existing jet
    for (int i = 0; i < SimulationParameters.Nodes; i++) {
        double z = i * SimulationParameters.dz;
        if (z < h0) {
            double hInit = sqrt(max(0.0, h0 * h0 - z * z));
            double vInit = v0;
            initValues[2 * i] = hInit;
            initValues[2 * i + 1] = vInit;
        } else {
            initValues[2 * i] = 0.0;
            initValues[2 * i + 1] = 0.0;
        }
    }
    return initValues;
}

// A detached child droplet, approximated as a point with volume, velocity and position.
struct Droplet {
    double volume;   // m^3
    double velocity; // m/s (axial)
    double position; // m (axial, centroid)
    double time;     // s (creation time)

    double radius() const { return cbrt(3.0 * volume / (4.0 * PI)); }
};

// Advance each detached droplet along the axial direction and handle coalescence (Van Hoeve et al., 2010, Sec. II)
void advanceDroplets(vector<Droplet>& droplets, double dt, double currentTime = -1.0, double z_tip = -1.0) {
    for (auto& d : droplets) {
        if (currentTime < 0.0 || d.time < currentTime) {
            d.position += d.velocity * dt;
        }
    }

    const double critDist = 0.005 * h0;

    // Remove recoiling microdroplets that re-enter the continuous jet tip
    if (z_tip > 0.0) {
        erase_if(droplets, [&](const Droplet& d) {
            return d.velocity < 0.0 && (d.position - d.radius()) <= (z_tip + critDist);
        });
    }

    // Coalesce overlapping micro/satellite droplets (excluding primary-satellite merging per Fig. 5 caption)
    const double primaryRadiusCutoff = 0.8 * h0;
    // switch to size_t since droplets.size() outputs type size_t
    for (size_t i = 0; i < droplets.size(); ++i) {
        for (size_t j = i + 1; j < droplets.size(); ) {
            double r1 = droplets[i].radius(), r2 = droplets[j].radius();
            bool bothPrimary = (r1 >= primaryRadiusCutoff && r2 >= primaryRadiusCutoff);
            bool bothSub = (r1 < primaryRadiusCutoff && r2 < primaryRadiusCutoff);
            if ((bothPrimary || bothSub) && abs(droplets[i].position - droplets[j].position) <= (r1 + r2 + critDist)) {
                double V = droplets[i].volume + droplets[j].volume;
                droplets[i].velocity = (droplets[i].volume * droplets[i].velocity + droplets[j].volume * droplets[j].velocity) / V;
                droplets[i].position = (droplets[i].volume * droplets[i].position + droplets[j].volume * droplets[j].position) / V;
                droplets[i].volume = V;
                droplets[i].time = min(droplets[i].time, droplets[j].time);
                droplets.erase(droplets.begin() + j);
            } else {
                ++j;
            }
        }
    }
}

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

    // The severed fragment is the entire head downstream of the neck
    int dropEnd = lastNode;

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

    if (volume <= 0.0) return false;

    // add the droplet to the list
    droplets.push_back({volume, momentum / volume, moment / volume, t});

    // close the jet at the new boundary: smooth tip from upstream node so it does not retain singular neck
    state[2 * newLast] = state[2 * (newLast - 1)];
    state[2 * newLast + 1] = state[2 * (newLast - 1) + 1];

    for (int i = newLast + 1; i < SimulationParameters.Nodes; i++) {
        state[2 * i] = 0.00; // radius
        state[2 * i + 1] = 0.00; // velocity
    }

    // update the last node in the jet domain to the new pinched off value.
    lastNode = newLast;
    return true;
}

// Builds the ODE system for the live domain [0, lastNode].
// Node 0 is held fixed, lastNode is an outflow boundary, nodes beyond lastNode are omitted (frozen at 0).
vector<pair<expression, expression>> buildODEs(int lastNode) {
    const int numActive = lastNode + 1;
    // Create state variables only for the active nodes [0, lastNode]
    vector<expression> h(numActive), v(numActive);

    // Initialize each node as its own expression variable
    for (int i = 0; i < numActive; i++) {
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
    vector<expression> hP(numActive), hPP(numActive);
    vector<expression> vP(numActive), vPP(numActive);
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
        auto sqrtTerm1 = sqrt(term1);
        return 1.0 / (u[i] * sqrtTerm1) - hpp / (term1 * sqrtTerm1);
    };

    // calculate spatial derivative of curvature (C')
    // uses the precomputed C vector (below) so boundary values C[0], C[lastNode] are respected
    vector<expression> C(numActive);
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
    odes.reserve(2 * numActive);

    // Boundary Node 0: nozzle vibrates continuously in time
    auto dh0_dt = h0 * amp * (2.0 * PI * f) * cos(2.0 * PI * f * heyoka::time);
    auto dv0_dt = -2.0 * v0 * (h0 * h0) / (h[0] * h[0] * h[0]) * dh0_dt;
    odes.push_back({h[0], dh0_dt});
    odes.push_back({v[0], dv0_dt});

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
        dhdt += smooth * hPP[i];

        // store node in vectors
        odes.push_back({h[i], dhdt});
        odes.push_back({v[i], dvdt});
    }

    // Boundary Node lastNode (convective outflow at constant jet speed v0 so capillary recoil v < 0 remains upwind-stable)
    odes.push_back({h[lastNode], -v0 * (h[lastNode] - h[lastNode - 1]) / SimulationParameters.dz});
    odes.push_back({v[lastNode], -v0 * (v[lastNode] - v[lastNode - 1]) / SimulationParameters.dz});

    return odes;
}

// Build the integrator (slow)
// Create integrator object
// taylor_adaptive<type>{expressions, initial values}
taylor_adaptive<double> buildIntegrator(const vector<pair<expression, expression>>& odes,
                                        const vector<double>& state, double t) {
    cout << "Building integrator" << endl;
    clock_t buildStart = clock();
    vector<double> activeState(state.begin(), state.begin() + odes.size());
    auto ta = taylor_adaptive<double>{odes, std::move(activeState), kw::tol = 1e-9, kw::compact_mode = true, kw::opt_level = 0u};
    ta.set_time(t);
    double buildSeconds = double(clock() - buildStart) / CLOCKS_PER_SEC;
    cout << "Build time: " << buildSeconds << " s, Taylor order: " << ta.get_order() << endl;
    return ta;
}

int main() {
    clock_t mainStart = clock();
    // domain starts at the hemispherical droplet (z < h0) and advances as the jet flows
    int lastNode = 0;
    while ((lastNode + 1) * SimulationParameters.dz < h0 && lastNode < N - 1) {
        lastNode++;
    }
    double z_tip = lastNode * SimulationParameters.dz;
    vector<Droplet> droplets;
    auto odes = buildODEs(lastNode);
    // get initial values for integrator
    vector<double> currentState = getInitialValues();

    auto ta = buildIntegrator(odes, currentState, 0.0);

    // Integrate
    double initialTime = 0; // initial time
    double finalTime = 400 * 1e-6; // 400 microseconds, final time
    double sampleTime = 1 * 1e-6; // sample every microsecond

    // set initial time
    ta.set_time(initialTime);
    double t = ta.get_time();

    // get date
    time_t now = std::time(0);
    tm* ltm = localtime(&now);
    string date = to_string(ltm->tm_mon + 1) + "-" + to_string(ltm->tm_mday) + "-" + to_string(ltm->tm_year + 1900);

    // ---------------------------
    // Integrate until final time
    // ---------------------------

    // Check for any existing files:
    bool existing = false;
    ifstream testFile("microjet" + date + ".csv"); 
    if (testFile.is_open()) {
        existing = true;
        testFile.close();
    }

    // Create file
    string filename;    
    if (existing == false){
        filename = "microjet" + date + ".csv";
    }
    else {
        filename = "microjet" + date + "_" + to_string(clock()) + ".csv";
    }

    initCSV(filename);

    // open CSV once (state is ordered h_0, v_0, h_1, v_1, ...)
    ofstream file(filename, ios::app);
    file << setprecision(16);

    // write initial state
    file << t;
    for (double x : currentState) file << "," << x;
    file << "\n";

    // Step to each sample time and write that row right away
    // (instead of storing every adaptive step in memory)
    clock_t stepStart = clock();
    double nextSample = initialTime + sampleTime;
    while (t < finalTime) {
        ta.propagate_until(nextSample); // work until next checkpoint
        t = ta.get_time();
        const auto& activeState = ta.get_state();
        copy(activeState.begin(), activeState.end(), currentState.begin());

        // Check minimum radius to detect state prior to breakup
        // check both h and v for non-finite values across active nodes [0, lastNode]
        double min_h = currentState[2];
        bool has_nan = false;
        for (int i = 0; i <= lastNode; i++) {
            double hi = currentState[2 * i];
            double vi = currentState[2 * i + 1];
            if (!std::isfinite(hi) || !std::isfinite(vi)) { has_nan = true; break; }
            if (i >= 1 && i <= lastNode && hi < min_h) min_h = hi;
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
        const double pinchThreshold = 1.5e-6;
        if (min_h <= pinchThreshold) {
            cout << "Breakup: h_min = " << min_h * 1e6 << " µm at t = " << t * 1e6 << " µs." << endl;
            // perform pinchoff operations
            // 1. stop time integration (propagate_until has returned; nothing is stepping)
            // 2-3. redistribute neck volume into a child droplet (point: volume, velocity, position)
            // close the slender jet at the new boundary at the pinch-off point
            if (performPinchoff(currentState, lastNode, pinchThreshold, t, droplets)) {
                z_tip = lastNode * SimulationParameters.dz;
                const Droplet& d = droplets.back();
                cout << "Droplet #" << droplets.size() << ": V = " << d.volume
                     << " m^3, v = " << d.velocity << " m/s, z = " << d.position * 1e6
                     << " µm. Jet now ends at node " << lastNode << endl;
                // 4. recalculate the integrator from the current state and time
                odes = buildODEs(lastNode);
                ta = buildIntegrator(odes, currentState, t);
                // 5. continue stepping (the loop continues with the new integrator)
            }
        } else if (lastNode < N - 1) {
            // Advance jet front if not yet at the domain outflow
            double v_tip = currentState[2 * lastNode + 1];
            z_tip += v_tip * sampleTime;
            bool advanced = false;
            while (z_tip >= (lastNode + 1) * SimulationParameters.dz && lastNode < N - 1) {
                lastNode++;
                currentState[2 * lastNode] = currentState[2 * (lastNode - 1)];
                currentState[2 * lastNode + 1] = currentState[2 * (lastNode - 1) + 1];
                advanced = true;
            }
            if (advanced) {
                odes = buildODEs(lastNode);
                ta = buildIntegrator(odes, currentState, t);
            }
        }

        // Advance each droplet based upon its current velocity and check for coalescence
        advanceDroplets(droplets, sampleTime, t, z_tip);

        nextSample += sampleTime;
    }
    clock_t stepEnd = clock();
    double stepSeconds = double(stepEnd - stepStart) / CLOCKS_PER_SEC;
    cout << "Total integration time: " << stepSeconds << " s" << endl;
    cout << "Total elapsed time: " << double(clock() - mainStart) / CLOCKS_PER_SEC << " s" << endl;
    return 0;
}