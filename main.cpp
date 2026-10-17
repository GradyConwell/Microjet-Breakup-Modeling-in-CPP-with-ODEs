// Grady Conwell, 2026
// 
// Rayleigh breakup of a 18.5 µm microjet, after Van Hoeve et al. (2010), "Breakup of diminutive Rayleigh jets",
// using the slender-jet (lubrication) equations of Eggers & Dupont, J. Fluid Mech. 262, 205 (1994).
//
// lengths in h0, time in the capillary time tc = sqrt(rho h0^3 / gamma).
// The jet is written in terms of a = h^2, which keeps the equations regular at a rounded tip (h -> 0):
//   a_t + (a v)'  = 0
//   (a v)_t + (a v^2)' = G' + 3 Oh (a v')',   G = h^2 K = a^2 (1 + a''/2) / (a + a'^2/4)^(3/2)
// G is the axial surface-tension force on a cross-section (divided by pi gamma h0), so both equations
// conserve mass and momentum exactly.
//
// N cells spread evenly from the nozzle (z = 0) to the jet tip (z = L), so z = xi L(t).
// Cell centres hold a, cell faces hold v (staggered), and the tip moves with the fluid: dL/dt = v_tip.
// Pinch-off: when the neck radius drops below 100 nm, the part downstream becomes a droplet
// and the remaining jet is redistributed over the N cells. A free droplet's centre of mass moves at constant
// velocity (no gravity or air interaction is modeled), so detached drops are carried as (volume, momentum). Objects that touch
// (gap < 0.005 h0) while approaching coalesce. A drop caught by the jet tip is absorbed into the tip cells.
//
// VARIABLE DICTIONARY   [dimensionless] unless SI units are given. Scales: length h0, time tc, velocity uc.
// ---------------------------------------------------------------------------------------------------------
// Physical inputs (SI)
//   rho        Liquid density [kg/m^3]
//   gam        Surface tension gamma [N/m]
//   mu         Dynamic viscosity [Pa s]
//   h0         Unperturbed jet radius [m]; the length scale
//   vjet       Jet velocity at the nozzle [m/s]
//   delta      Nozzle radius modulation amplitude, as a fraction of h0 (paper: 0.005)
//
// Scales and dimensionless groups
//   tc         Capillary time sqrt(rho h0^3 / gam) [s]; the time scale (~10 µs)
//   uc         Velocity scale h0 / tc [m/s]
//   Oh         Ohnesorge number mu / sqrt(rho gam h0): viscosity vs. inertia and surface tension (~0.10)
//   U          Dimensionless jet velocity vjet / uc = sqrt(We)
//   omega      Angular nozzle forcing frequency 2 pi f, with f = vjet / lambda_opt
//   PI         pi
//
// Numerical settings
//   N          Number of grid cells from the nozzle to the jet tip
//   hPinch     Pinch-off threshold radius (100 nm / h0)
//   gapMerge   Gap below which approaching objects coalesce (0.005 h0)
//   dtOut      Time between output frames (1 µs)
//   tEnd       End time of the simulation (set from argv[1] in µs)
//
// Grid and state (z = xi L, xi in [0, 1]; cell j spans faces j and j + 1)
//   s, x       State vector: x are its heyoka symbols, s its numerical values (s[i] is the value of x[i])
//                s[0 .. N-1]     a_j: squared radius h^2 at cell centres
//                s[iV .. iV+N-1] v_1 .. v_N: axial velocity at cell faces (v_N is the tip velocity)
//                s[iL]           L: jet length, nozzle to tip
//   iV, iL     Offsets of the velocities and of L in the state vector
//   a, a(j)    Squared jet radius h^2 (a(-1) and a(N) are ghost cells for the nozzle and tip)
//   h          Jet radius sqrt(a)
//   v, v(k)    Axial velocity at face k (v(0) is the nozzle inflow velocity vn)
//   L          Jet length (position of the tip)
//   Ldot       dL/dt, the tip velocity (= v_N)
//   dz         Cell width L / N
//   j, k       Cell index, face index
//
// Nozzle boundary
//   hn         Nozzle radius 1 + delta sin(omega t)
//   an         Nozzle area term hn^2
//   vn         Nozzle velocity U / an (keeps the flow rate constant, paper Eq. 9)
//   nozzleA()  an evaluated numerically at a given time
//
// Equation terms (buildJet)
//   ap, app    a' and a'' (first and second z-derivatives of a)
//   d          a + a'^2 / 4, the denominator term of the curvature
//   G[j]       Capillary tension h^2 K at cell j (axial surface-tension force / (pi gam h0))
//   S[j]       Viscous tension a v' at cell j (multiplied by 3 Oh in the momentum equation)
//   F[k]       Mass flux through face k, relative to the moving grid
//   adv        Momentum advection term (v - grid velocity) v'
//   force      Capillary + viscous force per unit mass at a face
//   pinch      Terminal events that stop the integrator when any a_j falls to hPinch^2
//   odes       (variable, d/dt) pairs passed to heyoka
//   ta         heyoka Taylor integrator for the jet
//
// Drops (struct Drop)
//   V          Drop volume
//   P          Drop momentum V v
//   z          Drop centre of mass at time t
//   t          Reference time for z
//   v()        Drop velocity P / V (constant: free drops feel no net force)
//   R()        Equivalent sphere radius (3V / 4pi)^(1/3)
//   at(time)   Centre of mass at the given time
//   drops      All detached drops in flight
//   p, q       Neighbouring drops (p behind q) checked for coalescence
//
// Pinch-off, absorption and remapping
//   findNeck() Index of the thinnest local minimum below hMax that has a real bulge downstream, else -1
//   hMax       Radius below which findNeck() accepts a neck
//   bulge      Largest a seen downstream of the cell being tested
//   best       Index of the thinnest neck found so far
//   neck       Cell index of the current pinch-off neck (-1 if none)
//   caught     Index in drops of a drop the jet tip has caught (-1 if none)
//   pinches    Number of pinch-offs so far (event index in neck.csv)
//   remap()    Redistributes the N cells over a new jet length by linear interpolation
//   Lnew       New length of the slender jet after a pinch-off
//   aAt, vAt   Interpolate a (cell centres) and v (faces) of the old grid at a position z
//   p          Fractional grid index used for interpolation
//   dV         Volume of one cell added to a new drop
//   n          Number of tip cells that receive an absorbed drop
//   dA         Increase in a for each of those cells
//   aOld       a of a tip cell before absorption (to weight the face velocity)
//
// Output and bookkeeping
//   um, us     Conversion factors to µm and µs for output
//   tOut       Time of the next output frame
//   onStep     Callback after every integrator step: logs the neck radius and checks for drops the tip has caught
//   dir        Output folder results/<YYYYMMDD_HHMMSS>
//   jetCsv     jet.csv: t, L and h at every cell centre, one row per frame
//   dropCsv    drops.csv: t, z, R, v of every drop in flight, per frame
//   eventCsv   events.csv: every pinch-off and absorption
//   neckCsv    neck.csv: minimum neck radius after every step (for the Fig. 6 power law)
// ---------------------------------------------------------------------------------------------------------

#include <algorithm>
#include <chrono>
#include <ctime>
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <numbers>
#include <vector>
#include <heyoka/heyoka.hpp>

using namespace std;
using namespace heyoka;
constexpr double PI = numbers::pi;

// 40% w/w glycerol-saline microjet, Van Hoeve et al. (2010) Sec. III.A
const double rho = 1098.0, gam = 0.0679, mu = 3.65e-3; // kg/m^3, N/m, Pa s
const double h0 = 18.5e-6, vjet = 5.4;                 // m, m/s
const double delta = 0.005;                            // nozzle modulation amplitude, delta / h0

// Scales and dimensionless groups
const double tc = sqrt(rho * h0 * h0 * h0 / gam); // capillary time, ~10 µs
const double uc = h0 / tc;                        // velocity scale
const double Oh = mu / sqrt(rho * gam * h0);      // ~0.10
const double U = vjet / uc;                       // sqrt(We) ~ 2.95
const double omega = U / sqrt(2.0);               // 2 pi f with f = vjet / lambda_opt, lambda_opt = 2 sqrt(2) pi h0

// Simulation parameters
const int N = 1000;                       // grid points, nozzle to tip (as in the paper). Changed from a global system with nodes that are fixed in space.
const double hPinch = 100e-9 / h0;        // pinch-off threshold, 100 nm. (much lower than before, singularity issues addressed)
const double gapMerge = 0.005;            // coalescence distance, 0.005 h0
const double dtOut = 1e-6 / tc;           // 1 µs frames (Data output resolution)
double tEnd = 400e-6 / tc;                // override with ./microjet <t_end in µs>

struct Drop {
    double V, P, z, t; // volume, momentum (V v), centre at time t
    double v() const { return P / V; }
    double R() const { return cbrt(3.0 * V / (4.0 * PI)); }
    double at(double time) const { return z + v() * (time - t); }
};

// State layout: a_0..a_{N-1} (cells), v_1..v_N (faces, v_N = tip), L.
const int iV = N, iL = 2 * N;
double nozzleA(double t) { double hn = 1.0 + delta * sin(omega * t); return hn * hn; }

taylor_adaptive<double> buildJet(vector<double> init) {
    vector<expression> x(2 * N + 1);
    for (int i = 0; i <= 2 * N; i++) x[i] = make_vars("x" + to_string(i));

    // Nozzle: h = h0 (1 + delta sin wt) with constant flow rate, v = vjet h0^2 / h^2 (paper, Eqs. 8-9)
    auto hn = 1.0 + delta * sin(omega * heyoka::time);
    auto an = hn * hn, vn = U / an;
    auto L = x[iL], Ldot = x[iV + N - 1], dz = L / double(N);
    // Ghost cells pin a = a_nozzle at z = 0 and a = 0 at the tip
    auto a = [&](int j) { return j < 0 ? 2.0 * an - x[0] : j >= N ? -x[N - 1] : x[j]; };
    auto v = [&](int k) { return k == 0 ? vn : x[iV + k - 1]; };

    vector<expression> G(N), S(N), F(N + 1);
    for (int j = 0; j < N; j++) {
        auto ap = (a(j + 1) - a(j - 1)) / (2.0 * dz), app = (a(j + 1) - 2.0 * a(j) + a(j - 1)) / (dz * dz);
        auto d = a(j) + 0.25 * ap * ap;
        G[j] = a(j) * a(j) * (1.0 + 0.5 * app) / (d * sqrt(d)); // capillary tension
        S[j] = a(j) * (v(j + 1) - v(j)) / dz;                     // viscous tension / 3 Oh
    }
    F[0] = an * vn; // inflow
    F[N] = expression(0.0); // the tip moves with the fluid
    for (int k = 1; k < N; k++) F[k] = 0.5 * (a(k - 1) + a(k)) * (v(k) - (double(k) / N) * Ldot);

    vector<pair<expression, expression>> odes;
    for (int j = 0; j < N; j++) odes.push_back({x[j], -(F[j + 1] - F[j]) / dz - a(j) * Ldot / L});
    for (int k = 1; k < N; k++) {
        auto adv = (v(k) - (double(k) / N) * Ldot) * (v(k + 1) - v(k - 1)) / (2.0 * dz);
        auto force = (G[k] - G[k - 1] + 3.0 * Oh * (S[k] - S[k - 1])) / (0.5 * (a(k - 1) + a(k)) * dz);
        odes.push_back({x[iV + k - 1], force - adv});
    }
    // Tip: the last half cell is pulled back by the tension of its neighbour (zero area, zero force at the tip)
    odes.push_back({x[iV + N - 1], -(G[N - 1] + 3.0 * Oh * S[N - 1]) / (0.5 * a(N - 1) * dz)});
    odes.push_back({L, Ldot});

    // Terminal events stop the integrator exactly when a cell reaches the pinch-off radius
    vector<taylor_adaptive<double>::t_event_t> pinch;
    for (int j = 3; j < N - 3; j++) pinch.emplace_back(x[j] - hPinch * hPinch, kw::direction = event_direction::negative);
    return taylor_adaptive<double>{odes, std::move(init), kw::tol = 1e-10, kw::compact_mode = true, kw::t_events = pinch};
}

// Redistribute the grid over the jet [0, Lnew] (linear interpolation).
void remap(double* s, double t, double Lnew) {
    const double L = s[iL], dz = L / N;
    vector<double> a(s, s + N), v(N + 1);
    v[0] = U / nozzleA(t);
    copy(s + iV, s + iV + N, v.begin() + 1);
    auto aAt = [&](double z) { // cell centres, a = 0 at the tip
        double p = z / dz - 0.5;
        if (p <= 0) return a[0];
        if (p >= N - 1) return max(0.0, a[N - 1] * (L - z) / (0.5 * dz));
        int j = int(p); return a[j] + (p - j) * (a[j + 1] - a[j]);
    };
    auto vAt = [&](double z) { // faces
        double p = min(z, L) / dz; int k = min(int(p), N - 1);
        return v[k] + (p - k) * (v[k + 1] - v[k]);
    };
    for (int j = 0; j < N; j++) s[j] = aAt((j + 0.5) * Lnew / N);
    for (int k = 1; k <= N; k++) s[iV + k - 1] = vAt(k * Lnew / N);
    s[iL] = Lnew;
}

// Neck: a local minimum of h below the threshold with a real bulge (drop) downstream of it.
int findNeck(const double* s, double hMax) {
    int best = -1;
    double bulge = 0.0;
    for (int j = N - 1; j > 2; j--) {
        bulge = max(bulge, s[j]);
        if (s[j] < hMax * hMax && s[j] <= s[j - 1] && s[j] <= s[j + 1] && j < N - 3 && bulge > 100.0 * hPinch * hPinch
            && (best < 0 || s[j] < s[best]))
            best = j;
    }
    return best;
}

int main(int argc, char** argv) {
    cout << unitbuf;
    if (argc > 1) tEnd = stod(argv[1]) * 1e-6 / tc;
    char stamp[32];
    time_t now = std::time(nullptr);
    strftime(stamp, sizeof stamp, "%Y%m%d_%H%M%S", localtime(&now));
    auto dir = filesystem::path("results") / stamp;
    filesystem::create_directories(dir);
    ofstream jetCsv(dir / "jet.csv"), dropCsv(dir / "drops.csv"), eventCsv(dir / "events.csv"), neckCsv(dir / "neck.csv");
    jetCsv << "t_us,L_um";
    for (int j = 0; j < N; j++) jetCsv << ",h" << j;
    jetCsv << "\n";
    dropCsv << "t_us,z_um,R_um,v_ms\n";
    eventCsv << "t_us,type,z_um,R_um,v_ms\n";
    neckCsv << "event,t_us,hmin_um" << endl;
    jetCsv.flush(), dropCsv.flush(), eventCsv.flush(); // never leave an empty file (synced folders replace them)
    const double um = h0 * 1e6, us = tc * 1e6;

    cout << format("Oh = {:.4f}, We = {:.3f}, tc = {:.3f} us, f = {:.2f} kHz, N = {}\nOutput: {}\n",
                   Oh, U * U, us, omega / (2 * PI) / tc * 1e-3, N, dir.string());

    // Initial condition: hemispherical drop h = sqrt(h0^2 - z^2) at the nozzle, moving at vjet (paper, Sec. II)
    vector<double> init(2 * N + 1, U);
    for (int j = 0; j < N; j++) init[j] = 1.0 - pow((j + 0.5) / N, 2);
    init[iL] = 1.0;

    auto start = chrono::steady_clock::now();
    auto ta = buildJet(init);
    cout << format("Integrator built in {:.1f} s (order {})\n",
                   chrono::duration<double>(chrono::steady_clock::now() - start).count(), ta.get_order());

    vector<Drop> drops;
    int pinches = 0, neck = -1, caught = -1;
    auto onStep = [&](taylor_adaptive<double>& ta) {
        const double* s = ta.get_state().data();
        int j = findNeck(s, 0.5);
        if (j >= 0) neckCsv << pinches << "," << ta.get_time() * us << "," << sqrt(s[j]) * um << "\n";
        for (size_t i = 0; i < drops.size() && caught < 0; i++)
            if (drops[i].at(ta.get_time()) - drops[i].R() - s[iL] < gapMerge && s[iV + N - 1] > drops[i].v()) caught = int(i);
        return caught < 0;
    };

    for (double tOut = 0.0; tOut <= tEnd;) {
        ta.propagate_until(tOut, kw::callback = onStep, kw::max_delta_t = dtOut);
        double t = ta.get_time(), *s = ta.get_state_data();
        neck = findNeck(s, hPinch * (1.0 + 1e-6));
        if (!all_of(s, s + 2 * N + 1, [](double x) { return isfinite(x); }) || *min_element(s, s + N) <= 0.0) {
            cout << format("Solver failed at t = {:.2f} us\n", t * us);
            break;
        }
        if (neck >= 0) { // pinch-off: cells downstream of the neck become a drop
            const double dz = s[iL] / N;
            Drop d{0, 0, 0, t};
            for (int j = neck + 1; j < N; j++) {
                double dV = PI * s[j] * dz;
                d.V += dV;
                d.P += dV * 0.5 * (s[iV + j - 1] + s[iV + j]);
                d.z += dV * (j + 0.5) * dz;
            }
            d.z /= d.V;
            drops.push_back(d);
            remap(s, t, (neck + 0.5) * dz);
            eventCsv << format("{},pinch,{},{},{}", t * us, d.z * um, d.R() * um, d.v() * uc) << endl;
            cout << format("t = {:7.2f} us  pinch-off at z = {:6.1f} um   R = {:5.2f} um  v = {:.2f} m/s\n",
                           t * us, (neck + 0.5) * dz * um, d.R() * um, d.v() * uc);
            pinches++, neck = -1;
            continue;
        }
        if (caught >= 0) { // the jet tip has reached a slower drop: add its mass and momentum to the tip cells
            const Drop d = drops[caught];
            drops.erase(drops.begin() + caught);
            const double dz = s[iL] / N;
            const int n = clamp(int(2.0 * d.R() / dz), 1, N - 1);
            const double dA = d.V / (PI * n * dz);
            for (int j = N - n; j < N; j++) s[j] += dA;
            for (int k = N - n; k <= N; k++) { // mass-weighted face velocity
                double aOld = s[min(k, N - 1)] - dA;
                s[iV + k - 1] = (aOld * s[iV + k - 1] + dA * d.v()) / (aOld + dA);
            }
            eventCsv << format("{},absorbed,{},{},{}", t * us, d.at(t) * um, d.R() * um, d.v() * uc) << endl;
            caught = -1;
            continue;
        }

        // Output time: coalesce touching drops (momentum conserving), then write the frame
        sort(drops.begin(), drops.end(), [&](auto& p, auto& q) { return p.at(t) < q.at(t); });
        for (size_t i = 1; i < drops.size();) {
            Drop &p = drops[i - 1], &q = drops[i];
            if (q.at(t) - q.R() - p.at(t) - p.R() < gapMerge && p.v() > q.v()) {
                p = {p.V + q.V, p.P + q.P, (p.V * p.at(t) + q.V * q.at(t)) / (p.V + q.V), t};
                drops.erase(drops.begin() + i);
            } else i++;
        }
        jetCsv << t * us << "," << s[iL] * um;
        for (int j = 0; j < N; j++) jetCsv << "," << sqrt(s[j]) * um;
        jetCsv << endl;
        for (auto& d : drops) dropCsv << format("{},{},{},{}\n", t * us, d.at(t) * um, d.R() * um, d.v() * uc);
        dropCsv.flush(), neckCsv.flush();
        tOut += dtOut;
        // Simulation Heartbeat: log jet properties every 10 output steps
        if(fmod(tOut, dtOut) == 0){
            cout << format("t = {:7.2f} us  total jet length = {:6.1f} um\n",
                           t * us, s[iL] * um);
        }
    }
    cout << format("Done in {:.1f} s, {} pinch-offs\n",
                   chrono::duration<double>(chrono::steady_clock::now() - start).count(), pinches);
}
