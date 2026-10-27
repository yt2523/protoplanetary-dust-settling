#pragma once
// Faithful C++17 port of the protoplanetary-dust-settling NumPy solvers.
// Scope (per PLAN.md):
//   - PureAdvection1D_v2  (Gaussian, signed, periodic, first-order upwind)
//   - RealDiskDustSolver  (laminar, scalar loop, velocity UNCLAMPED)
//   - DustSolver_Turbulence (vectorized weighted diffusion, velocity capped)
// Excluded: teaching DustSolver1D and legacy scalar turbulence.
//
// All numerical behaviour mirrors the original module line-by-line; no silent
// corrections.  Diagnostics record a pre-step boundary-flux integral so mass
// balance can be checked as |M_final - M_initial + sum(dt*(F_right-F_left))|.

#include <vector>
#include <string>
#include <cmath>
#include <limits>

namespace dust {

inline double nan_d() { return std::numeric_limits<double>::quiet_NaN(); }
inline double inf_d() { return std::numeric_limits<double>::infinity(); }

// Shared hard limits (mirrored by the CLI and used by the library at run/init).
constexpr long long MIN_GRID = 2;                 // N must be >= 2 (stencils index +-1)
constexpr long long MAX_GRID = 1000000LL;         // bounded memory
constexpr long long MAX_STEPS = 100000LL;         // safety guard against silent huge runs

// physical / disk constants shared by the two disk solvers
struct NaturalUnits {
    double omega = 0.0;
    double temperature = 0.0;
    double sound_speed = 0.0;
    double scale_height = 0.0;
    double rho_g0 = 0.0;
};

double year_to_second(double year);
NaturalUnits compute_natural_units(double R_AU);  // AU constants inside

// Validate a spatial grid size: 2 <= N <= 1,000,000, finite and exactly integral.
void check_grid_size(double N, const char* name);

// Single checked step-count helper: requires finite T >= 0, finite dt > 0, and a
// finite computed ratio T/dt <= MAX_STEPS, all verified BEFORE any cast.  The
// returned value is the exact (long long)truncation of T/dt, which is what both
// the Python `int(T/dt)` and the existing stencil use.  Called by every
// constructor / choose_time_step so the run loop never sees an invalid count.
long long checked_step_count(double T, double dt, const char* name);

// Per-run scalars recorded for validation.
struct Diagnostics {
    std::string mode;
    long long N = 0;
    double L = nan_d();
    double h = 0.0;
    double dt = 0.0;
    long long Nt = 0;              // planned number of steps = int(T/dt)
    double t_actual = 0.0;         // actual simulated time = actual_steps * dt
    long long actual_steps = 0;    // steps actually executed (<= Nt on early stop)
    double t_actual_years = nan_d();  // physical actual time in years
    double mass_initial = 0.0;
    double mass_final = 0.0;
    double cumulative_boundary_flux = 0.0;
    double mass_balance_error = 0.0;
    double min_rho = 0.0;          // min over ALL steps + initial (never clipped)
    double max_rho = 0.0;
    bool finite = true;            // finite over every step + initial
    double t_formation = nan_d();
    double t_formation_years = nan_d();
    double v_max = 0.0;
    double alpha = nan_d();
    double z_real_max = nan_d();
    double real_z_max = nan_d();
    double rho_initial_max = nan_d();
    double scale_height = nan_d();
    long long cap_count = 0;   // turbsolcapped velocities (v_cap_factor), from raw unclipped v
};

// ---------------- PureAdvection1D_v2 ----------------
class Advection {
public:
    Advection(double v, double dt, double T, long long Nx);

    void initial();
    std::vector<double> upwind(const std::vector<double>& r);
    void run();

    double v, dt, T;
    long long Nx;
    double h = 0.0;
    long long Nt = 0;
    std::vector<double> x, rho, rhofinal;
    Diagnostics diag;
};

// ---------------- RealDiskDustSolver (laminar) ----------------
class RealDiskDustSolver {
public:
    RealDiskDustSolver(double total_time, long long Nz, double a_lim,
                       double R_AU, double a_micron, double zmax, double dust_density = 1e3);

    void calculate_natural_unit();
    void choose_time_step();
    void initial_dust_points();
    // returns new rho; also reports the boundary fluxes used by this step
    // (F_right_top is the flux out of cell N-1, F_left_bottom is the flux into cell 0)
    std::vector<double> upwind(const std::vector<double>& r, double* F_right_top, double* F_left_bottom);
    void run(bool stop_when_cross);
    double get_t_formation_years() const;

    double total_time;
    long long Nz;
    double a_lim;
    double R_AU, a_micron, zmax, dust_density;
    double omega = 0.0, temperature = 0.0, sound_speed = 0.0, scale_height = 0.0, rho_g0 = 0.0;
    double T = 0.0, L = 0.0, h = 0.0, dt = 0.0;
    long long Nt = 0;
    std::vector<double> z, rho, rho_initial, rhofinal, v, z_real, real_rho;
    double t_formation = nan_d();
    double cumulative_boundary_flux = 0.0;
    Diagnostics diag;
};

// ---------------- DustSolver_Turbulence ----------------
class TurbulenceSolver {
public:
    TurbulenceSolver(double total_time, long long Nz, double a_lim, double d_lim,
                     double R_AU, double a_micron, double zmax, double D = 1e-5,
                     double dust_density = 1e3);

    void calculate_natural_unit();
    void choose_time_step();
    void initial_dust_points();
    std::vector<double> upwind(const std::vector<double>& r, double* F_right_top, double* F_left_bottom);
    void run(bool stop_when_cross);
    double get_t_formation_years() const;

    double total_time;
    long long Nz;
    double a_lim, d_lim;
    double R_AU, a_micron, zmax, dust_density;
    double omega = 0.0, temperature = 0.0, sound_speed = 0.0, scale_height = 0.0, rho_g0 = 0.0;
    double alpha = 0.0;
    double T = 0.0, L = 0.0, h = 0.0, dt = 0.0;
    long long Nt = 0;
    std::vector<double> z, rho, rho_initial, rhofinal, v, z_real, real_rho;
    double D;
    double t_formation = nan_d();
    double v_cap_factor = 1.0;
    double cumulative_boundary_flux = 0.0;
    long long cap_count = 0;
    Diagnostics diag;
};

// ---------------- output helpers ----------------
// coord / aux are physical columns; headers carry the mode-specific meaning.
void write_profile_csv(const std::string& mode, const std::string& outdir,
                       const std::vector<double>& coord,
                       const std::vector<double>& rho_final,
                       const std::vector<double>& rho_initial,
                       const std::vector<double>& aux);
void write_diagnostics(const std::string& outdir, const Diagnostics& d);

}  // namespace dust
