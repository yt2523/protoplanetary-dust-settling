#include "dust_settling.hpp"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <stdexcept>

namespace dust {

namespace {
constexpr double G = 6.67430e-11;
constexpr double M_SUN = 1.98847e30;
constexpr double AU = 1.495978707e11;
constexpr double KB = 1.380649e-23;
constexpr double M_H = 1.6735575e-27;
constexpr double PI = 3.14159265358979323846;

std::string f17(double x) {
    if (std::isnan(x)) return "nan";
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.17g", x);
    return std::string(buf);
}
}  // namespace

double year_to_second(double year) {
    return year * 365.25 * 24.0 * 3600.0;
}

void check_grid_size(double N, const char* name) {
    if (std::isnan(N) || std::isinf(N) || N >= 9.007199254740992e15) {
        throw std::runtime_error(std::string(name) + " is not a finite integral grid size");
    }
    long long Ni = (long long)N;
    if ((double)Ni != N) {
        throw std::runtime_error(std::string(name) + " must be an integer");
    }
    if (Ni < MIN_GRID || Ni > MAX_GRID) {
        throw std::runtime_error(std::string(name) + " must be in [2, 1000000]");
    }
}

long long checked_step_count(double T, double dt, const char* name) {
    if (std::isnan(T) || std::isinf(T) || T < 0.0) {
        throw std::runtime_error(std::string(name) + ": T must be finite and >= 0");
    }
    if (std::isnan(dt) || std::isinf(dt) || dt <= 0.0) {
        throw std::runtime_error(std::string(name) + ": dt must be finite and > 0");
    }
    double ratio = T / dt;
    if (std::isnan(ratio) || std::isinf(ratio)) {
        throw std::runtime_error(std::string(name) + ": T/dt is not finite");
    }
    if (ratio > (double)MAX_STEPS) {
        throw std::runtime_error(std::string(name) +
                                 ": step count exceeds safety guard (100000)");
    }
    return (long long)ratio;
}

NaturalUnits compute_natural_units(double R_AU) {
    NaturalUnits u;
    u.omega = std::sqrt(G * M_SUN / std::pow(R_AU * AU, 3.0));
    u.temperature = 264.0 * std::pow(R_AU, -0.5);
    u.sound_speed = std::sqrt(KB * u.temperature / (2.0 * M_H));
    u.scale_height = u.sound_speed / u.omega;
    u.rho_g0 = 5e-7 * std::pow(R_AU, -3.0);
    return u;
}

// ----------------------------------------------------------------- Advection
Advection::Advection(double v_, double dt_, double T_, long long Nx_)
    : v(v_), dt(dt_), T(T_), Nx(Nx_), Nt(0) {
    check_grid_size((double)Nx_, "Nx");
    if (std::isnan(v_) || std::isinf(v_)) throw std::runtime_error("v must be finite");
    Nt = checked_step_count(T_, dt_, "advection");
}

void Advection::initial() {
    double Length = 1.0;
    long long N = Nx;
    h = Length / (double)N;
    if (h <= 0.0) {
        throw std::runtime_error("Spatial step size invalid, check initial_function_points");
    }
    x.resize(N);
    double step = Length / (double)N;
    for (long long i = 0; i < N; ++i) x[i] = step * (double)i;

    double x0 = 0.5, sigma = 0.07;
    rho.resize(N);
    for (long long i = 0; i < N; ++i) {
        double d = (x[i] - x0) / sigma;
        rho[i] = std::exp(-0.5 * d * d);
    }
    double mx = 0.0;
    for (long long i = 0; i < N; ++i) mx = std::max(mx, rho[i]);
    for (long long i = 0; i < N; ++i) rho[i] /= mx;
}

std::vector<double> Advection::upwind(const std::vector<double>& r) {
    long long n = (long long)r.size();
    double a = v * dt / h;
    double abs_a = std::fabs(a);
    if (abs_a < 1e-13) {
        return r;  // copy, no change
    }
    if (abs_a > 1.0) {
        throw std::runtime_error("invalid advective number, upwind_method do not work.");
    } else if (abs_a > 0.75) {
        std::fprintf(stderr, "advective number %.3f too large, may cause big error.\n", abs_a);
    } else {
        std::fprintf(stderr, "advective number %.3f acceptable, proceed simulation.\n", abs_a);
    }

    std::vector<double> rn(n, 0.0);
    if (a < 0.0) {
        for (long long i = 0; i < n - 1; ++i) {
            rn[i] = r[i] * (1.0 - abs_a) + r[i + 1] * abs_a;
        }
        rn[n - 1] = r[n - 1] * (1.0 - abs_a) + r[0] * abs_a;
    }
    if (a > 0.0) {
        rn[0] = r[0] * (1.0 - abs_a) + r[n - 1] * abs_a;
        for (long long i = 1; i < n; ++i) {
            rn[i] = r[i] * (1.0 - abs_a) + r[i - 1] * abs_a;
        }
    }
    return rn;
}

void Advection::run() {
    initial();
    std::vector<double> r = rho;  // copy
    long long steps = Nt;

    // Track min/max/finite over EVERY step including the initial condition.
    double rmin = rho[0], rmax = rho[0];
    bool fin = true;
    double rmax_initial = rho[0];
    for (long long i = 0; i < Nx; ++i) {
        rmin = std::min(rmin, rho[i]);
        rmax = std::max(rmax, rho[i]);
        rmax_initial = std::max(rmax_initial, rho[i]);
        if (!std::isfinite(rho[i])) fin = false;
    }

    for (long long s = 0; s < steps; ++s) {
        r = upwind(r);
        for (long long i = 0; i < Nx; ++i) {
            rmin = std::min(rmin, r[i]);
            rmax = std::max(rmax, r[i]);
            if (!std::isfinite(r[i])) fin = false;
        }
    }
    rhofinal = r;

    diag.mode = "advection";
    diag.N = Nx;
    diag.L = 1.0;
    diag.h = h;
    diag.dt = dt;
    diag.Nt = Nt;
    diag.actual_steps = Nt;
    diag.t_actual = (double)diag.actual_steps * dt;
    diag.t_actual_years = nan_d();

    double s0 = 0.0, s1 = 0.0;
    for (long long i = 0; i < Nx; ++i) s0 += rho[i];
    for (long long i = 0; i < Nx; ++i) s1 += rhofinal[i];
    diag.mass_initial = h * s0;
    diag.mass_final = h * s1;
    diag.cumulative_boundary_flux = 0.0;
    diag.mass_balance_error = diag.mass_final - diag.mass_initial;

    diag.min_rho = rmin;
    diag.max_rho = rmax;
    diag.finite = fin;
    diag.t_formation = nan_d();
    diag.t_formation_years = nan_d();
    diag.v_max = std::fabs(v);
    diag.alpha = nan_d();
    diag.z_real_max = nan_d();
    diag.real_z_max = nan_d();
    diag.rho_initial_max = rmax_initial;
    diag.scale_height = nan_d();
}

// ------------------------------------------------- RealDiskDustSolver (laminar)
RealDiskDustSolver::RealDiskDustSolver(double total_time_, long long Nz_, double a_lim_,
                                       double R_AU_, double a_micron_, double zmax_,
                                       double dust_density_)
    : total_time(total_time_), Nz(Nz_), a_lim(a_lim_), R_AU(R_AU_),
      a_micron(a_micron_), zmax(zmax_), dust_density(dust_density_) {
    check_grid_size((double)Nz_, "Nz");
    if (!std::isfinite(total_time_) || total_time_ < 0.0)
        throw std::runtime_error("total_time must be finite and >= 0");
}

void RealDiskDustSolver::calculate_natural_unit() {
    NaturalUnits u = compute_natural_units(R_AU);
    omega = u.omega;
    temperature = u.temperature;
    sound_speed = u.sound_speed;
    scale_height = u.scale_height;
    rho_g0 = u.rho_g0;

    double real_total_time_s = year_to_second(total_time);
    T = real_total_time_s * omega;

    double max_L = 4.0;
    L = (zmax * AU) / scale_height;
    L = std::min(L, max_L);

    h = L / (double)Nz;
    if (h <= 0.0) {
        throw std::runtime_error("Spatial step size is not valid.check initial_function_points");
    }
}

void RealDiskDustSolver::choose_time_step() {
    double v_max = 0.0;
    for (double vi : v) v_max = std::max(v_max, std::fabs(vi));
    if (v_max == 0.0) {
        throw std::runtime_error("All velocities are zero, cannot choose dt from CFL.");
    }
    double max_dt = a_lim * h / v_max;
    dt = max_dt;
    Nt = checked_step_count(T, dt, "laminar");
}

void RealDiskDustSolver::initial_dust_points() {
    double real_z_max = zmax * AU;
    z.resize(Nz);
    z_real.resize(Nz);
    double z_step = L / (double)Nz;
    double zr_step = real_z_max / (double)Nz;
    for (long long i = 0; i < Nz; ++i) {
        z[i] = z_step * (double)i;
        z_real[i] = zr_step * (double)i;
    }

    rho_initial.resize(Nz);
    real_rho.resize(Nz);
    for (long long i = 0; i < Nz; ++i) {
        rho_initial[i] = std::exp(-0.5 * z[i] * z[i]);
        double q = z_real[i] / scale_height;
        real_rho[i] = rho_g0 * std::exp(-0.5 * q * q);
    }
    rho.resize(Nz);
    for (long long i = 0; i < Nz; ++i) {
        rho[i] = 0.01 * rho_initial[i];
    }

    v.resize(Nz);
    double real_dust_a = a_micron * 1e-6;
    for (long long i = 0; i < Nz; ++i) {
        double ts = std::sqrt(PI / 8.0) * dust_density * real_dust_a / (real_rho[i] * sound_speed);
        double v_real = -(omega * omega) * z_real[i] * ts;
        // NOTE: laminar velocity is intentionally NOT clipped (quirk preserved).
        v[i] = v_real / (omega * scale_height);
    }
}

std::vector<double> RealDiskDustSolver::upwind(const std::vector<double>& r,
                                               double* F_right_top, double* F_left_bottom) {
    long long n = (long long)r.size();
    if (n == 0) throw std::runtime_error("rho array is empty.");
    if ((long long)v.size() != n) {
        throw std::runtime_error("velocity array length not equal to rho array length");
    }
    std::vector<double> rn(n, 0.0);
    double dt_over_h = dt / h;
    double F_right, F_left;

    // i = 0
    if (v[1] >= 0.0) {
        F_right = v[1] * r[0];
    } else {
        F_right = v[1] * r[1];
    }
    F_left = 0.0;
    rn[0] = r[0] - dt_over_h * (F_right - F_left);

    // i = 1 .. n-2
    for (long long i = 1; i < n - 1; ++i) {
        if (v[i + 1] >= 0.0) {
            F_right = v[i + 1] * r[i];
        } else {
            F_right = v[i + 1] * r[i + 1];
        }
        if (v[i] >= 0.0) {
            F_left = v[i] * r[i - 1];
        } else {
            F_left = v[i] * r[i];
        }
        rn[i] = r[i] - dt_over_h * (F_right - F_left);
    }

    // i = n-1
    if (v[n - 1] >= 0.0) {
        F_right = v[n - 1] * r[n - 1];
        F_left = v[n - 1] * r[n - 2];
    } else {
        F_right = 0.0;
        F_left = v[n - 1] * r[n - 1];
    }
    rn[n - 1] = r[n - 1] - dt_over_h * (F_right - F_left);

    *F_right_top = F_right;
    *F_left_bottom = 0.0;
    return rn;
}

void RealDiskDustSolver::run(bool stop_when_cross) {
    calculate_natural_unit();
    initial_dust_points();
    std::vector<double> r = rho;  // dust initial
    choose_time_step();

    cumulative_boundary_flux = 0.0;
    t_formation = nan_d();
    long long steps_done = Nt;

    double rmin = rho[0], rmax = rho[0];
    bool fin = true;
    for (long long i = 0; i < Nz; ++i) {
        rmin = std::min(rmin, rho[i]);
        rmax = std::max(rmax, rho[i]);
        if (!std::isfinite(rho[i])) fin = false;
    }

    for (long long n = 0; n < Nt; ++n) {
        double Frt = 0.0, Flb = 0.0;
        r = upwind(r, &Frt, &Flb);
        cumulative_boundary_flux += dt * (Frt - Flb);
        for (long long i = 0; i < Nz; ++i) {
            rmin = std::min(rmin, r[i]);
            rmax = std::max(rmax, r[i]);
            if (!std::isfinite(r[i])) fin = false;
        }

        if (std::isnan(t_formation)) {
            bool any = false;
            for (long long i = 0; i < Nz; ++i) {
                if (r[i] / rho_initial[i] > 1.0) {
                    any = true;
                    break;
                }
            }
            if (any) {
                t_formation = (double)(n + 1) * dt;
                if (stop_when_cross) {
                    steps_done = n + 1;
                    break;
                }
            }
        }
    }

    rhofinal = r;
    double s0 = 0.0, s1 = 0.0;
    for (long long i = 0; i < Nz; ++i) s0 += rho[i];
    for (long long i = 0; i < Nz; ++i) s1 += rhofinal[i];

    diag.mode = "laminar";
    diag.N = Nz;
    diag.L = L;
    diag.h = h;
    diag.dt = dt;
    diag.Nt = Nt;
    diag.actual_steps = steps_done;
    diag.t_actual = (double)steps_done * dt;
    diag.t_actual_years = (double)steps_done * dt / omega / (365.25 * 24.0 * 3600.0);
    diag.mass_initial = h * s0;
    diag.mass_final = h * s1;
    diag.cumulative_boundary_flux = cumulative_boundary_flux;
    diag.mass_balance_error = diag.mass_final - diag.mass_initial + cumulative_boundary_flux;
    diag.min_rho = rmin;
    diag.max_rho = rmax;
    diag.finite = fin;
    diag.t_formation = t_formation;
    diag.t_formation_years = get_t_formation_years();
    double vmx = 0.0;
    for (double vi : v) vmx = std::max(vmx, std::fabs(vi));
    diag.v_max = vmx;
    diag.alpha = nan_d();
    double zrmx = z_real[0], rimx = rho_initial[0];
    for (long long i = 0; i < Nz; ++i) {
        zrmx = std::max(zrmx, z_real[i]);
        rimx = std::max(rimx, rho_initial[i]);
    }
    diag.z_real_max = zrmx;
    diag.real_z_max = zmax * AU;
    diag.rho_initial_max = rimx;
    diag.scale_height = scale_height;
}

double RealDiskDustSolver::get_t_formation_years() const {
    if (std::isnan(t_formation)) return nan_d();
    double t_seconds = t_formation / omega;
    return t_seconds / (365.25 * 24.0 * 3600.0);
}

// ------------------------------------------------- DustSolver_Turbulence
TurbulenceSolver::TurbulenceSolver(double total_time_, long long Nz_, double a_lim_,
                                   double d_lim_, double R_AU_, double a_micron_,
                                   double zmax_, double D_, double dust_density_)
    : total_time(total_time_), Nz(Nz_), a_lim(a_lim_), d_lim(d_lim_), R_AU(R_AU_),
      a_micron(a_micron_), zmax(zmax_), dust_density(dust_density_), D(D_) {
    check_grid_size((double)Nz_, "Nz");
    if (!std::isfinite(total_time_) || total_time_ < 0.0)
        throw std::runtime_error("total_time must be finite and >= 0");
}

void TurbulenceSolver::calculate_natural_unit() {
    NaturalUnits u = compute_natural_units(R_AU);
    omega = u.omega;
    temperature = u.temperature;
    sound_speed = u.sound_speed;
    scale_height = u.scale_height;
    rho_g0 = u.rho_g0;

    alpha = D * omega / (sound_speed * sound_speed);

    double real_total_time_s = year_to_second(total_time);
    T = real_total_time_s * omega;

    double max_L = 4.0;
    L = (zmax * AU) / scale_height;
    L = std::min(L, max_L);

    h = L / (double)Nz;
    if (h <= 0.0) {
        throw std::runtime_error("Spatial step size is not valid.check initial_function_points");
    }
}

void TurbulenceSolver::choose_time_step() {
    double v_max = 0.0;
    for (double vi : v) v_max = std::max(v_max, std::fabs(vi));
    if (v_max == 0.0) {
        throw std::runtime_error("All velocities are zero, cannot choose dt from CFL.");
    }
    double max_dt_adv = a_lim * h / v_max;

    double density_max = rho_initial[0];
    for (double ri : rho_initial) density_max = std::max(density_max, ri);

    // D == 0 makes alpha == 0, so the Python expression raises a divide-by-zero
    // RuntimeWarning and yields +inf.  Reproduce that explicitly.
    double max_dt_diff;
    if (alpha * density_max == 0.0) {
        max_dt_diff = inf_d();
        std::fprintf(stderr, "warning: diffusion dt denominator zero (D=0), dt_diff=inf\n");
    } else {
        max_dt_diff = d_lim * h * h / (alpha * density_max);
    }

    dt = std::min(max_dt_adv, max_dt_diff);
    Nt = checked_step_count(T, dt, "turbulence");
}

void TurbulenceSolver::initial_dust_points() {
    z.resize(Nz);
    z_real.resize(Nz);
    double z_step = L / (double)Nz;
    for (long long i = 0; i < Nz; ++i) {
        z[i] = z_step * (double)i;
        z_real[i] = z[i] * scale_height;  // NOTE: z * H, NOT linspace (quirk)
    }

    rho_initial.resize(Nz);
    real_rho.resize(Nz);
    for (long long i = 0; i < Nz; ++i) {
        rho_initial[i] = std::exp(-0.5 * z[i] * z[i]);
        double q = z_real[i] / scale_height;
        real_rho[i] = rho_g0 * std::exp(-0.5 * q * q);
    }
    rho.resize(Nz);
    for (long long i = 0; i < Nz; ++i) {
        rho[i] = 0.01 * rho_initial[i];
    }

    v.resize(Nz);
    double real_dust_a = a_micron * 1e-6;
    long long cap_count = 0;
    for (long long i = 0; i < Nz; ++i) {
        double ts = std::sqrt(PI / 8.0) * dust_density * real_dust_a / (real_rho[i] * sound_speed);
        double v_real = -(omega * omega) * z_real[i] * ts;
        double v_cap = v_cap_factor * sound_speed;
        if (v_real < -v_cap) { v_real = -v_cap; ++cap_count; }
        else if (v_real > v_cap) { v_real = v_cap; ++cap_count; }
        v[i] = v_real / (omega * scale_height);
    }
    this->cap_count = cap_count;
}

std::vector<double> TurbulenceSolver::upwind(const std::vector<double>& r,
                                             double* F_right_top, double* F_left_bottom) {
    long long n = (long long)r.size();
    if (n == 0) throw std::runtime_error("rho array is empty.");
    if ((long long)v.size() != n) {
        throw std::runtime_error("velocity array length not equal to rho array length");
    }

    std::vector<double> ratio(n, 0.0);
    for (long long i = 0; i < n; ++i) ratio[i] = r[i] / rho_initial[i];

    std::vector<double> rn(n, 0.0);
    double dt_over_h = dt / h;

    double F_right_adv, F_left_adv, F_right_diff, F_left_diff;
    double right_slope, left_slope, rho_g_right, rho_g_left;
    double F_right_total, F_left_total;

    // i = 0
    if (v[1] >= 0.0) F_right_adv = v[1] * r[0];
    else F_right_adv = v[1] * r[1];
    F_left_adv = 0.0;
    right_slope = (ratio[1] - ratio[0]) / h;
    rho_g_right = 0.5 * (rho_initial[1] + rho_initial[0]);
    F_right_diff = -alpha * rho_g_right * right_slope;
    F_left_diff = 0.0;
    F_right_total = F_right_adv + F_right_diff;
    F_left_total = F_left_adv + F_left_diff;
    rn[0] = r[0] - dt_over_h * (F_right_total - F_left_total);

    // interior
    if (n > 2) {
        for (long long i = 1; i < n - 1; ++i) {
            double v_right = v[i + 1];
            double v_left = v[i];
            F_right_adv = (v_right >= 0.0) ? v_right * r[i] : v_right * r[i + 1];
            F_left_adv = (v_left >= 0.0) ? v_left * r[i - 1] : v_left * r[i];

            right_slope = (ratio[i + 1] - ratio[i]) / h;
            left_slope = (ratio[i] - ratio[i - 1]) / h;
            rho_g_right = 0.5 * (rho_initial[i + 1] + rho_initial[i]);
            rho_g_left = 0.5 * (rho_initial[i] + rho_initial[i - 1]);
            F_right_diff = -alpha * rho_g_right * right_slope;
            F_left_diff = -alpha * rho_g_left * left_slope;

            F_right_total = F_right_adv + F_right_diff;
            F_left_total = F_left_adv + F_left_diff;
            rn[i] = r[i] - dt_over_h * (F_right_total - F_left_total);
        }
    }

    // i = n-1
    long long il = n - 1;
    if (v[il] >= 0.0) {
        F_right_adv = v[il] * r[il];
        F_left_adv = v[il] * r[il - 1];
    } else {
        F_right_adv = 0.0;
        F_left_adv = v[il] * r[il];
    }
    right_slope = (0.0 - ratio[il]) / h;
    rho_g_right = 0.5 * (0.0 + rho_initial[il]);
    left_slope = (ratio[il] - ratio[il - 1]) / h;
    rho_g_left = 0.5 * (rho_initial[il] + rho_initial[il - 1]);
    F_right_diff = -alpha * rho_g_right * right_slope;
    F_left_diff = -alpha * rho_g_left * left_slope;
    F_right_total = F_right_adv + F_right_diff;
    F_left_total = F_left_adv + F_left_diff;
    rn[il] = r[il] - dt_over_h * (F_right_total - F_left_total);

    *F_right_top = F_right_total;
    *F_left_bottom = 0.0;
    return rn;
}

void TurbulenceSolver::run(bool stop_when_cross) {
    // v_cap_factor is public/mutable; reject non-finite / non-positive values at
    // run-time before initial_dust_points() uses it to clip velocities (a <=0 cap
    // would invert the clip direction).  Mirrors the CLI check in main.cpp.
    if (!(v_cap_factor > 0.0) || std::isnan(v_cap_factor) || std::isinf(v_cap_factor)) {
        throw std::runtime_error("v_cap_factor must be finite and > 0");
    }
    calculate_natural_unit();
    initial_dust_points();
    std::vector<double> r = rho;
    choose_time_step();

    cumulative_boundary_flux = 0.0;
    t_formation = nan_d();
    long long steps_done = Nt;

    double rmin = rho[0], rmax = rho[0];
    bool fin = true;
    for (long long i = 0; i < Nz; ++i) {
        rmin = std::min(rmin, rho[i]);
        rmax = std::max(rmax, rho[i]);
        if (!std::isfinite(rho[i])) fin = false;
    }

    for (long long n = 0; n < Nt; ++n) {
        double Frt = 0.0, Flb = 0.0;
        r = upwind(r, &Frt, &Flb);
        cumulative_boundary_flux += dt * (Frt - Flb);
        for (long long i = 0; i < Nz; ++i) {
            rmin = std::min(rmin, r[i]);
            rmax = std::max(rmax, r[i]);
            if (!std::isfinite(r[i])) fin = false;
        }

        if (std::isnan(t_formation)) {
            bool any = false;
            for (long long i = 0; i < Nz; ++i) {
                if (r[i] / rho_initial[i] > 1.0) {
                    any = true;
                    break;
                }
            }
            if (any) {
                t_formation = (double)(n + 1) * dt;
                if (stop_when_cross) {
                    steps_done = n + 1;
                    break;
                }
            }
        }
    }

    rhofinal = r;
    double s0 = 0.0, s1 = 0.0;
    for (long long i = 0; i < Nz; ++i) s0 += rho[i];
    for (long long i = 0; i < Nz; ++i) s1 += rhofinal[i];

    diag.mode = "turbulence";
    diag.N = Nz;
    diag.L = L;
    diag.h = h;
    diag.dt = dt;
    diag.Nt = Nt;
    diag.actual_steps = steps_done;
    diag.t_actual = (double)steps_done * dt;
    diag.t_actual_years = (double)steps_done * dt / omega / (365.25 * 24.0 * 3600.0);
    diag.mass_initial = h * s0;
    diag.mass_final = h * s1;
    diag.cumulative_boundary_flux = cumulative_boundary_flux;
    diag.mass_balance_error = diag.mass_final - diag.mass_initial + cumulative_boundary_flux;
    diag.min_rho = rmin;
    diag.max_rho = rmax;
    diag.finite = fin;
    diag.t_formation = t_formation;
    diag.t_formation_years = get_t_formation_years();
    double vmx = 0.0;
    for (double vi : v) vmx = std::max(vmx, std::fabs(vi));
    diag.v_max = vmx;
    diag.alpha = alpha;
    double zrmx = z_real[0], rimx = rho_initial[0];
    for (long long i = 0; i < Nz; ++i) {
        zrmx = std::max(zrmx, z_real[i]);
        rimx = std::max(rimx, rho_initial[i]);
    }
    diag.z_real_max = zrmx;
    diag.real_z_max = zmax * AU;
    diag.rho_initial_max = rimx;
    diag.scale_height = scale_height;
    diag.cap_count = cap_count;
}

double TurbulenceSolver::get_t_formation_years() const {
    if (std::isnan(t_formation)) return nan_d();
    double t_seconds = t_formation / omega;
    return t_seconds / (365.25 * 24.0 * 3600.0);
}

// ---------------- output helpers ----------------
void write_profile_csv(const std::string& mode, const std::string& outdir,
                       const std::vector<double>& coord,
                       const std::vector<double>& rho_final,
                       const std::vector<double>& rho_initial,
                       const std::vector<double>& aux) {
    if (coord.size() != rho_final.size() || coord.size() != rho_initial.size() ||
        coord.size() != aux.size()) {
        throw std::runtime_error("profile array size mismatch");
    }
    std::ofstream f(outdir + "/profile.csv");
    if (!f.is_open()) {
        throw std::runtime_error("cannot open output file: " + outdir + "/profile.csv");
    }
    f.exceptions(std::ios::badbit | std::ios::failbit);
    // Exactly 4 comma-separated fields per row; units carry no commas so the
    // field count is unambiguous (review: header comma-in-unit mismatch).
    if (mode == "advection") {
        f << "x_over_L,rho_final_normalized,rho_initial_normalized,x_over_L_aux\n";
    } else {
        // 4th column is the physical z coordinate (z_real, in metres).
        f << "z_over_H,rho_d_over_rho_g0,rho_g_over_rho_g0,z_coordinate_m\n";
    }
    f.precision(17);
    for (size_t i = 0; i < coord.size(); ++i) {
        f << coord[i] << "," << rho_final[i] << "," << rho_initial[i] << "," << aux[i] << "\n";
    }
    f.flush();
    if (!f.good()) {
        throw std::runtime_error("write failed for: " + outdir + "/profile.csv");
    }
}

void write_diagnostics(const std::string& outdir, const Diagnostics& d) {
    std::ofstream f(outdir + "/diagnostics.txt");
    if (!f.is_open()) {
        throw std::runtime_error("cannot open output file: " + outdir + "/diagnostics.txt");
    }
    f.exceptions(std::ios::badbit | std::ios::failbit);
    f.precision(17);
    f << "mode " << d.mode << "\n";
    f << "N " << d.N << "\n";
    f << "L " << f17(d.L) << "\n";
    f << "h " << f17(d.h) << "\n";
    f << "dt " << f17(d.dt) << "\n";
    f << "Nt " << d.Nt << "\n";
    f << "actual_steps " << d.actual_steps << "\n";
    f << "t_actual " << f17(d.t_actual) << "\n";
    f << "t_actual_years " << f17(d.t_actual_years) << "\n";
    f << "mass_initial " << f17(d.mass_initial) << "\n";
    f << "mass_final " << f17(d.mass_final) << "\n";
    f << "cumulative_boundary_flux " << f17(d.cumulative_boundary_flux) << "\n";
    f << "mass_balance_error " << f17(d.mass_balance_error) << "\n";
    f << "min_rho " << f17(d.min_rho) << "\n";
    f << "max_rho " << f17(d.max_rho) << "\n";
    f << "finite " << (d.finite ? 1 : 0) << "\n";
    f << "t_formation " << f17(d.t_formation) << "\n";
    f << "t_formation_years " << f17(d.t_formation_years) << "\n";
    f << "v_max " << f17(d.v_max) << "\n";
    f << "alpha " << f17(d.alpha) << "\n";
    f << "z_real_max " << f17(d.z_real_max) << "\n";
    f << "real_z_max " << f17(d.real_z_max) << "\n";
    f << "rho_initial_max " << f17(d.rho_initial_max) << "\n";
    f << "scale_height " << f17(d.scale_height) << "\n";
    f << "cap_count " << d.cap_count << "\n";
    f.flush();
    if (!f.good()) {
        throw std::runtime_error("write failed for: " + outdir + "/diagnostics.txt");
    }
}

}  // namespace dust
