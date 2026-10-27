// CLI for the faithful C++17 port.
//
// Usage:
//   dust_settling advection   v dt T Nx outdir
//   dust_settling laminar     total_time Nz a_lim R_AU a_micron zmax dust_density stop outdir
//   dust_settling turbulence  total_time Nz a_lim d_lim R_AU a_micron zmax D dust_density v_cap stop outdir
//
// Writes <outdir>/profile.csv and <outdir>/diagnostics.txt.
// All numeric arguments are validated strictly: exact argc per mode, finite
// doubles, integer-valued args with finite/range check BEFORE cast, positivity
// where physically required, bounded N, and a hard max-step safety guard.

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

#include "dust_settling.hpp"

using namespace dust;

namespace {
// Shared limits live in dust_settling.hpp (dust::MIN_GRID / MAX_GRID /
// MAX_STEPS); no local copies that could drift.

double parse_double(const std::string& s, const char* name) {
    if (s.empty()) throw std::runtime_error(std::string("missing value for ") + name);
    char* end = nullptr;
    errno = 0;
    double v = std::strtod(s.c_str(), &end);
    if (end == s.c_str() || *end != '\0') {
        throw std::runtime_error(std::string("invalid number for ") + name + ": '" + s + "'");
    }
    if (std::isnan(v) || std::isinf(v)) {
        throw std::runtime_error(std::string("non-finite value for ") + name);
    }
    return v;
}

long long parse_int(const std::string& s, const char* name, bool positive_only) {
    if (s.empty()) throw std::runtime_error(std::string("missing value for ") + name);
    char* end = nullptr;
    errno = 0;
    long long v = std::strtoll(s.c_str(), &end, 10);
    if (end == s.c_str() || *end != '\0' || errno == ERANGE) {
        throw std::runtime_error(std::string("invalid integer for ") + name + ": '" + s + "'");
    }
    if (positive_only && v <= 0) {
        throw std::runtime_error(std::string(name) + " must be positive");
    }
    if (v < 0 && !positive_only) {
        throw std::runtime_error(std::string(name) + " must be >= 0");
    }
    return v;
}

bool parse_bool(const std::string& s, const char* name) {
    long long v = parse_int(s, name, false);
    if (v > 1) throw std::runtime_error(std::string(name) + " must be 0 or 1");
    return v != 0;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            std::cerr << "usage: dust_settling <mode> <params...> <outdir>\n";
            return 2;
        }
        std::string mode = argv[1];

        if (mode == "advection") {
            if (argc != 7) {
                std::cerr << "advection usage: v dt T Nx outdir\n";
                return 2;
            }
            double v = parse_double(argv[2], "v");
            double dt = parse_double(argv[3], "dt");
            double T = parse_double(argv[4], "T");
            long long Nx = parse_int(argv[5], "Nx", true);
            if (dt <= 0.0) throw std::runtime_error("dt must be > 0");
            if (T < 0.0) throw std::runtime_error("T must be >= 0");
            if (Nx < MIN_GRID || Nx > MAX_GRID) throw std::runtime_error("Nx must be in [2, 1000000]");
            // checked_step_count enforces finite T/dt <= MAX_STEPS before any cast;
            // the resulting Advection ctor re-validates, so no silent huge run.
            (void)checked_step_count(T, dt, "advection");
            Advection s(v, dt, T, Nx);
            s.run();
            std::filesystem::create_directories(argv[6]);
            write_profile_csv("advection", argv[6], s.x, s.rhofinal, s.rho, s.x);
            write_diagnostics(argv[6], s.diag);
        } else if (mode == "laminar") {
            if (argc != 11) {
                std::cerr << "laminar usage: total_time Nz a_lim R_AU a_micron zmax dust_density stop outdir\n";
                return 2;
            }
            double total_time = parse_double(argv[2], "total_time");
            long long Nz = parse_int(argv[3], "Nz", true);
            double a_lim = parse_double(argv[4], "a_lim");
            double R_AU = parse_double(argv[5], "R_AU");
            double a_micron = parse_double(argv[6], "a_micron");
            double zmax = parse_double(argv[7], "zmax");
            double dust_density = parse_double(argv[8], "dust_density");
            bool stop = parse_bool(argv[9], "stop");
            if (total_time < 0.0) throw std::runtime_error("total_time must be >= 0");
            if (a_lim <= 0.0) throw std::runtime_error("a_lim must be > 0");
            if (R_AU <= 0.0) throw std::runtime_error("R_AU must be > 0");
            if (a_micron <= 0.0) throw std::runtime_error("a_micron must be > 0");
            if (zmax <= 0.0) throw std::runtime_error("zmax must be > 0");
            if (dust_density <= 0.0) throw std::runtime_error("dust_density must be > 0");
            if (Nz < MIN_GRID || Nz > MAX_GRID) throw std::runtime_error("Nz must be in [2, 1000000]");
            RealDiskDustSolver s(total_time, Nz, a_lim, R_AU, a_micron, zmax, dust_density);
            s.run(stop);
            std::filesystem::create_directories(argv[10]);
            write_profile_csv("laminar", argv[10], s.z, s.rhofinal, s.rho_initial, s.z_real);
            write_diagnostics(argv[10], s.diag);
        } else if (mode == "turbulence") {
            if (argc != 14) {
                std::cerr << "turbulence usage: total_time Nz a_lim d_lim R_AU a_micron zmax D dust_density v_cap stop outdir\n";
                return 2;
            }
            double total_time = parse_double(argv[2], "total_time");
            long long Nz = parse_int(argv[3], "Nz", true);
            double a_lim = parse_double(argv[4], "a_lim");
            double d_lim = parse_double(argv[5], "d_lim");
            double R_AU = parse_double(argv[6], "R_AU");
            double a_micron = parse_double(argv[7], "a_micron");
            double zmax = parse_double(argv[8], "zmax");
            double D = parse_double(argv[9], "D");
            double dust_density = parse_double(argv[10], "dust_density");
            double v_cap = parse_double(argv[11], "v_cap");
            bool stop = parse_bool(argv[12], "stop");
            if (total_time < 0.0) throw std::runtime_error("total_time must be >= 0");
            if (a_lim <= 0.0) throw std::runtime_error("a_lim must be > 0");
            if (d_lim <= 0.0) throw std::runtime_error("d_lim must be > 0");
            if (R_AU <= 0.0) throw std::runtime_error("R_AU must be > 0");
            if (a_micron <= 0.0) throw std::runtime_error("a_micron must be > 0");
            if (zmax <= 0.0) throw std::runtime_error("zmax must be > 0");
            if (D < 0.0) throw std::runtime_error("D must be >= 0");
            if (dust_density <= 0.0) throw std::runtime_error("dust_density must be > 0");
            if (Nz < MIN_GRID || Nz > MAX_GRID) throw std::runtime_error("Nz must be in [2, 1000000]");
            if (!(v_cap > 0.0) || std::isnan(v_cap) || std::isinf(v_cap))
                throw std::runtime_error("v_cap must be finite and > 0");
            TurbulenceSolver s(total_time, Nz, a_lim, d_lim, R_AU, a_micron, zmax, D, dust_density);
            s.v_cap_factor = v_cap;
            s.run(stop);
            std::filesystem::create_directories(argv[13]);
            write_profile_csv("turbulence", argv[13], s.z, s.rhofinal, s.rho_initial, s.z_real);
            write_diagnostics(argv[13], s.diag);
        } else {
            std::cerr << "unknown mode: " << mode << "\n";
            return 2;
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 2;
    }
    return 0;
}
