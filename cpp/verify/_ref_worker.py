#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Subprocess worker: run ONE original-Python reference case in isolation.

The acceptance harness runs the original Python solver in a bounded subprocess
(review fix #5) so a single case cannot hang the whole run; the parent enforces a
30 s wall-clock timeout with NO retry.  This worker does NOT reimplement any
evolution logic — it imports ``verify``'s ``run_python_*`` wrappers, which wrap the
*original* upstream ``upwind_method`` (see ``_instrument_upwind``).

Usage:
    python _ref_worker.py <output_prefix> <source_repo> <mode> <p1> <p2> ...

Writes:
    <output_prefix>.json   p_diag (JSON, default=str)
    <output_prefix>.npz    p_prof arrays (coord, rho_final, rho_init, aux)
"""
import os
import sys
import json

import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
if _here not in sys.path:
    sys.path.insert(0, _here)

import verify  # noqa: E402


def _load_float(s):
    return float(s)


def main():
    if len(sys.argv) < 4:
        sys.stderr.write("usage: _ref_worker.py <prefix> <source_repo> <mode> <args...>\n")
        return 2
    prefix = sys.argv[1]
    source_repo = sys.argv[2]
    mode = sys.argv[3]
    args = sys.argv[4:]

    verify._import_sources(source_repo)
    if mode == "advection":
        v, dt, T, Nx = _load_float(args[0]), _load_float(args[1]), _load_float(args[2]), int(args[3])
        p_diag, p_prof = verify.run_python_advection(v, dt, T, Nx)
    elif mode == "laminar":
        total_time, Nz, a_lim, R_AU, a_mic, zmax, dd, stop = args
        p_diag, p_prof = verify.run_python_laminar(
            _load_float(total_time), int(Nz), _load_float(a_lim), _load_float(R_AU),
            _load_float(a_mic), _load_float(zmax), _load_float(dd), bool(int(stop)))
    elif mode == "turbulence":
        total, Nz, a, d, R_AU, amic, zmax, D, dd, vcap, stop = args
        p_diag, p_prof = verify.run_python_turbulence(
            _load_float(total), int(Nz), _load_float(a), _load_float(d), _load_float(R_AU),
            _load_float(amic), _load_float(zmax), _load_float(D), _load_float(dd),
            _load_float(vcap), bool(int(stop)))
    else:
        sys.stderr.write("unknown mode %r\n" % mode)
        return 2

    with open(prefix + ".json", "w", encoding="utf-8") as f:
        # NaN / +/-inf are not valid JSON; serialize them as null and restore on
        # the parent side (all non-finite p_diag floats are diagnostic scalars).
        def _default(o):
            if isinstance(o, float) and (o != o or o in (float("inf"), float("-inf"))):
                return None
            raise TypeError("cannot serialize %r" % (type(o),))
        json.dump(p_diag, f, indent=2, default=_default)
    np.savez(prefix + ".npz", coord=p_prof[0], rho_final=p_prof[1],
             rho_init=p_prof[2], aux=p_prof[3])
    return 0


if __name__ == "__main__":
    sys.exit(main())
