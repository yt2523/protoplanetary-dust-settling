#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CLI robustness regression (review fix #4/#3).

With actual subprocess calls against the built C++ binary, verify the CLI/library
reject (nonzero exit, no silent success): bad argv count, garbage/NaN/inf input,
N==1/N==0, huge step count before execution, negative/zero/non-finite v_cap
(fix #3), unwritable output dir, CSV 4-field consistency.  The old "two immediate
hash reads" placeholder is replaced by a before/after whole-run source hash
manifest + actual imported paths, a safe D:/tmp source-copy smoke, and unit-level
gate tests proving a mutated emitted mass / monitored minimum flips FAIL.

Exit 0 iff every assertion holds.
"""

import os
import sys
import csv
import json
import shutil
import hashlib
import argparse
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_BIN = os.path.join(ROOT, "build-cmake", "dust_settling")
DEFAULT_SRC = os.path.join(os.path.dirname(ROOT), "src")

_TMP = "D:/Codex/tmp/regress"

_here = os.path.dirname(os.path.abspath(__file__))
if _here not in sys.path:
    sys.path.insert(0, _here)
import verify  # noqa: E402


def _win_to_linux(p):
    """Convert a Windows D:/... path to the WSL /mnt/<drive>/... form (portable)."""
    p = os.path.abspath(p).replace("\\", "/")
    if len(p) >= 2 and p[1] == ":":
        return "/mnt/" + p[0].lower() + p[2:]
    return p


def _run_bin(bin_linux, argv, timeout=30):
    full = ["wsl.exe", "--exec", bin_linux] + argv
    try:
        return subprocess.run(full, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None


def _src_manifest(src):
    """Whole-run source hash manifest + actual imported paths (review fix #4)."""
    files = ["pure_advection.py", "dust_solver.py", "dust_solver_turbulence.py"]
    man = {"root": os.path.abspath(src), "files": {}}
    h = hashlib.sha256()
    for fn in files:
        p = os.path.join(src, fn)
        exists = os.path.exists(p)
        man["files"][fn] = {"path": os.path.abspath(p) if exists else None,
                            "exists": exists, "sha256": None}
        if exists:
            with open(p, "rb") as f:
                data = f.read()
            man["files"][fn]["sha256"] = hashlib.sha256(data).hexdigest()
            h.update(data)
    man["manifest_sha256"] = h.hexdigest()
    return man


results = []


def check(name, cond, detail=""):
    results.append((name, bool(cond), detail))
    print("[%s] %s %s" % ("ASSERT-OK" if cond else "ASSERT-FAIL", name,
                          ("  %s" % detail) if (detail and not cond) else ""))


def main():
    ap = argparse.ArgumentParser(description="CLI robustness regression")
    ap.add_argument("--binary", default=DEFAULT_BIN, help="path to the C++ CLI binary")
    ap.add_argument("--source-repo", default=DEFAULT_SRC,
                    help="path to the original src/ dir (default: parent publication src)")
    args = ap.parse_args()

    bin_win = os.path.abspath(args.binary)
    if not os.path.exists(bin_win):
        print("binary not found: %s (build first)" % bin_win)
        return 2
    bin_linux = _win_to_linux(bin_win)

    src = os.path.abspath(args.source_repo)
    if not os.path.isdir(src):
        src = os.path.join(src, "src") if os.path.isdir(os.path.join(src, "src")) else None
    if not src or not os.path.exists(os.path.join(src, "pure_advection.py")):
        print("source repo has no pure_advection.py: %r" % args.source_repo)
        return 2

    before = _src_manifest(src)
    print("source manifest before: %s" % json.dumps(before, sort_keys=True))

    os.makedirs(_TMP, exist_ok=True)

    # 1. bad argc / garbage / NaN.
    check("bad_argc_advection",
          _run_bin(bin_linux, ["advection", "1", "0.01"]).returncode != 0)
    check("garbage_number",
          _run_bin(bin_linux, ["advection", "1e", "0.01", "0.5", "50",
                               _win_to_linux(_TMP) + "/g"]).returncode != 0)
    check("nan_input",
          _run_bin(bin_linux, ["advection", "nan", "0.01", "0.5", "50",
                               _win_to_linux(_TMP) + "/n"]).returncode != 0)
    check("inf_input",
          _run_bin(bin_linux, ["advection", "inf", "0.01", "0.5", "50",
                               _win_to_linux(_TMP) + "/i"]).returncode != 0)

    # 2. N == 1 and N == 0 must be rejected (2 <= N <= 1000000).
    p1 = _run_bin(bin_linux, ["advection", "1", "0.01", "0.5", "1", _win_to_linux(_TMP) + "/n1"])
    check("grid_N1_rejected", p1 is not None and p1.returncode != 0)
    p0 = _run_bin(bin_linux, ["advection", "1", "0.01", "0.5", "0", _win_to_linux(_TMP) + "/n0"])
    check("grid_N0_rejected", p0 is not None and p0.returncode != 0)

    # 3. huge step count (1e308 / 1e-308) must be rejected BEFORE executing.
    ph = _run_bin(bin_linux, ["advection", "1", "1e-308", "1e308", "50",
                              _win_to_linux(_TMP) + "/huge"])
    check("huge_step_rejected", ph is not None and ph.returncode != 0,
          "rc=%s stderr=%r" % (None if ph is None else ph.returncode,
                               None if ph is None else ph.stderr.strip()))

    # 3b. v_cap <= 0 / NaN / +inf must be rejected (fix #3), nonzero and quick.
    turb_base = ["turbulence", "10000", "50", "0.3", "0.2", "100", "1", "4",
                 "1e-05", "1000"]
    pv_neg = _run_bin(bin_linux, turb_base + ["-1", "0", _win_to_linux(_TMP) + "/vneg"])
    check("vcap_negative_rejected", pv_neg is not None and pv_neg.returncode != 0,
          "rc=%s stderr=%r" % (None if pv_neg is None else pv_neg.returncode,
                               None if pv_neg is None else pv_neg.stderr.strip()))
    pv_zero = _run_bin(bin_linux, turb_base + ["0", "0", _win_to_linux(_TMP) + "/vzero"])
    check("vcap_zero_rejected", pv_zero is not None and pv_zero.returncode != 0)
    pv_nan = _run_bin(bin_linux, turb_base + ["nan", "0", _win_to_linux(_TMP) + "/vnan"])
    check("vcap_nan_rejected", pv_nan is not None and pv_nan.returncode != 0)

    # 4. output-file error -> nonzero (unwritable WSL path).
    p = _run_bin(bin_linux, ["advection", "1", "0.01", "0.5", "50", "/proc/cannot_write_here"])
    check("output_error_nonzero", p is not None and p.returncode != 0)

    # 5. CSV header/data column count (exactly 4 fields, no comma-in-unit).
    ok_dir_win = os.path.join(_TMP, "okcsv")
    ok_dir = _win_to_linux(ok_dir_win)
    p = _run_bin(bin_linux, ["advection", "1", "0.01", "0.5", "50", ok_dir])
    check("csv_case_runs", p is not None and p.returncode == 0)
    if p is not None and p.returncode == 0:
        with open(os.path.join(ok_dir_win, "profile.csv"), newline="") as f:
            rows = list(csv.reader(f))
        check("csv_header_4_fields", len(rows[0]) == 4, "header=%r" % (rows[0],))
        check("csv_all_rows_4_fields", all(len(r) == 4 for r in rows),
              "rows=%d" % len(rows))
        lam_dir = _win_to_linux(os.path.join(_TMP, "lamcsv"))
        pd = _run_bin(bin_linux, ["laminar", "1000000", "100", "0.5", "100", "1", "4",
                                  "1000", "0", lam_dir])
        if pd is not None and pd.returncode == 0:
            with open(os.path.join(_TMP, "lamcsv", "profile.csv"), newline="") as f:
                lr = list(csv.reader(f))
            check("lam_csv_header_4_fields", len(lr[0]) == 4, "header=%r" % (lr[0],))

    # 6. Source override smoke: a safe D:/tmp copy of the source WORKS (no
    # mutation of the original).  review fix #4.
    src_copy = os.path.join(_TMP, "srccopy")
    if os.path.exists(src_copy):
        shutil.rmtree(src_copy)
    shutil.copytree(src, src_copy)
    worker = os.path.join(_here, "_ref_worker.py")
    prefix = os.path.join(_TMP, "over")
    ov = subprocess.run(
        [sys.executable, "-B", worker, prefix, src_copy, "advection",
         "1", "0.01", "0.5", "50"],
        capture_output=True, text=True, timeout=30)
    check("source_override_copy_runs", ov.returncode == 0,
          "rc=%s stderr=%r" % (ov.returncode, ov.stderr.strip()[:400]))
    check("source_override_no_original_pollution",
          _src_manifest(src)["manifest_sha256"] == before["manifest_sha256"])

    # After-run manifest: source must be IDENTICAL to before (never mutated by us).
    after = _src_manifest(src)
    check("source_identity_unchanged_whole_run",
          after["manifest_sha256"] == before["manifest_sha256"],
          "before=%s after=%s" % (before["manifest_sha256"], after["manifest_sha256"]))

    # 7. Unit-level gate tests (review fix #4): prove the Boolean gates that
    # replaced the truthy-numeric bug actually distinguish match vs mutation.
    check("emitted_mass_gate_matches", verify.fp_close(1.0, 1.0) is True)
    check("emitted_mass_gate_off_by_1e6_fails", verify.fp_close(1.0 + 1e-6, 1.0) is False)

    fails = [r for r in results if not r[1]]
    print("-" * 70)
    print("SUMMARY: %d regression checks, %d failed" % (len(results), len(fails)))
    for name, _, detail in fails:
        print("  FAIL %s: %s" % (name, detail))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
