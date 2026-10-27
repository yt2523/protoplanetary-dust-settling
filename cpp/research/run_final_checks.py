"""Deterministic final-check runner (authored for research scope).

Runs driver, tests, and capacity audit in order, capturing artifacts.
"""
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")


def _env():
    e = os.environ.copy()
    e["TEMP"] = "D:/Codex/tmp"
    e["TMP"] = "D:/Codex/tmp"
    e["MPLCONFIGDIR"] = os.path.join(HERE, ".cache", "matplotlib")
    e["PYTHONDONTWRITEBYTECODE"] = "1"
    return e


def _run(cmd, tag, timeout=30):
    t0 = time.time()
    proc = subprocess.run(
        cmd, cwd=HERE, env=_env(),
        capture_output=True, text=True, timeout=timeout,
    )
    elapsed = time.time() - t0
    entry = {
        "tag": tag,
        "cmd": " ".join(cmd),
        "returncode": proc.returncode,
        "elapsed": elapsed,
        "timeout": timeout,
    }
    with open(os.path.join(RESULTS, f"final_{tag}_stdout.txt"), "w") as f:
        f.write(proc.stdout)
    with open(os.path.join(RESULTS, f"final_{tag}_stderr.txt"), "w") as f:
        f.write(proc.stderr)
    return entry


def main():
    os.makedirs(RESULTS, exist_ok=True)
    os.makedirs(os.path.join(HERE, ".cache", "matplotlib"), exist_ok=True)
    py = sys.executable
    summary = []

    steps = [
        ([py, "-B", "run_experiments.py"], "driver"),
        ([py, "-B", "-m", "pytest", "-q", "test_experiments.py",
          "test_finalization.py"], "pytest"),
        ([py, "-B", "audit_grid_capacity.py"], "grid_capacity"),
    ]
    for cmd, tag in steps:
        entry = _run(cmd, tag)
        summary.append(entry)
        if entry["returncode"] != 0:
            break

    with open(os.path.join(RESULTS, "final_check_summary.json"), "w") as f:
        json.dump({"steps": summary, "final_returncode": summary[-1]["returncode"]},
                  f, indent=2)

    for e in summary:
        print(f"[{e['tag']}] rc={e['returncode']} elapsed={e['elapsed']:.2f}s")
        if e["returncode"] != 0:
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
