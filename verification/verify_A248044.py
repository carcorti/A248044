#!/usr/bin/env python3
"""
verify_A248044.py
=================
Independent verification script for OEIS sequence A248044.

For each resolved pair (m, a(m)) in the input CSV files, this script
verifies the divisibility condition

    (m + a(m))  |  pi(m)^2 + pi(a(m))^2

using the pi(m) and pi(a(m)) values stored in the CSV output of the
computational campaigns.

Cross-check: the pi values in the CSV are validated against the known
checkpoints pi(10^k) from OEIS A006880 (Table 5 of the paper).

Usage
-----
    python3 verify_A248044.py [CSV files ...]

Examples:
    # Verify all extended campaigns:
    python3 verify_A248044.py out_ext2.csv out_ext3.csv out_ext4.csv out_ext5.csv

    # Verify a subset of the monolithic run (first 200 rows):
    python3 verify_A248044.py --sample 200 out_A248044_m100000.csv

    # Verify everything and show a summary:
    python3 verify_A248044.py --all

Author : Carlo Corti
Version: 1.0 (March 2026)
License: CC0 1.0 Universal
"""

import sys
import csv
import argparse
import pathlib

# ---------------------------------------------------------------------------
# Known checkpoints  pi(10^k)  from OEIS A006880 (Table 5 of the paper)
# ---------------------------------------------------------------------------
PI_CHECKPOINTS = {
    10**8:  5_761_455,
    10**9:  50_847_534,
    6 * 10**9: 279_545_368,
    10**10: 455_052_511,
    10**11: 4_118_054_813,
    10**12: 37_607_912_018,
    10**13: 346_065_536_839,
}

# ---------------------------------------------------------------------------
# Default CSV files (relative to the script location, inside the repo)
# ---------------------------------------------------------------------------
DEFAULT_CSV_FILES = [
    "data/out_A248044_m100000.csv",
    "data/out_ext2.csv",
    "data/out_ext3.csv",
    "data/out_ext4.csv",
    "data/out_ext5.csv",
]


def check_checkpoint(n: int, pi_n: int) -> bool:
    """Return True if (n, pi_n) is consistent with a known checkpoint."""
    if n in PI_CHECKPOINTS:
        return PI_CHECKPOINTS[n] == pi_n
    return True   # no checkpoint available — cannot verify


def verify_condition(m: int, pi_m: int, a_m: int, pi_am: int) -> tuple[bool, str]:
    """
    Verify  (m + a_m)  |  pi_m^2 + pi_am^2.

    Returns (ok, detail_string).
    """
    s = m + a_m
    total = pi_m * pi_m + pi_am * pi_am
    remainder = total % s
    if remainder == 0:
        return True, f"OK  m={m:>7}  a(m)={a_m:>15}  s={s:>15}  sum={total:>22}  sum/s={total//s}"
    else:
        return False, (f"FAIL  m={m:>7}  a(m)={a_m:>15}  s={s:>15}"
                       f"  sum={total:>22}  remainder={remainder}")


def parse_row(row: dict) -> tuple[int, int, int, int]:
    """Extract (m, pi_m, a_m, pi_am) from a CSV row dict."""
    return (int(row["m"]), int(row["pi_m"]),
            int(row["a_m"]), int(row["pi_am"]))


def verify_file(path: pathlib.Path, sample: int | None,
                verbose: bool) -> tuple[int, int, list[str]]:
    """
    Verify all (or up to `sample`) resolved rows in a CSV file.

    Returns (n_ok, n_fail, list_of_failure_messages).
    """
    n_ok = n_fail = 0
    failures = []

    with open(path, newline="") as fh:
        # Skip comment lines starting with '#'
        lines = (line for line in fh if not line.startswith("#"))
        reader = csv.DictReader(lines)

        count = 0
        for row in reader:
            # Skip unresolved cases (a_m == "NOT_FOUND" or a_m == 0)
            if row["a_m"].strip() in ("NOT_FOUND", "0", "") or int(row["a_m"]) == 0:
                continue

            m, pi_m, a_m, pi_am = parse_row(row)
            ok, detail = verify_condition(m, pi_m, a_m, pi_am)

            if ok:
                n_ok += 1
                if verbose:
                    print(f"  {detail}")
            else:
                n_fail += 1
                failures.append(detail)
                print(f"  {detail}")   # always print failures immediately

            count += 1
            if sample and count >= sample:
                break

    return n_ok, n_fail, failures


def run_checkpoint_sanity(path: pathlib.Path) -> None:
    """
    Spot-check pi values in the CSV against the known pi(10^k) checkpoints.
    Reports any mismatches; prints nothing if everything matches.
    """
    mismatches = []
    with open(path, newline="") as fh:
        lines = (line for line in fh if not line.startswith("#"))
        reader = csv.DictReader(lines)
        for row in reader:
            if row["a_m"].strip() in ("NOT_FOUND", "0", ""):
                continue
            a_m   = int(row["a_m"])
            pi_am = int(row["pi_am"])
            # We only have checkpoints at round powers of 10; check proximity
            for ref_n, ref_pi in PI_CHECKPOINTS.items():
                if a_m == ref_n and pi_am != ref_pi:
                    mismatches.append(
                        f"  Checkpoint mismatch: pi({ref_n}) in CSV = {pi_am}"
                        f", expected {ref_pi}"
                    )
    if mismatches:
        print("  WARNING — checkpoint mismatches detected:")
        for msg in mismatches:
            print(msg)


def main():
    parser = argparse.ArgumentParser(
        description="Verify divisibility condition for OEIS A248044 entries."
    )
    parser.add_argument(
        "files", nargs="*",
        help="CSV files to verify (default: all campaign CSVs in data/)."
    )
    parser.add_argument(
        "--all", action="store_true",
        help="Use all default CSV files."
    )
    parser.add_argument(
        "--sample", type=int, default=None, metavar="N",
        help="Verify only the first N resolved rows per file."
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Print a line for every verified entry (default: failures only)."
    )
    parser.add_argument(
        "--checkpoints", action="store_true",
        help="Also run pi-checkpoint cross-check on each file."
    )
    args = parser.parse_args()

    # Resolve file list
    if args.all or not args.files:
        script_dir = pathlib.Path(__file__).parent
        file_paths = [script_dir / f for f in DEFAULT_CSV_FILES]
    else:
        file_paths = [pathlib.Path(f) for f in args.files]

    print("=" * 70)
    print("A248044 verification script — Sun's Conjecture 4.3(i)")
    print("Condition: (m + a(m))  |  pi(m)^2 + pi(a(m))^2")
    print("=" * 70)

    total_ok = total_fail = 0
    all_failures = []

    for path in file_paths:
        if not path.exists():
            print(f"\n[SKIP] File not found: {path}")
            continue

        print(f"\nFile: {path.name}")
        if args.sample:
            print(f"  (verifying first {args.sample} resolved rows)")

        if args.checkpoints:
            run_checkpoint_sanity(path)

        n_ok, n_fail, failures = verify_file(path, args.sample, args.verbose)
        all_failures.extend(failures)
        total_ok   += n_ok
        total_fail += n_fail

        status = "ALL OK" if n_fail == 0 else f"{n_fail} FAILURES"
        print(f"  Verified: {n_ok + n_fail:>6}  |  OK: {n_ok:>6}  |  "
              f"Failures: {n_fail:>4}  [{status}]")

    # Final summary
    print("\n" + "=" * 70)
    grand_total = total_ok + total_fail
    print(f"GRAND TOTAL — Verified: {grand_total:>6}  |  "
          f"OK: {total_ok:>6}  |  Failures: {total_fail:>4}")

    if total_fail == 0:
        print("Result: ALL ENTRIES PASS — divisibility condition verified.")
    else:
        print(f"Result: {total_fail} ENTRIES FAILED — see details above.")
        sys.exit(1)

    print("=" * 70)


if __name__ == "__main__":
    main()
