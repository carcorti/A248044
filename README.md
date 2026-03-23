# A248044 — Computational Verification of Sun's Conjecture 4.3(i)

This repository contains the source code, data, and verification tools
accompanying the paper:

> Carlo Corti, *Computational Verification of Sun's Conjecture 4.3(i)
> and Exceptional Structure of Sequence A248044*,
> preprint, 2026.

---

## The conjecture

Zhi-Wei Sun (2013/2019) conjectured that for every positive integer *m*,
there exists a positive integer *n* such that

&nbsp;&nbsp;&nbsp;&nbsp;(*m* + *n*) divides π(*m*)² + π(*n*)²,

where π(*x*) denotes the prime-counting function.

The sequence *a*(*m*) = min{ *n* ∈ ℤ⁺ : (*m*+*n*) | π(*m*)²+π(*n*)² }
is [OEIS A248044](https://oeis.org/A248044).

---

## Main results

| Metric | Value |
|--------|-------|
| Range verified | *m* = 1, …, 100 000 |
| Cases resolved | 99 998 / 100 000 |
| Final search bound | *n* ≤ 10¹³ |
| Unresolved cases | *m* = 19 623 and *m* = 19 624 (lower bound *a*(*m*) > 10¹³) |
| Largest *a*(*m*) found | *a*(83 115) = 4 629 736 117 663 |
| Largest Cramér-type ratio | *a*(40 963)/*m* ≈ 6.47 × 10⁷ |

---

## Repository structure

```
A248044/
├── README.md
├── paper/
│   ├── A248044_rev5.tex          # LaTeX source of the paper (amsart class, first version)
│   └── A248044_rev6.tex          # LaTeX source of the paper (amsart class,current version)
├── src/
│   ├── sun_A248044_v5.c          # Monolithic search (campaigns M1 and Ext-1)
│   └── sun_A248044_targeted_v3.c # Segmented targeted search (campaigns Ext-2..5)
├── data/
│   ├── out_A248044_m100000.csv   # Results of campaign M1 (bound n ≤ 10⁹)
│   ├── out_ext2.csv              # 36 cases resolved in Ext-2
│   ├── out_ext3.csv              # 20 cases resolved in Ext-3
│   ├── out_ext4.csv              #  4 cases resolved in Ext-4
│   ├── out_ext5.csv              #  3 cases resolved in Ext-5
│   └── a248044.txt               # OEIS a-file: all 100 000 entries (? for gaps)
├── logs/
│   ├── log_m100000.txt           # Full log of campaign M1 (bound n ≤ 10⁹)
│   ├── log_ext1.txt              # Log of Ext-1 (monolithic, bound n ≤ 6 × 10⁹)
│   ├── log_ext2.txt              # Summary log of Ext-2
│   ├── log_ext3.txt              # Summary log of Ext-3
│   ├── log_ext4.txt              # Summary log of Ext-4
│   └── log_ext5.txt              # Summary log of Ext-5
├── figures/
│   └── a248044.png               # Logarithmic scatterplot of A248044(n), n = 1..100000
└── verification/
    └── verify_A248044.py         # Python verification script
```

**Notes on data files:**

- `out_A248044_m100000.csv` covers the initial monolithic search (bound *n* ≤ 10⁹),
  leaving 214 unresolved cases. Ext-1 extended this to *n* ≤ 6 × 10⁹ (resolving 149
  of those cases) but its CSV is not included separately, as its results are fully
  incorporated in `a248044.txt`.
- `out_ext2.csv` through `out_ext5.csv` are the targeted campaign outputs for the
  65 cases still open after Ext-1; they resolve 63 of those 65.
- `a248044.txt` is the authoritative OEIS a-file consolidating all campaigns,
  with `?` at *m* = 19 623 and *m* = 19 624 (lower bound *a*(*m*) > 10¹³).
- The OEIS b-file for A248044 (terms *n* = 1..10 000, authored by Chai Wah Wu)
  is not reproduced here as it is not our work; it is available at
  https://oeis.org/A248044/b248044.txt.

---

## Computational campaigns

| Campaign | Tool | Search range | Newly resolved | Still open | Wall time |
|----------|------|-------------|----------------|------------|-----------|
| M1 (monolithic) | `sun_A248044_v5.c` | *n* ≤ 10⁹ | 99 786 | 214 | ~10 min |
| Ext-1 (monolithic extended) | `sun_A248044_v5.c` | *n* ≤ 6 × 10⁹ | 149 | 65 | ~3 h |
| Ext-2 (targeted) | `sun_A248044_targeted_v3.c` | 6×10⁹ → 10¹⁰ | 36 | 29 | 137 s |
| Ext-3 (targeted) | `sun_A248044_targeted_v3.c` | 10¹⁰ → 10¹¹ | 20 | 9 | 1 132 s |
| Ext-4 (targeted) | `sun_A248044_targeted_v3.c` | 10¹¹ → 10¹² | 4 | 5 | 9 452 s |
| Ext-5 (targeted) | `sun_A248044_targeted_v3.c` | 10¹² → 10¹³ | 3 | **2** | 86 098 s |

Hardware: AMD Ryzen 9 7940HS, 64 GB DDR5 RAM, Linux Mint 22.3.  
Compiler: GCC 13.2.0, flags `-O3 -march=znver4 -fopenmp`.

---

## Building and running the C programs

Both programs require GCC with OpenMP support.

### Monolithic search (campaign M1)

```bash
gcc -O3 -march=native -fopenmp -o sun_A248044_v5 src/sun_A248044_v5.c
./sun_A248044_v5
```

> **Warning:** the monolithic search allocates approximately 30 GB of RAM
> for the global π-table (*N*_max = 6 × 10⁹). Ensure sufficient memory
> before running.

### Targeted segmented search (campaigns Ext-2..5)

```bash
gcc -O3 -march=native -fopenmp -o sun_A248044_targeted src/sun_A248044_targeted_v3.c
./sun_A248044_targeted <targets_file> <n_min> <n_max>
```

The targets file lists the unresolved *m*-values, one per line.
Example (reproducing Ext-2):

```bash
./sun_A248044_targeted targets.txt 6000000001 10000000000
```

Memory per segment: approximately 5 GB. The program writes results to
standard output in CSV format.

---

## Verifying the results

The script `verification/verify_A248044.py` checks the divisibility
condition (*m* + *a*(*m*)) | π(*m*)² + π(*a*(*m*))² for every resolved
entry in the CSV files, using the π values stored in the CSV itself.

```bash
# Verify all extended campaign results:
python3 verification/verify_A248044.py data/out_ext2.csv data/out_ext3.csv \
    data/out_ext4.csv data/out_ext5.csv

# Verify a sample of 500 rows from the monolithic run:
python3 verification/verify_A248044.py --sample 500 data/out_A248044_m100000.csv

# Verify everything with checkpoint cross-check against OEIS A006880:
python3 verification/verify_A248044.py --all --checkpoints
```

Expected output (all campaigns, no `--verbose`):

```
======================================================================
A248044 verification script — Sun's Conjecture 4.3(i)
Condition: (m + a(m))  |  pi(m)^2 + pi(a(m))^2
======================================================================

File: out_A248044_m100000.csv
  Verified:  99786  |  OK:  99786  |  Failures:    0  [ALL OK]

File: out_ext2.csv
  Verified:     36  |  OK:     36  |  Failures:    0  [ALL OK]

...

GRAND TOTAL — Verified:  99849  |  OK:  99849  |  Failures:    0
Result: ALL ENTRIES PASS — divisibility condition verified.
======================================================================
```

Python 3.10 or later required; no external dependencies.

---

## CSV format

All campaign CSV files share the same column structure:

| Column | Description |
|--------|-------------|
| `m` | Input value |
| `pi_m` | π(*m*) — prime-counting function at *m* |
| `a_m` | *a*(*m*) — least *n* satisfying the condition (`NOT_FOUND` if unresolved within the campaign bound) |
| `pi_am` | π(*a*(*m*)) |
| `sum` | *m* + *a*(*m*) |
| `quotient` | (π(*m*)² + π(*a*(*m*))²) / (*m* + *a*(*m*)) |
| `ratio_am_m` | *a*(*m*) / *m* — Cramér-type ratio |

Lines beginning with `#` are comments.

## OEIS a-file format

`a248044.txt` follows the standard OEIS format: one entry per line as `n a(n)`,
with `?` for the two unknown values:

```
...
19622 2676187477
19623 ?
19624 ?
19625 2676187476
...
```

---

## π checkpoints (OEIS A006880)

These values are used by the verification script and appear in Table 5
of the paper.

| *n* | π(*n*) |
|-----|--------|
| 10⁸ | 5 761 455 |
| 10⁹ | 50 847 534 |
| 6 × 10⁹ | 279 545 368 |
| 10¹⁰ | 455 052 511 |
| 10¹¹ | 4 118 054 813 |
| 10¹² | 37 607 912 018 |
| 10¹³ | 346 065 536 839 |

---

## Companion sequence

The companion sequence [OEIS A247975](https://oeis.org/A247975) encodes
Sun's Conjecture 4.1(i), which replaces π(*m*) with *p*_*m* (the *m*-th
prime). Its computational study is available at:
[github.com/carcorti/A247975](https://github.com/carcorti/A247975).

---

## Permanent archive

A permanent Zenodo archive with DOI will be linked here upon submission
of the paper.

---

## License

The source code (`src/`) and verification script (`verification/`) are
released under the [MIT License](https://opensource.org/licenses/MIT).  
The dataset (`data/`, `logs/`) is released under
[CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/).  
The paper source (`paper/`) is © Carlo Corti, all rights reserved
pending journal publication.

---

## Citation

If you use this code or data, please cite the paper:

```
@article{Corti2026A248044,
  author  = {Carlo Corti},
  title   = {Computational Verification of {Sun}'s Conjecture~4.3(i)
             and Exceptional Structure of Sequence~{A248044}},
  journal = {Mathematics of Computation},
  year    = {2026},
  note    = {Preprint available at \url{https://github.com/carcorti/A248044}}
}
```

---

## Contact

Carlo Corti — [carlo.corti@outlook.com](mailto:carlo.corti@outlook.com)
