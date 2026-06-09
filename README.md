# A248044 — Computational Verification of Sun's Conjecture 4.3(i)

Repository for the manuscript:

**Computational Verification of Sun's Conjecture 4.3(i): Exceptional Structure of OEIS Sequence A248044**  
Carlo Corti, independent researcher  
March 2026

- GitHub repository: <https://github.com/carcorti/A248044>
- Zenodo DOI: <https://doi.org/10.5281/zenodo.19182594>
- OEIS sequence: <https://oeis.org/A248044>

---

## Overview

OEIS A248044 is the sequence defined by

\[
a(m) = \min\{ n \in \mathbb{Z}^{+} : (m+n) \mid \pi(m)^2 + \pi(n)^2 \},
\]

where \(\pi(x)\) is the prime-counting function.

This repository contains code, data, selected campaign and validation logs, and the accompanying paper for a large-scale computational verification of Sun's Conjecture 4.3(i), which asserts that such an \(n\) exists for every positive integer \(m\).

---

## Main results

| Item | Value |
|---|---:|
| Range investigated | \(m = 1, \ldots, 100000\) |
| Resolved cases | 99998 / 100000 |
| Final search bound | \(n \le 10^{13}\) |
| Unresolved cases | \(m = 19623\), \(m = 19624\) |
| Lower bound for unresolved cases | \(a(m) > 10^{13}\) |
| Largest value found | \(a(83115) = 4629736117663\) |
| Largest ratio found | \(a(40963)/40963 \approx 6.47 \times 10^7\) |
| Median value among resolved cases | \(52846\) |

The two unresolved cases have \(\pi(m)=2225=5^2\times 89\). Their prime-counting value has only prime factors congruent to \(1 \pmod 4\), so the implemented obstruction filter operates at full scope for these cases.

---

## Repository structure

```text
A248044/
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── paper/
│   ├── A248044.tex
│   └── A248044.pdf
├── src/
│   ├── sun_A248044_v5.c
│   └── sun_A248044_targeted_v3.c
├── data/
│   └── ...
├── logs/
│   └── ...
├── figures/
│   └── a248044.png
└── verification/
    ├── verify_A248044.py
    └── primecount_verification_A248044.md
```

The paper source is `paper/A248044.tex`. The compiled manuscript, if included in the repository, should be `paper/A248044.pdf`.

The files in `data/` contain the resolved values, campaign outputs, and the machine-readable representation of the two unresolved values. The files in `logs/` are selected campaign and validation logs supporting the computations. Do not interpret the repository as containing every temporary or intermediate execution artifact unless explicitly stated.

---

## Computational method

The computation uses two C programs:

- `src/sun_A248044_v5.c` — monolithic search for all \(m \le 100000\) up to \(n \le 6 \times 10^9\);
- `src/sun_A248044_targeted_v3.c` — targeted segmented search for the remaining unresolved cases up to \(n \le 10^{13}\).

Both programs use a sieve of Eratosthenes for prime enumeration and an obstruction filter based on the quadratic non-residuacity of \(-1\) modulo primes congruent to \(3 \pmod 4\). The targeted program performs segmented computation and is used for the extended ranges.

---

## Computational campaigns

| Campaign | Search range | Newly resolved | Still open | Wall time |
|---|---:|---:|---:|---:|
| M1, monolithic | \(n \le 6 \times 10^9\) | 99935 | 65 | about 10 min |
| Ext-2 | \(6 \times 10^9 \to 10^{10}\) | 36 | 29 | 137 s |
| Ext-3 | \(10^{10} \to 10^{11}\) | 20 | 9 | 1132 s |
| Ext-4 | \(10^{11} \to 10^{12}\) | 4 | 5 | 9452 s |
| Ext-5 | \(10^{12} \to 10^{13}\) | 3 | 2 | 86098 s |

The extended campaigns required about 26.9 hours of wall-clock time. Including the monolithic phase, the full computation required about 27 hours of wall-clock time, or about 432 CPU-hours on 16 threads.

---

## Hardware and software environment

All computations reported in the paper were performed on:

- CPU: AMD Ryzen 9 7940HS, 16 threads;
- RAM: 64 GB DDR5;
- operating system: Linux Mint 22.3;
- compiler: GCC with `-O3 -march=znver4 -fopenmp`.

The monolithic program requires approximately 30 GB of RAM for \(N_{\max}=6\times 10^9\). The targeted segmented program uses about 5 GB per segment.

---

## Building

Both C programs require GCC with OpenMP support.

```bash
gcc -O3 -march=native -fopenmp -o sun_A248044_v5 src/sun_A248044_v5.c
gcc -O3 -march=native -fopenmp -o sun_A248044_targeted src/sun_A248044_targeted_v3.c
```

Example targeted run:

```bash
./sun_A248044_targeted targets.txt 6000000001 10000000000
```

The targets file should contain one unresolved \(m\)-value per line. The program writes CSV-formatted output to standard output.

---

## Verification

The divisibility condition

\[
(m+a(m)) \mid \pi(m)^2 + \pi(a(m))^2
\]

was verified for the resolved entries.

The extended-campaign witnesses from Ext-2 through Ext-5 were independently checked by recomputing \(\pi(a(m))\) with `primecount` 7.10, using the Deleglise--Rivat algorithm. All 63 recomputed values matched the CSV values exactly. The congruence condition was then checked with Python arbitrary-precision integers, with zero discrepancies.

The verification script is:

```text
verification/verify_A248044.py
```

Example use:

```bash
python3 verification/verify_A248044.py data/out_ext2.csv data/out_ext3.csv data/out_ext4.csv data/out_ext5.csv
```

Python 3.10 or later is recommended. The verification script has no external Python dependencies.

---

## Data availability

The dataset, campaign outputs, selected machine-readable verification logs, code, and paper are available from:

- GitHub: <https://github.com/carcorti/A248044>
- Zenodo: <https://doi.org/10.5281/zenodo.19182594>
- OEIS: <https://oeis.org/A248044>

The paper cites the GitHub repository and the permanent Zenodo archive as the public sources for code and data.

---

## Open problems

The two values

```text
m = 19623
m = 19624
```

remain unresolved after the search up to \(n \le 10^{13}\). Equivalently, the computation establishes

\[
a(19623) > 10^{13}, \qquad a(19624) > 10^{13}.
\]

The paper gives heuristic estimates for their eventual resolution beyond \(10^{13}\), but these estimates are probabilistic and are not proofs.

---

## Citation

If you use this repository, please cite both the paper and the archived repository:

```text
Carlo Corti.
Computational Verification of Sun's Conjecture 4.3(i):
Exceptional Structure of OEIS Sequence A248044.
GitHub repository and Zenodo archive, 2026.
DOI: 10.5281/zenodo.19182594
```

Machine-readable citation metadata are provided in `CITATION.cff`.

---

## License

This repository is released under the MIT License. See `LICENSE`.

---

## AI use disclosure

The paper states that the computational implementation, including the C search programs and segmented sieve, was developed with extensive assistance from Anthropic's Claude language models (`claude-sonnet-4-6` and `claude-opus-4-6`). The author verified and approved the final text and numerical results and bears full scientific responsibility for the work.

---

## Contact

Carlo Corti  
<carlo.corti@outlook.com>
