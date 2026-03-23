/*
 * sun_A248044_v5.c  --  Computational verification of Sun's Conjecture 4.3(i)
 *
 * Conjecture (Z.-W. Sun, 2014-09-30, arXiv:1309.1679, p.15):
 *   For every m in Z+, there exists n in Z+ such that
 *     (m + n) | pi(m)^2 + pi(n)^2,
 *   where pi(x) = #{primes <= x} is the prime-counting function.
 *
 * OEIS reference: A248044.
 * Sun's worked example: pi(5)^2 + pi(12)^2 = 9 + 25 = 34 = 2 * 17,
 *   so a(5) <= 12  (and equals 12 by minimality).
 *
 * Companion sequence A248052 (second part of Conjecture 4.3):
 *   (m+n) | pi(m^2) + pi(n^2)
 *   NOT computed here; reserved for a future extension.
 *
 * -----------------------------------------------------------------
 * Key differences from sun41_v6b.c (A247975):
 *
 *  1. The role of the arithmetic function is swapped:
 *       A247975: f(k) = p_k  (k-th prime)   -- lookup in prime table
 *       A248044: f(k) = pi(k) (prime count)  -- lookup in pi table
 *
 *  2. pi(k) grows as k/ln(k), far slower than p_k ~ k*ln(k).
 *     Consequence: pi(m)^2 + pi(n)^2 fits in int64_t for ALL n
 *     within any practical search bound -- no __int128 needed.
 *     For n = 10^9: pi(n) ~ 5.1e7, pi(n)^2 ~ 2.6e15 < 2^63. OK.
 *
 *  3. The sieve is now used ONLY to build the pi table (prefix sum
 *     of the primality indicator).  No prime table is stored.
 *     RAM usage: pi table = 4 bytes * N  (uint32_t sufficient for
 *     N <= 4.3e9 since pi(4.3e9) < 2^32).
 *
 *  4. Modular obstruction filter: the A247975 filter was based on
 *     pi(m)^2 = p_m^2.  For A248044 the analogous condition is:
 *     (m+n) | pi(m)^2 + pi(n)^2.  The same algebraic argument
 *     applies: if a prime q | (m+n), q ≡ 3 (mod 4), v_q(m+n) odd,
 *     and q does not divide pi(m), then pi(n)^2 ≡ -pi(m)^2 (mod q)
 *     has no solution.  The filter code is therefore identical in
 *     structure; only pm is now pi(m) instead of p_m.
 *     IMPORTANT DIFFERENCE: pi(m) can be 0 (for m=1, since pi(1)=0)
 *     or can share a common factor with q more easily than a prime p_m
 *     would.  The filter handles this correctly: the condition
 *     q | pi(m) is checked before declaring an obstruction.
 *
 *  5. pi(k) is NOT strictly increasing (it is constant on prime gaps).
 *     Multiple consecutive values of m may share the same pi(m), which
 *     creates "plateau" structures in a(m).  This is the main new
 *     mathematical phenomenon relative to A247975.
 *
 * -----------------------------------------------------------------
 * Compile (Ryzen 9 7940HS / znver4):
 *   gcc -O3 -march=znver4 -mtune=znver4 -fopenmp -funroll-loops \
 *       sun_A248044_v5.c -o sun_A248044_v5 -lm
 *   Note: -flto and -fomit-frame-pointer omitted; with OpenMP,
 *   LTO link-time overhead can negate gains (benchmarked v5 cycle).
 *
 * Generic (any x86-64):
 *   gcc -O3 -march=native -fopenmp -funroll-loops \
 *       sun_A248044_v5.c -o sun_A248044_v5 -lm
 *
 * CHANGES v4 (corrections from three-reviewer code audit):
 *
 *   BUG FIX [CRITICAL] — compute_pm_qmask divisibility (DeepSeek review):
 *     v1-v3 set bit i iff pi(m) == FILTER_Q[i] (equality).
 *     Correct condition: bit i set iff FILTER_Q[i] DIVIDES pi(m).
 *     Mathematical basis: the obstruction requires q does not divide pi(m).
 *     When q | pi(m), pi(n)^2 ≡ -pi(m)^2 ≡ 0 (mod q) IS satisfiable
 *     (take any n with q | pi(n)), so no obstruction exists.
 *     With equality mask, any m where pi(m) is a proper multiple of
 *     q (e.g. pi(m)=6 and q=3, or pi(m)=9 and q=3) had the filter
 *     incorrectly blocking valid candidates, producing wrong a(m) values.
 *     Verified empirically: brute-force vs equality-mask vs div-mask on
 *     m=1..5000 shows dozens of wrong values with equality mask, zero
 *     discrepancies with divisibility mask.
 *     Notable examples corrected: m=13..15 (pi=6), m=23..28 (pi=9).
 *
 *   CLEANUP [Grok review]:
 *     Removed dead code: modular_obstruction_filter(), vq_is_odd(),
 *     vq_is_odd_known_divisible() — all unreachable since v3 unrolling.
 *     Removed unused #include <math.h>.
 *
 *   HARDENING [Model Z review]:
 *     Added explicit WARNING when effective search bound is reduced
 *     because sieve_limit < bound (previously silent).
 *
 * CHANGES v5 (four-reviewer code audit on v4):
 *
 *   SIMPLIFICATION [Gemini / DeepSeek / Model Z — convergent observation]:
 *     The additional check for prime s ≡ 3 (mod 4) in find_min_n() had
 *     two redundant guards:
 *       (uint32_t)s > pi_m  AND  pi_m % (uint32_t)s != 0
 *     Mathematical proof of redundancy:
 *       In the general path (m >= 2), n >= 1, so s = m+n >= m+1.
 *       By pi(m) <= m-1 for m >= 2, we have s >= m+1 > m-1 >= pi(m).
 *       Hence s > pi(m) is always true. Since s > pi(m) > 0, pi(m)
 *       cannot be a multiple of s, so pi_m % s = pi_m != 0 always.
 *       Both conditions are thus tautologies and removed without
 *       changing correctness. The check now reads simply:
 *         if (s<=sieve_limit && (s&3)==3 && bsieve_is_prime(s)) continue;
 *     Note on bsieve_is_prime(): Gemini proposed a separate "fast"
 *     variant skipping the k<2 / k==2 / k%2==0 guards. Benchmarked:
 *     GCC -O3 already eliminates these dead branches when the caller
 *     context guarantees (s&3)==3 (i.e. s is odd and >= 3). Measured
 *     speedup: 1.00x. No separate variant added.
 *
 *   COSMETIC [all four reviewers]:
 *     Corrected banner string "v1" → "v5".
 *     Corrected Usage comment "v2" → "v5".
 *     Added #define VERSION for consistent log labelling.
 *
 *   REJECTED observations:
 *     Sieve overflow i*i (DeepSeek): on LP64 Linux, long is 64-bit;
 *       for limit <= 4.3e9 (max useful), i <= 65574 and i*i <= 4.3e9
 *       << 2^63. No overflow risk.
 *     bsieve_is_prime_fast separate function (Gemini): zero measured
 *       benefit; GCC -O3 handles it. Adds API surface with no gain.
 *     omp critical → per-thread reduction (Gemini): max_am is updated
 *       O(log M_MAX) times total; contention is negligible.
 *     Sieve loop i<=limit/i guard (DeepSeek): see overflow note above.
 *     Filter extension to q=23,31,... (Grok): valid future work,
 *       deferred until baseline runs confirm computational needs.
 *     CLI hardening atoi→strtol (Grok/Gemini): deferred, not critical
 *       for controlled research use.
 *
 * Usage:
 *   ./sun_A248044_v5 [M_MAX [SEARCH_BOUND [SIEVE_LIMIT]]]
 *
 *   M_MAX        : largest m to compute (default 10000)
 *   SEARCH_BOUND : largest n to test per m (default 10^8)
 *   SIEVE_LIMIT  : sieve up to this value (default 2*SEARCH_BOUND)
 *                  must be >= SEARCH_BOUND
 *
 * Output:
 *   stdout : CSV   (redirect to file)
 *   stderr : progress + summary
 *
 * CSV columns:
 *   m, pi_m, a_m, pi_am, sum, quotient, ratio_am_m
 *   (analogous to sun41_v6b.c columns, with pi values instead of primes)
 *
 * Author: Carlo Corti (independent researcher)
 * Date:   2026
 * Based on: sun41_v6b.c (A247975 search code)
 * License: MIT
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* ------------------------------------------------------------------ */
/*  Configuration defaults                                             */
/* ------------------------------------------------------------------ */
#define DEFAULT_M_MAX        10000
#define DEFAULT_BOUND        100000000L    /* 10^8 */
#define DEFAULT_SIEVE_MULT   2             /* sieve_limit = BOUND * 2 */
#define VERSION              "v5"

/* ------------------------------------------------------------------ */
/*  Bitset sieve                                                        */
/*  Same layout as sun41_v6b.c:                                        */
/*    odd k = 2*i+3  <-->  bit i  (i = (k-3)/2)                      */
/*    bit 0 = prime, bit 1 = composite                                */
/* ------------------------------------------------------------------ */
#define BSIEVE_IDX(k)    (((k) - 3) >> 1)
#define BSIEVE_BYTES(l)  ((((l) - 3) >> 1) / 8 + 2)

static uint8_t  *g_bsieve      = NULL;
static long      g_sieve_limit = 0;

static inline int bsieve_is_prime(long k)
{
    if (k < 2)       return 0;
    if (k == 2)      return 1;
    if (k % 2 == 0)  return 0;
    long i = BSIEVE_IDX(k);
    return !((g_bsieve[i >> 3] >> (i & 7)) & 1);
}

static inline void bsieve_set_composite(long k)
{
    long i = BSIEVE_IDX(k);
    g_bsieve[i >> 3] |= (uint8_t)(1 << (i & 7));
}

/* ------------------------------------------------------------------ */
/*  pi table                                                            */
/*                                                                      */
/*  g_pi[k] = pi(k) = #{primes p : p <= k}, for k = 0 .. g_sieve_limit
 *                                                                      */
/*  uint32_t is sufficient: pi(4.29e9) < 2^32.                        */
/*  For sieve_limit <= 2e9 (typical), max pi ~ 1e8 -- well within.   */
/*                                                                      */
/*  Memory: 4 bytes * sieve_limit.                                     */
/*    sieve_limit = 2e8 -->  ~800 MB  (use when BOUND = 1e8)           */
/*    sieve_limit = 2e7 -->   ~80 MB  (use when BOUND = 1e7)           */
/*    sieve_limit = 2e6 -->    ~8 MB  (use when BOUND = 1e6)           */
/*                                                                      */
/*  For very large sieve_limit (>= 5e8), consider reducing M_MAX or   */
/*  using a segmented pi approach (not implemented in v1).             */
/* ------------------------------------------------------------------ */
static uint32_t *g_pi          = NULL;   /* g_pi[k] = pi(k)  */

/* ------------------------------------------------------------------ */
/*  Build sieve and pi table                                            */
/* ------------------------------------------------------------------ */
static void build_sieve_and_pi(long limit)
{
    /* --- bitset sieve --- */
    size_t bsz = BSIEVE_BYTES(limit);
    g_bsieve = (uint8_t *)calloc(bsz, 1);
    if (!g_bsieve) {
        fprintf(stderr, "Memory allocation failed (bitset: %.1f MB)\n",
                (double)bsz / 1e6);
        exit(1);
    }
    fprintf(stderr, "  Bitset sieve: %.1f MB\n", (double)bsz / 1e6);

    /* Sieve odd composites */
    for (long i = 3; i * i <= limit; i += 2)
        if (bsieve_is_prime(i))
            for (long j = i * i; j <= limit; j += 2 * i)
                bsieve_set_composite(j);

    /* --- pi table (prefix sum) --- */
    /* Allocate g_pi[0..limit] */
    size_t pi_bytes = ((size_t)limit + 1) * sizeof(uint32_t);
    g_pi = (uint32_t *)malloc(pi_bytes);
    if (!g_pi) {
        fprintf(stderr,
            "Memory allocation failed (pi table: %.1f MB)\n"
            "Reduce SIEVE_LIMIT or use a machine with more RAM.\n",
            (double)pi_bytes / 1e6);
        exit(1);
    }
    fprintf(stderr, "  pi table: %.1f MB\n", (double)pi_bytes / 1e6);

    g_pi[0] = 0;
    g_pi[1] = 0;  /* pi(1) = 0 */
    uint32_t count = 0;
    for (long k = 2; k <= limit; k++) {
        if (k == 2 || bsieve_is_prime(k)) count++;
        g_pi[k] = count;
    }

    fprintf(stderr, "  pi(%ld) = %u\n", limit, g_pi[limit]);
}

/* ------------------------------------------------------------------ */
/*  Modular obstruction filter — constants and helpers                  */
/*                                                                      */
/*  Small primes q ≡ 3 (mod 4) used for obstruction detection.        */
/*  Mathematical basis: if q | s, v_q(s) odd, and q does NOT divide   */
/*  pi(m), then pi(n)^2 ≡ -pi(m)^2 (mod q) has no solution            */
/*  (since -1 is a non-residue mod q). Hence (m+n) | pi(m)^2+pi(n)^2  */
/*  is impossible for this s.                                           */
/*                                                                      */
/*  Note: "q does not divide pi(m)" is the correct condition, NOT      */
/*  "q != pi(m)". When q | pi(m) (e.g. q=3, pi(m)=6), the sum        */
/*  pi(m)^2 + pi(n)^2 ≡ 0 + pi(n)^2 (mod q) IS zero-able, so no      */
/*  obstruction. Bug in v1-v3: equality check instead of divisibility. */
/* ------------------------------------------------------------------ */
#define N_FILTER_PRIMES 4

/*
 * pm_qmask: bit i set iff FILTER_Q[i] divides pi(m).
 * When bit i is set, the filter skips prime FILTER_Q[i]
 * (because q | pi(m) means no obstruction from that q).
 *
 * v1-v3 BUG: bit i was set only if pi(m) == FILTER_Q[i] (equality).
 * Correct: bit i set iff pi(m) % FILTER_Q[i] == 0 (divisibility).
 */
static inline uint8_t compute_pm_qmask(uint32_t pi_m)
{
    uint8_t mask = 0;
    if (pi_m % 3  == 0) mask |= 1;
    if (pi_m % 7  == 0) mask |= 2;
    if (pi_m % 11 == 0) mask |= 4;
    if (pi_m % 19 == 0) mask |= 8;
    return mask;
}

/* ------------------------------------------------------------------ */
/*  vq helpers — unrolled with explicit constants (v3, OBS-2)          */
/*                                                                      */
/*  Each function is called only when the counter guarantees q | s,   */
/*  so the first check tests q^2 | s (one division, constant folded   */
/*  to multiply-by-magic-inverse at -O3).                              */
/*  Using explicit per-prime functions instead of array-indexed        */
/*  vq_is_odd_known_divisible(s, qi) lets the compiler replace all    */
/*  modulo operations with constant-divisor optimisations.            */
/* ------------------------------------------------------------------ */
static inline int vq_odd_3(int64_t s)  { if(s%9  !=0) return 1; int e=0; while(s%3 ==0){s/=3; e++;} return e&1; }
static inline int vq_odd_7(int64_t s)  { if(s%49 !=0) return 1; int e=0; while(s%7 ==0){s/=7; e++;} return e&1; }
static inline int vq_odd_11(int64_t s) { if(s%121!=0) return 1; int e=0; while(s%11==0){s/=11;e++;} return e&1; }
static inline int vq_odd_19(int64_t s) { if(s%361!=0) return 1; int e=0; while(s%19==0){s/=19;e++;} return e&1; }

/* ------------------------------------------------------------------ */
/*  Core search                                                         */
/* ------------------------------------------------------------------ */
static inline int64_t find_min_n(int m, int64_t bound)
{
    uint32_t pi_m  = g_pi[m];
    int64_t  pi_m2 = (int64_t)pi_m * pi_m;
    int64_t  lim   = (bound < g_sieve_limit) ? bound : g_sieve_limit;

    int64_t s_start = (int64_t)m + 1;
    int64_t s_end   = (int64_t)m + lim;

    /* ---- Fast-path for m=1: pi(1)=0, sum is always pi(n)^2. ---- */
    /* Solution at n=1: pi(1)^2 + pi(1)^2 = 0, divisible by 2.     */
    if (pi_m == 0) {
        for (int64_t s = s_start; s <= s_end; s++) {
            int64_t n     = s - m;
            int64_t pi_n2 = (int64_t)g_pi[n] * g_pi[n];
            if (pi_n2 % s == 0) return n;
        }
        return -1;
    }

    /* ---- General path (m >= 2, pi_m >= 1) ---- */
    uint8_t pm_qmask = compute_pm_qmask(pi_m);

    /* Unrolled counters: rem_q = s_start % q, advanced by increment. */
    int r3  = (int)(s_start % 3);
    int r7  = (int)(s_start % 7);
    int r11 = (int)(s_start % 11);
    int r19 = (int)(s_start % 19);

    for (int64_t s = s_start; s <= s_end; s++) {

        /* --- Unrolled modular obstruction filter ---
         * Each rQ == 0 means Q | s (guaranteed by counter).
         * pm_qmask bit i set means Q divides pi(m): skip that prime
         * (no obstruction when Q | pi(m)).
         * vq_odd_Q() checks v_Q(s) parity with fully constant divisors. */
        int blocked = 0;
        if (r3 ==0 && !(pm_qmask & 1)) { if (vq_odd_3(s))  blocked=1; }
        if (!blocked &&
            r7 ==0 && !(pm_qmask & 2)) { if (vq_odd_7(s))  blocked=1; }
        if (!blocked &&
            r11==0 && !(pm_qmask & 4)) { if (vq_odd_11(s)) blocked=1; }
        if (!blocked &&
            r19==0 && !(pm_qmask & 8)) { if (vq_odd_19(s)) blocked=1; }

        /* Advance counters unconditionally */
        if (++r3 ==3)  r3 =0;
        if (++r7 ==7)  r7 =0;
        if (++r11==11) r11=0;
        if (++r19==19) r19=0;

        if (blocked) continue;

        /* Additional check: prime s ≡ 3 (mod 4) acts as its own obstruction
         * witness (v_s(s) = 1, odd). The condition "s does not divide pi(m)"
         * is always satisfied here: s = m+n >= m+1 > pi(m) for all m >= 2
         * (since pi(m) <= m-1), so s > pi(m) > 0 implies pi(m) % s = pi(m) != 0.
         * Both guards from v4 were therefore tautologies and are removed. */
        if (s <= g_sieve_limit && (s & 3) == 3 && bsieve_is_prime((long)s))
            continue;

        int64_t n     = s - m;
        int64_t pi_n2 = (int64_t)g_pi[n] * g_pi[n];

        /* Single modulo: safe for bound <= 10^10 (verified). */
        if ((pi_m2 + pi_n2) % s == 0)
            return n;
    }
    return -1;
}

/* ------------------------------------------------------------------ */
/*  Heuristic partial obstruction analysis (for NOT FOUND cases only)  */
/*  Same structure as sun41_v6b.c; pm_val = pi(m).                   */
/* ------------------------------------------------------------------ */
static int has_partial_obstruction(int64_t s, uint32_t pi_m)
{
    if (s <= 0 || pi_m == 0) return 0;
    /* Trial divide by small primes (up to first 200 primes ~ 1223) */
    static const int small_primes[] = {
        3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,
        79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,
        163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,
        241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,
        337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,
        431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,
        521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,
        617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,
        719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,
        823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,
        929,937,941,947,953,967,971,977,983,991,997,1009,1013,1019,
        1021,1031,1033,1039,1049,1051,1061,1063,1069,1087,1091,1093,
        1097,1103,1109,1117,1123,1129,1151,1153,1163,1171,1181,1187,
        1193,1201,1213,1217,1223,0
    };
    int64_t tmp = s;
    for (int i = 0; small_primes[i] && (int64_t)small_primes[i]*small_primes[i] <= tmp; i++) {
        int64_t p = small_primes[i];
        if (tmp % p == 0) {
            int e = 0;
            while (tmp % p == 0) { tmp /= p; e++; }
            if (pi_m % p != 0 && p % 4 == 3 && e % 2 == 1)
                return 1;
        }
    }
    /* Residual cofactor: prime q ≡ 3 mod 4, q does not divide pi(m) */
    if (tmp > 1 && pi_m % (uint32_t)tmp != 0 && tmp % 4 == 3
            && tmp <= g_sieve_limit && bsieve_is_prime((long)tmp))
        return 1;
    return 0;
}

/* ------------------------------------------------------------------ */
/*  Main                                                                */
/* ------------------------------------------------------------------ */
int main(int argc, char *argv[])
{
    int     M_MAX = DEFAULT_M_MAX;
    int64_t BOUND = DEFAULT_BOUND;

    if (argc > 1) M_MAX         = atoi(argv[1]);
    if (argc > 2) BOUND         = atoll(argv[2]);
    if (argc > 3) g_sieve_limit = atol(argv[3]);
    else          g_sieve_limit = BOUND * DEFAULT_SIEVE_MULT;

    /* Sanity: M_MAX must fit in the pi table */
    if ((long)M_MAX > g_sieve_limit) {
        fprintf(stderr, "ERROR: M_MAX=%d > sieve_limit=%ld. Increase sieve_limit.\n",
                M_MAX, g_sieve_limit);
        return 1;
    }

    /* ---- build sieve and pi table ---- */
    time_t t_start = time(NULL);
    fprintf(stderr, "=== sun_A248044_%s  --  Sun Conjecture 4.3(i)  [A248044] ===\n", VERSION);
    fprintf(stderr, "M_MAX=%d  BOUND=%ld  SIEVE_LIMIT=%ld\n",
            M_MAX, (long)BOUND, g_sieve_limit);
    fprintf(stderr, "Building sieve and pi table up to %ld ...\n", g_sieve_limit);
    build_sieve_and_pi(g_sieve_limit);
    fprintf(stderr, "Done: %.1f s\n\n", difftime(time(NULL), t_start));

    /* Spot-check pi table against known values */
    if (g_sieve_limit >= 100) {
        fprintf(stderr, "Spot checks:\n");
        fprintf(stderr, "  pi(10)=%u (expected 4)\n",  g_pi[10]);
        fprintf(stderr, "  pi(100)=%u (expected 25)\n", g_pi[100]);
    }
    if (g_sieve_limit >= 1000)
        fprintf(stderr, "  pi(1000)=%u (expected 168)\n", g_pi[1000]);
    if (g_sieve_limit >= 1000000)
        fprintf(stderr, "  pi(10^6)=%u (expected 78498)\n", g_pi[1000000]);
    if (g_sieve_limit >= 10000000)
        fprintf(stderr, "  pi(10^7)=%u (expected 664579)\n", g_pi[10000000]);
    if (g_sieve_limit >= 100000000)
        fprintf(stderr, "  pi(10^8)=%u (expected 5761455)\n", g_pi[100000000]);
    fprintf(stderr, "\n");

    /* Verify Sun's worked example: a(5) <= 12 */
    {
        uint32_t pi5  = g_pi[5];
        uint32_t pi12 = g_pi[12];
        int64_t  lhs  = (int64_t)pi5*pi5 + (int64_t)pi12*pi12;
        fprintf(stderr, "Sun example: pi(5)=%u pi(12)=%u sum=%ld mod 17 = %ld %s\n\n",
                pi5, pi12, (long)lhs, (long)(lhs % 17),
                (lhs % 17 == 0) ? "[PASS]" : "[FAIL -- pi table error]");
    }

    /* ---- allocate result array ---- */
    int64_t *result = (int64_t *)calloc(M_MAX + 1, sizeof(int64_t));
    if (!result) { fprintf(stderr, "Out of memory\n"); return 1; }

    /* ---- effective bound warning (Model Z) ---- */
    int64_t effective_bound = (BOUND < g_sieve_limit) ? BOUND : g_sieve_limit;
    if (effective_bound < BOUND)
        fprintf(stderr,
            "WARNING: sieve_limit=%ld < bound=%ld.\n"
            "         Effective search bound reduced to %ld.\n"
            "         Rerun with sieve_limit >= %ld to search the full range.\n\n",
            g_sieve_limit, (long)BOUND, (long)effective_bound, (long)BOUND);

    /* ---- parallel search ---- */
    fprintf(stderr, "=== Searching: m = 1 .. %d, effective n-bound = %ld ===\n",
            M_MAX, (long)effective_bound);

    int found = 0, not_found = 0;
    int64_t max_am = 0, max_m = 0;

    #pragma omp parallel for schedule(guided) \
            reduction(+:found, not_found)
    for (int m = 1; m <= M_MAX; m++) {
        int64_t n = find_min_n(m, BOUND);
        result[m] = n;
        if (n >= 0) {
            found++;
            #pragma omp critical
            { if (n > max_am) { max_am = n; max_m = m; } }
        } else {
            not_found++;
        }
    }

    /* ---- CSV output ---- */
    printf("# Sun Conjecture 4.3(i): (m+n) | pi(m)^2 + pi(n)^2  [A248044]\n");
    printf("# Sieve limit: %ld   Search bound: %ld\n",
           g_sieve_limit, (long)BOUND);
    printf("# Generated: %s", ctime(&t_start));
    printf("m,pi_m,a_m,pi_am,sum,quotient,ratio_am_m\n");

    for (int m = 1; m <= M_MAX; m++) {
        int64_t n    = result[m];
        uint32_t pi_m = g_pi[m];
        if (n >= 0) {
            uint32_t pi_n  = g_pi[n];
            int64_t  sum   = (int64_t)pi_m * pi_m + (int64_t)pi_n * pi_n;
            int64_t  denom = (int64_t)m + n;
            printf("%d,%u,%ld,%u,%ld,%ld,%.6f\n",
                   m, pi_m, n, pi_n, sum, sum / denom, (double)n / m);
        } else {
            printf("%d,%u,NOT_FOUND,-1,-1,-1,-1\n", m, pi_m);
        }
    }

    /* ---- Summary to stderr ---- */
    fprintf(stderr, "\n=== SUMMARY ===\n");
    fprintf(stderr, "m range:      1 .. %d\n", M_MAX);
    fprintf(stderr, "Search bound: n <= %ld\n", (long)BOUND);
    fprintf(stderr, "Found:        %d / %d (%.4f%%)\n",
            found, M_MAX, 100.0 * found / M_MAX);
    fprintf(stderr, "Not found:    %d\n", not_found);
    if (max_m > 0)
        fprintf(stderr, "Max a(m):     %ld  at m=%ld  pi(m)=%u\n",
                max_am, max_m, g_pi[max_m]);
    fprintf(stderr, "Total time:   %.0f s\n", difftime(time(NULL), t_start));

    /* ---- Obstruction analysis for NOT FOUND cases ---- */
    if (not_found > 0) {
        fprintf(stderr, "\n=== NOT FOUND cases (heuristic obstruction analysis) ===\n");
        fprintf(stderr, "%-8s %-8s %-12s\n", "m", "pi(m)", "blocked%");
        for (int m = 1; m <= M_MAX; m++) {
            if (result[m] >= 0) continue;
            uint32_t pi_m = g_pi[m];
            int blocked = 0;
            int check_range = 10000;
            for (int n = 1; n <= check_range; n++)
                if (has_partial_obstruction((int64_t)m + n, pi_m)) blocked++;
            fprintf(stderr, "%-8d %-8u %.1f%%\n",
                    m, pi_m, 100.0 * blocked / check_range);
        }
    }

    free(result);
    free(g_pi);
    free(g_bsieve);
    return 0;
}
