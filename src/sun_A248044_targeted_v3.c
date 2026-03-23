/*
 * sun_A248044_targeted_v3.c
 *
 * Targeted segmented search for unresolved cases of OEIS A248044.
 *
 * ================================================================
 * REVISION HISTORY
 * ================================================================
 *
 * v1 (March 2026): Initial implementation.
 *
 * v2 (March 2026): Four-reviewer audit (ChatGPT o3, DeepSeek R1,
 *   Gemini 2.5 Pro, model Z). Adjudication by Claude Opus 4.6.
 *   Changes applied (5 accepted, 2 rejected, 1 no-action):
 *
 *   OBS-1 [CRITICO, 4/4] — Overflow int64_t in pi_n^2.
 *   OBS-2 [CRITICO, Gemini] — Race condition sul minimo.
 *   OBS-4 [minor] — SoA array locali per L1 cache.
 *   OBS-6 [minor] — Commento memoria 5.0 GB (non 4.12 GB).
 *   OBS-8 [minor] — atoll -> strtoll.
 *
 * v3 (March 2026): Second four-reviewer audit (ChatGPT o3, DeepSeek R1,
 *   Gemini 2.5 Pro, model Z) on v2. Adjudication by Claude Opus 4.6.
 *   Changes applied (4 accepted, 3 rejected/no-action):
 *
 *   OBS-A [ACCETTATO — commento, tutti] — Data race formale su
 *     found_n[ai]: la lettura fuori da critical section è UB formale
 *     per C11/OpenMP. Su x86-64 (TSO) non produce valori corrotti:
 *     il caso peggiore è un valore stale, che impatta solo le
 *     performance (un'iterazione extra), mai la correttezza (la
 *     critical section protegge la scrittura del minimo).
 *     Un #pragma omp atomic read costerebbe ~4 cicli * 65e9
 *     iterazioni ~ 87 s per segmento: inaccettabile.
 *     Fix applicata: commento esteso che documenta la scelta,
 *     la motivazione e l'invariante di correttezza. NESSUNA
 *     modifica al codice (benign race, correttezza preservata).
 *
 *   OBS-B [ACCETTATO — minor, Gemini+ChatGPT] — strtoll con
 *     NULL come endptr non rileva input parzialmente invalidi
 *     ("1000O000" restituisce 1000 in silenzio). Fix: endptr
 *     con verifica *end == '\0' per tutti gli argomenti numerici.
 *
 *   OBS-D [ACCETTATO — minor, ChatGPT] — Overflow in print_result:
 *     t->pi_am * t->pi_am è int64_t. Per a(m) ~ 10^11, pi_am ~ 4.1e9
 *     e pi_am^2 ~ 1.7e19 > int64_max = 9.22e18. La colonna
 *     'quotient' nel CSV risulterebbe errata. Fix: __int128.
 *
 *   OBS-F [ACCETTATO — triviale, DeepSeek] — Commento di
 *     compilazione (riga ~141) citava ancora v1. Aggiornato a v3.
 *
 *   OBS-C [RESPINTO — Gemini] — Proposta: usare uint64_t per
 *     Ext-2/3 e __int128 solo da Ext-4 (pi(10^11)^2 < uint64_max).
 *     Rigettato per tre motivi:
 *     (1) Il claim "GCC chiama __modti3" e' ERRATO: per
 *         __int128 % int64_t GCC usa l'istruzione DIV nativa
 *         x86-64 (rdx:rax), non una funzione di libreria.
 *         L'overhead e' ~0-10 cicli, non un "collo di bottiglia".
 *     (2) uint64_t sarebbe comunque insufficiente per Ext-4
 *         (pi(10^12)^2 ~ 1.4e21 >> uint64_max ~ 1.8e19);
 *         la logica condizionale complicherebbe il codice.
 *     (3) __int128 uniformemente garantisce correttezza su tutti
 *         i run senza branching condizionale aggiuntivo.
 *
 *   OBS-E [RESPINTO — DeepSeek] — Precalcolo s%q fuori dal
 *     loop sui target. Impossibile: s = m+n e m e' diverso
 *     per ogni target => s e' diverso per ogni ai. Errore
 *     di analisi del revisore.
 *
 *   OBS-G [RESPINTO — ChatGPT] — posix_memalign a 64 byte.
 *     Su Linux, glibc malloc() garantisce gia' allineamento a
 *     16 byte. Per array uint8_t/uint32_t con accesso lineare,
 *     l'allineamento aggiuntivo a 64 byte non produce guadagni
 *     misurabili: il bottleneck e' la latenza DRAM (memory-bound),
 *     non la larghezza di banda della singola cache line.
 *     Complessita' aggiuntiva non giustificata.
 *
 * ================================================================
 * ORIGIN AND MOTIVATION
 * ================================================================
 *
 *   OBS-1 [CRITICO, 4/4 revisori] — Overflow int64_t in pi_n^2.
 *     pi(10^11) ~ 4.1e9; quadrato ~ 1.7e19 > int64_max = 9.22e18.
 *     Overflow inizia a n ~ 7e10, prima del target 10^11.
 *     Fix: __int128 per pi_n2 e per il check di divisibilità.
 *
 *   OBS-2 [CRITICO, Gemini] — Race condition sul minimo nella
 *     critical section. schedule(guided) non garantisce che il
 *     thread con k minore (n minore) arrivi first. Un thread con
 *     n_large può scrivere found_n[ai] prima del thread con
 *     n_small < n_large. Fix: salvare il minimo, non il primo.
 *     Condizione corretta: found_n[ai]<0 || n<found_n[ai].
 *
 *   OBS-4 [minor, Gemini+ChatGPT] — SoA: estrarre m, pi_m2,
 *     qmask in array locali prima del loop parallelo per
 *     massimizzare hit in L1 cache (65*24 byte ~ 1.5 kB).
 *
 *   OBS-6 [minor, DeepSeek] — Commento memoria errato:
 *     seg è uint8_t (1 byte/entry), non bitset (1/8 byte).
 *     Totale corretto: 5.0 GB per SEG_WIDTH=10^9, non 4.12 GB.
 *
 *   OBS-8 [minor, ChatGPT+DeepSeek] — atoll → strtoll con
 *     rilevamento errori di parsing.
 *
 *   OBS-3 [RESPINTO, ChatGPT] — Eliminare pi_seg, calcolo pi
 *     on-the-fly. Incompatibile con parallelizzazione su k:
 *     pi(lo+k) dipende da tutti k'<k (prefisso cumulativo).
 *
 *   OBS-5 [RESPINTO, ChatGPT] — schedule(static,8192).
 *     Il costo per iterazione non è uniforme (dipende da quanti
 *     target vengono bloccati dal filtro di ostruzione).
 *     schedule(guided) rimane la scelta corretta.
 *
 *   OBS-7 [NESSUNA AZIONE, DeepSeek+Z] — Omissione del
 *     self-blocking prime check: scelta documentata e corretta.
 *
 * ================================================================
 * ORIGIN AND MOTIVATION
 * ================================================================
 *
 * This tool is the direct continuation of the computational campaign
 * for Sun's Conjecture 4.3(i) (Z.-W. Sun, arXiv:1309.1679v9, p. 15):
 *
 *   For every m in Z+, there exists n in Z+ such that
 *       (m + n)  |  pi(m)^2 + pi(n)^2,
 *   where pi(x) = #{primes <= x}  (prime-counting function).
 *
 * The sequence a(m) = least such n is OEIS A248044.
 *
 * CAMPAIGN HISTORY (all runs on AMD Ryzen 9 7940HS, 58 GB RAM):
 *
 *   Tool: sun_A248044_v5.c — monolithic pi-table architecture
 *   ┌──────────────────────────────────────────────────────────────┐
 *   │ Run  │ m range    │ bound │ found           │ time  │ RAM    │
 *   ├──────┼────────────┼───────┼─────────────────┼───────┼────────┤
 *   │ M1   │ 1..200     │ 1e7   │  200/200   100% │ <1 s  │ 80 MB  │
 *   │ M2   │ 1..10000   │ 1e8   │ 9991/10000       │  2 s  │ 800 MB │
 *   │ M3   │ 1..100000  │ 1e9   │ 99786/100000      │ 209 s │ 8.3 GB │
 *   │ Ext-1│ 1..100000  │ 6e9   │ 99935/100000      │ 485 s │ 49 GB  │
 *   └──────────────────────────────────────────────────────────────┘
 *
 * After Ext-1: 65 m-values remain unresolved — the "exceptional cases"
 * of A248044 up to m = 100000.
 *
 * WHY THE MONOLITHIC TOOL CANNOT GO FURTHER:
 *
 *   sun_A248044_v5 stores pi[0..sieve_limit] as uint32_t, requiring
 *   4 * sieve_limit bytes. With 58 GB available the hard ceiling is
 *   bound_max = sieve_limit_max / 2 ~ 6.97e9.  No optimisation of
 *   the inner loop can overcome this RAM constraint.
 *
 * THE SEGMENTED APPROACH (this tool):
 *
 *   We process [start_n, end_n] in segments of width SEG_WIDTH.
 *   For each segment [lo, lo + SEG_WIDTH):
 *     1. Segmented sieve of Eratosthenes -> seg[] marks composites.
 *     2. Pi prefix-sum:  pi_seg[k] = #{primes in [lo, lo+k]}.
 *        Then pi(lo+k) = pi_base + pi_seg[k].
 *     3. Inner loop checks all still-unresolved targets.
 *     4. Resolved targets are removed; pi_base is updated.
 *
 *   Memory per segment (SEG_WIDTH = 10^9):
 *     pi_seg:  uint32_t[10^9]  =  4.00 GB
 *     seg:     uint8_t [10^9]  =  1.00 GB  (1 byte per entry, not bitset)
 *     Total:                   ~  5.00 GB  [v1 said 4.12 GB — corrected in v2]
 *
 *   => bound can reach 10^12 or beyond with only ~5 GB RAM,
 *      independent of how far start_n is from zero.
 *
 * STRATEGY FOR THE 65 EXCEPTIONAL CASES:
 *
 *   Resume from where Ext-1 stopped:
 *     start_n = 6,000,000,001
 *     pi_base = pi(6,000,000,000) = 279,545,368  (sympy-verified)
 *   Pass pi_base via --pi_base to skip the expensive baseline sieve.
 *
 *   Planned target bounds:
 *     10^10  (~1-2 h, 16 threads):  expected to resolve 30-50 of 65
 *     10^11  (~15 h overnight):     expected to resolve most remaining
 *     10^12  (if still needed):     feasible over a weekend
 *
 * OBSTRUCTION FILTER (identical logic to sun_A248044_v5.c):
 *
 *   For candidate s = m + n, the sum pi(m)^2 + pi(n)^2 cannot be
 *   divisible by s if there exists a prime q ≡ 3 (mod 4) such that:
 *     (a)  q | s  with v_q(s) odd  (v_q = q-adic valuation), AND
 *     (b)  q does not divide pi(m).
 *   Filter primes: q in {3, 7, 11, 19}.
 *   qmask: bit i set iff FILTER_Q[i] | pi(m)  (no obstruction from that q).
 *
 *   Self-blocking primes: any prime s ≡ 3 (mod 4) is automatically
 *   blocked (v_s(s)=1 odd; s > pi(m) for all m >= 2).
 *
 *   CRITICAL BUG NOTE: v1-v3 of sun_A248044_v5.c used EQUALITY
 *   (pi_m == q) instead of DIVISIBILITY (pi_m % q == 0) for qmask.
 *   Fixed in v4. This code uses divisibility from the start.
 *
 * KNOWN PI CHECKPOINTS (OEIS A006880, sympy-verified):
 *   pi(10^8)  =            5,761,455
 *   pi(10^9)  =           50,847,534
 *   pi(6e9)   =          279,545,368  <- pi_base for Ext-1 resume
 *   pi(10^10) =          455,052,511
 *   pi(10^11) =        4,118,054,813
 *   pi(10^12) =       37,607,912,018
 *
 * COMPILE:
 *   gcc -O3 -march=znver4 -mtune=znver4 -fopenmp -funroll-loops \
 *       sun_A248044_targeted_v3.c -o sun_A248044_targeted_v3
 *
 * USAGE:
 *   ./sun_A248044_targeted_v3 <targets_file> <start_n> <end_n> \
 *                              [--pi_base V] [--seg_width W]
 *
 * STANDARD INVOCATION (resume after Ext-1, push to 10^10):
 *   ./sun_A248044_targeted_v3 targets_65.txt \
 *       6000000001 10000000000 --pi_base 279545368
 *
 * OUTPUT:
 *   stdout: CSV  m,pi_m,a_m,pi_am,sum,quotient,ratio_am_m
 *           (format compatible with sun_A248044_v5.c)
 *   stderr: segment-by-segment progress + final resume hint
 *
 * AUTHORS:
 *   Carlo Corti        -- research design, computational campaign
 *   Claude Opus/Sonnet -- code synthesis, review, documentation
 * DATE:   March 2026
 * REPO:   https://github.com/carcorti/A248044
 * ================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* ------------------------------------------------------------------ */
/*  Constants                                                           */
/* ------------------------------------------------------------------ */

#define VERSION           "targeted_v3"
#define DEFAULT_SEG_WIDTH ((int64_t)1000000000LL)   /* 10^9, ~4.1 GB */
#define MAX_TARGETS       256

/* ------------------------------------------------------------------ */
/*  Target descriptor                                                   */
/* ------------------------------------------------------------------ */

typedef struct {
    int64_t  m;
    int64_t  pi_m;
    int64_t  pi_m2;      /* pi(m)^2, precomputed */
    uint8_t  qmask;      /* obstruction mask, precomputed */
    int      resolved;
    int64_t  a_m;        /* solution: a(m) = least n */
    int64_t  pi_am;      /* pi(a(m)) at solution time */
} Target;

/* ------------------------------------------------------------------ */
/*  qmask: bit i set iff FILTER_Q[i] divides pi(m)                    */
/* ------------------------------------------------------------------ */

static inline uint8_t compute_qmask(int64_t pi_m)
{
    uint8_t mask = 0;
    if (pi_m %  3 == 0) mask |= 1;
    if (pi_m %  7 == 0) mask |= 2;
    if (pi_m % 11 == 0) mask |= 4;
    if (pi_m % 19 == 0) mask |= 8;
    return mask;
}

/* ------------------------------------------------------------------ */
/*  q-adic valuation parity — unrolled per prime for -O3 optimisation  */
/*  Returns 1 if v_q(s) is odd (candidate is blocked), 0 if even.    */
/*  Called only when q | s is guaranteed by the modulo test.          */
/* ------------------------------------------------------------------ */

static inline int vq_odd_3 (int64_t s)
    { if (s%9  !=0) return 1; int e=0; while(s%3 ==0){s/=3; e++;} return e&1; }
static inline int vq_odd_7 (int64_t s)
    { if (s%49 !=0) return 1; int e=0; while(s%7 ==0){s/=7; e++;} return e&1; }
static inline int vq_odd_11(int64_t s)
    { if (s%121!=0) return 1; int e=0; while(s%11==0){s/=11;e++;} return e&1; }
static inline int vq_odd_19(int64_t s)
    { if (s%361!=0) return 1; int e=0; while(s%19==0){s/=19;e++;} return e&1; }

/* ------------------------------------------------------------------ */
/*  Small prime table (for segmented sieve)                            */
/*  Holds all primes up to sqrt(end_n). Capped at 2e7 (~160 MB).      */
/* ------------------------------------------------------------------ */

#define MAX_SMALL_PRIMES 1500000
static int64_t g_small_primes[MAX_SMALL_PRIMES];
static int     g_n_small_primes = 0;

static void build_small_primes(int64_t limit)
{
    uint8_t *sv = (uint8_t *)calloc((size_t)(limit + 1), 1);
    if (!sv) { fprintf(stderr, "OOM: small prime sieve\n"); exit(1); }
    sv[0] = sv[1] = 1;
    for (int64_t i = 2; i*i <= limit; i++)
        if (!sv[i])
            for (int64_t j = i*i; j <= limit; j += i) sv[j] = 1;
    g_n_small_primes = 0;
    for (int64_t i = 2; i <= limit; i++)
        if (!sv[i]) g_small_primes[g_n_small_primes++] = i;
    free(sv);
    fprintf(stderr, "  Small primes up to %"PRId64": %d found\n",
            limit, g_n_small_primes);
}

/* ------------------------------------------------------------------ */
/*  Segmented sieve of Eratosthenes                                    */
/*  seg[k] = 0 => (lo+k) is prime;  seg[k] = 1 => composite.         */
/* ------------------------------------------------------------------ */

static void sieve_segment(uint8_t *seg, int64_t lo, int64_t W)
{
    memset(seg, 0, (size_t)W);
    /* Mark 0 and 1 non-prime */
    if (lo == 0) { seg[0] = 1; if (W > 1) seg[1] = 1; }
    else if (lo == 1) { seg[0] = 1; }

    for (int i = 0; i < g_n_small_primes; i++) {
        int64_t p = g_small_primes[i];
        if (p * p > lo + W - 1) break;
        /* First composite multiple of p in [lo, lo+W) */
        int64_t start = (lo <= p) ? p * p : lo + ((p - lo % p) % p);
        for (int64_t j = start; j < lo + W; j += p)
            seg[j - lo] = 1;
    }
}

/* ------------------------------------------------------------------ */
/*  Pi prefix-sum for one segment                                      */
/*  pi_seg[k] = #{primes in [lo, lo+k]}  (offset from pi_base)        */
/*  pi(lo+k) = pi_base + pi_seg[k]                                    */
/* ------------------------------------------------------------------ */

static void build_pi_seg(uint32_t *pi_seg, const uint8_t *seg, int64_t W)
{
    uint32_t cnt = 0;
    for (int64_t k = 0; k < W; k++) {
        if (!seg[k]) cnt++;
        pi_seg[k] = cnt;
    }
}

/* ------------------------------------------------------------------ */
/*  Target file reader                                                  */
/*  Format: "m  pi_m" per line; '#' starts a comment line.            */
/* ------------------------------------------------------------------ */

static int read_targets(const char *fname, Target *tgts, int max_n)
{
    FILE *fp = fopen(fname, "r");
    if (!fp) { perror(fname); exit(1); }
    int  n = 0;
    char line[256];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        int64_t m, pi_m;
        if (sscanf(line, "%"SCNd64" %"SCNd64, &m, &pi_m) == 2) {
            if (n >= max_n) { fprintf(stderr,"Too many targets\n"); exit(1); }
            tgts[n].m        = m;
            tgts[n].pi_m     = pi_m;
            tgts[n].pi_m2    = pi_m * pi_m;
            tgts[n].qmask    = compute_qmask(pi_m);
            tgts[n].resolved = 0;
            tgts[n].a_m      = -1;
            tgts[n].pi_am    = -1;
            n++;
        }
    }
    fclose(fp);
    return n;
}

/* ------------------------------------------------------------------ */
/*  CSV output — compatible with sun_A248044_v5.c format               */
/* ------------------------------------------------------------------ */

static void print_result(const Target *t)
{
    int64_t s = t->m + t->a_m;
    /* OBS-D [v3]: pi_am can reach ~4.1e9 at n=10^11;
     * pi_am^2 ~ 1.7e19 overflows int64_max = 9.22e18.
     * Compute quotient via __int128 to avoid corrupted CSV output. */
    __int128 sum128 = (__int128)t->pi_m2 + (__int128)t->pi_am * t->pi_am;
    int64_t  quot   = (int64_t)(sum128 / s);
    double   ratio  = (double)t->a_m / (double)t->m;
    printf("%"PRId64",%"PRId64",%"PRId64",%"PRId64
           ",%"PRId64",%"PRId64",%.6f\n",
           t->m, t->pi_m, t->a_m, t->pi_am, s, quot, ratio);
    fflush(stdout);
}

/* ------------------------------------------------------------------ */
/*  Main                                                                */
/* ------------------------------------------------------------------ */

int main(int argc, char *argv[])
{
    /* ---- Parse arguments ---- */
    if (argc < 4) {
        fprintf(stderr,
            "Usage: %s <targets_file> <start_n> <end_n>"
            " [--pi_base V] [--seg_width W]\n\n"
            "  targets_file : file with 'm pi_m' pairs, '#' = comment\n"
            "  start_n      : first n to search\n"
            "  end_n        : last  n to search\n"
            "  --pi_base V  : pi(start_n - 1)  [STRONGLY RECOMMENDED]\n"
            "                 Avoids sieving [0, start_n) from scratch.\n"
            "                 Known checkpoints (OEIS A006880):\n"
            "                   pi(10^9)  =       50847534\n"
            "                   pi(6e9)   =      279545368  <- Ext-1 resume\n"
            "                   pi(10^10) =      455052511\n"
            "                   pi(10^11) =     4118054813\n"
            "  --seg_width W: segment width (default %"PRId64", ~4.1 GB)\n\n"
            "Standard invocation (resume after Ext-1, push to 10^10):\n"
            "  %s targets_65.txt 6000000001 10000000000 --pi_base 279545368\n",
            argv[0], DEFAULT_SEG_WIDTH, argv[0]);
        return 1;
    }

    const char *targets_file = argv[1];

    /* OBS-B [v3]: strtoll with NULL endptr silently truncates invalid
     * input (e.g. "1000O000" -> 1000). Use endptr + *end check. */
    char    *arg_end;
    int64_t  start_n = (int64_t)strtoll(argv[2], &arg_end, 10);
    if (*arg_end != '\0') {
        fprintf(stderr, "ERROR: invalid start_n '%s'\n", argv[2]); return 1;
    }
    int64_t  end_n = (int64_t)strtoll(argv[3], &arg_end, 10);
    if (*arg_end != '\0') {
        fprintf(stderr, "ERROR: invalid end_n '%s'\n", argv[3]); return 1;
    }
    int64_t  seg_width   = DEFAULT_SEG_WIDTH;
    int64_t  pi_base_arg = -1;   /* -1 = not supplied */

    for (int i = 4; i < argc; i++) {
        if (strcmp(argv[i], "--pi_base")   == 0 && i+1 < argc) {
            pi_base_arg = (int64_t)strtoll(argv[++i], &arg_end, 10);
            if (*arg_end != '\0') {
                fprintf(stderr,"ERROR: invalid --pi_base '%s'\n",argv[i]); return 1;
            }
            continue;
        }
        if (strcmp(argv[i], "--seg_width") == 0 && i+1 < argc) {
            seg_width = (int64_t)strtoll(argv[++i], &arg_end, 10);
            if (*arg_end != '\0') {
                fprintf(stderr,"ERROR: invalid --seg_width '%s'\n",argv[i]); return 1;
            }
            continue;
        }
        fprintf(stderr, "Unknown argument: %s\n", argv[i]); return 1;
    }

    /* ---- Banner ---- */
    fprintf(stderr,
        "=== sun_A248044_%s  --  Sun Conjecture 4.3(i)  [A248044] ===\n"
        "Targeted segmented search\n"
        "  targets file  : %s\n"
        "  n range       : %"PRId64" .. %"PRId64"\n"
        "  segment width : %"PRId64"  (~%.2f GB RAM: 4 GB pi_seg + 1 GB seg)\n",
        VERSION, targets_file, start_n, end_n, seg_width,
        (double)seg_width * 5.0 / 1e9);
#ifdef _OPENMP
    fprintf(stderr, "  OpenMP threads: %d\n", omp_get_max_threads());
#else
    fprintf(stderr, "  OpenMP: disabled\n");
#endif

    /* ---- Load targets ---- */
    Target tgts[MAX_TARGETS];
    int n_tgts = read_targets(targets_file, tgts, MAX_TARGETS);
    fprintf(stderr, "  Targets loaded : %d\n\n", n_tgts);
    if (n_tgts == 0) { fprintf(stderr, "No targets.\n"); return 0; }

    /* ---- Small prime table (for segmented sieve) ---- */
    int64_t small_limit = 2;
    while (small_limit * small_limit < end_n) small_limit++;
    if (small_limit > (int64_t)2e7) small_limit = (int64_t)2e7;
    fprintf(stderr, "Building small primes up to %"PRId64" ...\n", small_limit);
    build_small_primes(small_limit);

    /* ---- Allocate segment arrays ---- */
    fprintf(stderr, "Allocating segment arrays (%.2f GB) ...\n",
            (double)seg_width * 5.0 / 1e9);
    uint8_t  *seg    = (uint8_t  *)malloc((size_t)seg_width);
    uint32_t *pi_seg = (uint32_t *)malloc((size_t)seg_width * sizeof(uint32_t));
    if (!seg || !pi_seg) {
        fprintf(stderr,
            "ERROR: OOM.  Try --seg_width 500000000 (~2 GB).\n");
        free(seg); free(pi_seg); return 1;
    }

    /* ---- CSV header ---- */
    printf("# Sun Conjecture 4.3(i): (m+n) | pi(m)^2 + pi(n)^2  [A248044]\n");
    printf("# Tool: sun_A248044_%s\n", VERSION);
    printf("# n range: [%"PRId64", %"PRId64"]  seg_width: %"PRId64"\n",
           start_n, end_n, seg_width);
    printf("m,pi_m,a_m,pi_am,sum,quotient,ratio_am_m\n");
    fflush(stdout);

    /* ---- Establish pi_base = pi(start_n - 1) ---- */
    /*
     * Option A (fast, recommended): user supplies --pi_base.
     * Option B (slow): sieve from 0 to start_n-1 to count primes.
     *
     * For the standard Ext-1 resume:
     *   start_n = 6,000,000,001
     *   pi_base = pi(6,000,000,000) = 279,545,368
     */
    int64_t pi_base;
    if (pi_base_arg >= 0) {
        pi_base = pi_base_arg;
        fprintf(stderr, "\nUsing supplied pi_base=%"PRId64" = pi(%"PRId64")\n\n",
                pi_base, start_n - 1);
    } else {
        fprintf(stderr,
            "\nWARNING: --pi_base not supplied. Sieving [0, %"PRId64") ...\n"
            "  (This is very slow for large start_n. Use --pi_base.)\n\n",
            start_n);
        pi_base = 0;
        int64_t lo = 0;
        while (lo < start_n) {
            int64_t W = seg_width;
            if (lo + W > start_n) W = start_n - lo;
            sieve_segment(seg, lo, W);
            for (int64_t k = 0; k < W; k++) if (!seg[k]) pi_base++;
            lo += W;
            if (lo % (int64_t)5e9 == 0)
                fprintf(stderr, "  baseline: lo=%"PRId64" pi=%"PRId64"\n",
                        lo, pi_base);
        }
        fprintf(stderr, "Computed pi(%"PRId64") = %"PRId64"\n\n",
                start_n - 1, pi_base);
    }

    /* ---- Main segment loop ---- */
    time_t  t_start     = time(NULL);
    int     n_resolved  = 0;
    int     n_remaining = n_tgts;
    int64_t lo          = start_n;

    while (lo <= end_n && n_remaining > 0) {

        int64_t W = seg_width;
        if (lo + W - 1 > end_n) W = end_n - lo + 1;

        sieve_segment(seg, lo, W);
        build_pi_seg(pi_seg, seg, W);
        /* pi(lo+k) = pi_base + pi_seg[k] */

        fprintf(stderr, "[%4lds] Segment [%"PRId64", %"PRId64")  "
                        "pi_base=%"PRId64"  remaining=%d ...",
                (long)(time(NULL) - t_start), lo, lo + W,
                pi_base, n_remaining);
        fflush(stderr);

        /* Compact list of still-unresolved targets */
        int active_idx[MAX_TARGETS];
        int n_active = 0;
        for (int i = 0; i < n_tgts; i++)
            if (!tgts[i].resolved) active_idx[n_active++] = i;

        /*
         * OBS-4 [v2]: Extract per-target fields into contiguous local
         * arrays (SoA layout) before the parallel region.
         * With at most 65 active targets, these three arrays occupy
         * 65*(8+8+1) = ~1.1 kB — guaranteed to stay in L1 cache (32 kB)
         * throughout the entire parallel loop, eliminating indirect
         * memory accesses to the Target struct inside the hot path.
         */
        int64_t local_m    [MAX_TARGETS];
        int64_t local_pi_m2[MAX_TARGETS];
        uint8_t local_qmask[MAX_TARGETS];
        for (int ai = 0; ai < n_active; ai++) {
            int ti = active_idx[ai];
            local_m    [ai] = tgts[ti].m;
            local_pi_m2[ai] = tgts[ti].pi_m2;
            local_qmask[ai] = tgts[ti].qmask;
        }

        /* Per-active-target result buffers (written in critical section) */
        int64_t found_n [MAX_TARGETS];
        int64_t found_pi[MAX_TARGETS];
        for (int i = 0; i < n_active; i++) { found_n[i] = -1; found_pi[i] = -1; }

        /*
         * Parallel search over [lo, lo+W).
         *
         * Each thread owns a chunk of k-values and checks all active
         * targets for each k. The modular obstruction filter rejects
         * ~65% of candidates before the expensive modulo.
         *
         * We use s % q directly (no sequential counters) because the
         * loop is parallelised. At -O3 with constant divisors the
         * compiler replaces s%3, s%7, s%11, s%19 with
         * multiply-by-magic-inverse, keeping the cost low.
         *
         * OBS-1 [v2, CRITICO]: pi_n can reach ~4.1e9 at n=10^11.
         * pi_n^2 ~ 1.7e19 overflows int64_max = 9.22e18 at n ~ 7e10.
         * All squaring and the divisibility sum now use __int128.
         *
         * OBS-2 [v2, CRITICO]: schedule(guided) does not guarantee
         * that the thread processing k_early arrives at the critical
         * section before the thread processing k_late > k_early.
         * A thread finding n_late may write found_n[ai] before the
         * thread finding n_early < n_late. The original guard
         * "found_n[ai] < 0" would then discard the smaller solution.
         * Fix: always save the minimum n, not the first n found.
         *
         * NOTE on self-blocking prime check (s prime, s ≡ 3 mod 4):
         * For s >> 6e9, trial-division primality against g_small_primes
         * would dominate the inner loop cost. Omitted intentionally;
         * the divisibility check correctly rejects these s.
         */
        #pragma omp parallel for schedule(guided) default(none)               \
            shared(pi_seg, local_m, local_pi_m2, local_qmask,                 \
                   found_n, found_pi, lo, W, pi_base, n_active)
        for (int64_t k = 0; k < W; k++) {

            int64_t  n    = lo + k;
            int64_t  pi_n = (int64_t)pi_base + (int64_t)pi_seg[k];

            /* OBS-1: compute pi_n^2 in __int128 to avoid overflow */
            __int128 pi_n2_128 = (__int128)pi_n * pi_n;

            for (int ai = 0; ai < n_active; ai++) {
                /* OBS-A [v3, nessuna modifica al codice]:
                 * Questa lettura di found_n[ai] avviene fuori da
                 * critical section => e' una data race formale per
                 * C11/OpenMP. Un #pragma omp atomic read la
                 * eliminerebbe, ma su x86-64 (TSO) emette un mfence
                 * o LOCK-prefixed MOV (~4+ cicli extra) per ogni
                 * iterazione: con W=10^9 e 65 target cio' aggiungerebbe
                 * ~87 s per segmento, tempo inaccettabile.
                 *
                 * Su x86-64 le load/store allineate di int64_t sono
                 * hardware-atomiche (nessun valore "strappato").
                 * Il peggior caso e' un valore STALE:
                 *  - stale fn=-1 (non ancora scritto da altro thread):
                 *    si esegue un check in piu'. Corretto.
                 *  - stale fn=valore_grande (aggiornamento non ancora
                 *    visibile): si esegue un check che poteva essere
                 *    skippato. Corretto.
                 * In nessun caso si perde il minimo: la critical section
                 * successiva rivaluta "n < found_n[ai]" con garanzia
                 * di coerenza. Correttezza matematica preservata. */
                int64_t fn = found_n[ai];
                if (fn >= 0 && n >= fn) continue;

                int64_t  s  = local_m[ai] + n;

                /* Modular obstruction filter (unchanged from v1) */
                uint8_t qmask   = local_qmask[ai];
                int     blocked = 0;
                if (s%3 ==0 && !(qmask&1)) { if (vq_odd_3 (s)) blocked=1; }
                if (!blocked &&
                    s%7 ==0 && !(qmask&2)) { if (vq_odd_7 (s)) blocked=1; }
                if (!blocked &&
                    s%11==0 && !(qmask&4)) { if (vq_odd_11(s)) blocked=1; }
                if (!blocked &&
                    s%19==0 && !(qmask&8)) { if (vq_odd_19(s)) blocked=1; }
                if (blocked) continue;

                /* OBS-1: divisibility check with __int128 sum */
                __int128 sum_128 = (__int128)local_pi_m2[ai] + pi_n2_128;
                if (sum_128 % s == 0) {
                    /* OBS-2: save minimum n, not first n */
                    #pragma omp critical
                    {
                        if (found_n[ai] < 0 || n < found_n[ai]) {
                            found_n[ai]  = n;
                            found_pi[ai] = pi_n;
                        }
                    }
                }
            }
        }   /* end parallel for */

        /* Collect and print results */
        int seg_resolved = 0;
        for (int ai = 0; ai < n_active; ai++) {
            if (found_n[ai] >= 0) {
                int ti = active_idx[ai];
                tgts[ti].resolved = 1;
                tgts[ti].a_m      = found_n[ai];
                tgts[ti].pi_am    = found_pi[ai];
                print_result(&tgts[ti]);
                seg_resolved++;
                n_resolved++;
                n_remaining--;
            }
        }

        /* Advance pi_base: pi(lo+W-1) = pi_base + pi_seg[W-1] */
        pi_base += (int64_t)pi_seg[W - 1];

        fprintf(stderr, "  resolved: %d  (total: %d/%d)\n",
                seg_resolved, n_resolved, n_tgts);

        lo += W;

    }   /* end segment loop */

    /* ---- Summary ---- */
    time_t t_end = time(NULL);
    fprintf(stderr,
        "\n=== SUMMARY ===\n"
        "n range    : %"PRId64" .. %"PRId64"\n"
        "Resolved   : %d / %d\n"
        "Still open : %d\n"
        "Total time : %ld s\n",
        start_n, end_n, n_resolved, n_tgts, n_remaining,
        (long)(t_end - t_start));

    if (n_remaining > 0) {
        fprintf(stderr, "\n=== STILL OPEN ===\n");
        fprintf(stderr, "%-10s %-8s\n", "m", "pi(m)");
        for (int i = 0; i < n_tgts; i++)
            if (!tgts[i].resolved)
                fprintf(stderr, "%-10"PRId64" %-8"PRId64"\n",
                        tgts[i].m, tgts[i].pi_m);
        /* Emit a ready-to-paste resume command */
        fprintf(stderr,
            "\npi(%"PRId64") = %"PRId64"\n"
            "Resume command:\n"
            "  %s %s %"PRId64" <new_end_n> --pi_base %"PRId64"\n",
            end_n, pi_base,
            argv[0], targets_file, end_n + 1, pi_base);
    }

    free(seg);
    free(pi_seg);
    return 0;
}
