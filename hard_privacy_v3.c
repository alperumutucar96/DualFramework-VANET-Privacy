/**
 * hard_privacy_v3.c  —  Hard Privacy v3, Real Network Simulation
 * Paillier-1024 HE | Schnorr ZKP | Shamir MPC | 6G 3GPP TR38.824
 * Paper equations: 11, 12, 13, 14, 15, 16
 *
 * Compile (MSYS2 MinGW-64 / Linux):
 *   gcc -O2 -o hard_v3.exe hard_privacy_v3.c crypto_primitives.c
 *       -lgmp -lssl -lcrypto -lm -lpthread
 *
 * Run:
 *   ./hard_v3.exe
 *   Output: hard_v3_results.json
 */

#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

#include <gmp.h>
#include "crypto_primitives.h"

/* ── simulation parameters ─────────────────────────────────────────── */
static const int VEHICLE_COUNTS[] = {5, 10, 25, 50, 100, 250, 500};
#define N_COUNTS  (int)(sizeof(VEHICLE_COUNTS)/sizeof(VEHICLE_COUNTS[0]))
#define K_RSUS    4
#define C_CORES   16

static const char MPC_PRIME_STR[] =
    "170141183460469231731687303715884105773";

static double now_ms(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1e6;
}

#define MAX_STEPS 12

typedef struct { char name[64]; double ms; } Step;

typedef struct {
    int    N, K, C;
    Step   steps[MAX_STEPS];
    int    n_steps;
    double total_ms;
} RunResult;

static void result_add(RunResult *r, const char *name, double ms)
{
    int i = r->n_steps++;
    snprintf(r->steps[i].name, sizeof(r->steps[i].name), "%s", name);
    r->steps[i].ms  = ms;
    r->total_ms    += ms;
}

/* ════════════════════════════════════════════════════════════════════
   Per-vehicle bundle (Step 1 output)
   ════════════════════════════════════════════════════════════════════ */
typedef struct {
    mpz_t    blind_sig;
    ZkpProof piP, piB, piS;
    mpz_t    ctP, ctB, ctS, ctI;   /* Paillier ciphertexts */
} VBundle;

/* ════════════════════════════════════════════════════════════════════
   Step 1 worker: BlindSig + ZKP.Prove × 3 + Paillier.Enc × 4
   ════════════════════════════════════════════════════════════════════ */
typedef struct {
    const RSACtx      *rsa_v;
    const PaillierCtx *he;
    Vehicle           *vehs;
    VBundle           *bundles;
    int                start, end;
} S1Arg;

static void *worker_step1(void *arg)
{
    S1Arg *a = (S1Arg *)arg;
    char msg[32];

    for (int i = a->start; i < a->end; i++) {
        Vehicle *v = &a->vehs[i];
        VBundle *b = &a->bundles[i];

        /* Eq 5: BlindSign */
        snprintf(msg, sizeof(msg), "%lld", (long long)v->D_CID);
        blind_sign(a->rsa_v, b->blind_sig,
                   (const uint8_t *)msg, strlen(msg));

        /* Eq 12: ZKP.Prove for CP, CB, CS */
        zkp_prove(&b->piP, v->D_CP);
        zkp_prove(&b->piB, v->D_CB);
        zkp_prove(&b->piS, v->D_CS);

        /* Eq 11: HE.Enc for CP, CB, CS, CID%10000 */
        mpz_t m; mpz_init(m);

        mpz_set_si(m, v->D_CP);
        paillier_enc(a->he, b->ctP, m);

        mpz_set_si(m, v->D_CB);
        paillier_enc(a->he, b->ctB, m);

        mpz_set_si(m, v->D_CS);
        paillier_enc(a->he, b->ctS, m);

        mpz_set_si(m, v->D_CID % 10000);
        paillier_enc(a->he, b->ctI, m);

        mpz_clear(m);
    }
    return NULL;
}

/* ════════════════════════════════════════════════════════════════════
   Step 3 worker: ZKP Verify  (K×C threads)
   ════════════════════════════════════════════════════════════════════ */
typedef struct {
    VBundle *bundles;
    int      start, end;
    int      all_ok;    /* output */
} S3Arg;

static void *worker_step3(void *arg)
{
    S3Arg *a = (S3Arg *)arg;
    a->all_ok = 1;
    for (int i = a->start; i < a->end; i++) {
        VBundle *b = &a->bundles[i];
        /* Eq 16: Verify(π) = true */
        if (!zkp_verify(&b->piP) ||
            !zkp_verify(&b->piB) ||
            !zkp_verify(&b->piS))
            a->all_ok = 0;
    }
    return NULL;
}

/* ════════════════════════════════════════════════════════════════════
   Step 4: parallel Paillier tree-reduction for one field
   Each level halves the list; pairs are computed in parallel threads.
   ════════════════════════════════════════════════════════════════════ */

/* One pair-reduction job */
typedef struct {
    const PaillierCtx *he;
    mpz_t             *a, *b;   /* input pair  */
    mpz_t             *out;     /* result cell */
} PairJob;

static void *worker_pair(void *arg)
{
    PairJob *j = (PairJob *)arg;
    paillier_add(j->he, *j->out, *j->a, *j->b);
    return NULL;
}

/*
 * tree_reduce: in-place binary tree reduction of cts[0..n)
 * Returns the aggregate in cts[0].  All other slots undefined after.
 * Uses up to max_par parallel threads per level.
 */
static void tree_reduce(const PaillierCtx *he,
                        mpz_t *cts, int n, int max_par)
{
    while (n > 1) {
        int pairs = n / 2;
        int nt    = pairs < max_par ? pairs : max_par;

        PairJob   *jobs = malloc(pairs * sizeof(PairJob));
        pthread_t *tids = malloc(nt    * sizeof(pthread_t));

        for (int p = 0; p < pairs; p++) {
            jobs[p].he  = he;
            jobs[p].a   = &cts[2*p];
            jobs[p].b   = &cts[2*p + 1];
            jobs[p].out = &cts[2*p];    
        }

        int launched = 0;
        while (launched < pairs) {
            int batch = pairs - launched < nt ? pairs - launched : nt;
            for (int t = 0; t < batch; t++)
                pthread_create(&tids[t], NULL, worker_pair,
                               &jobs[launched + t]);
            for (int t = 0; t < batch; t++)
                pthread_join(tids[t], NULL);
            launched += batch;
        }

        if (n & 1) {
            mpz_set(cts[pairs], cts[n-1]);
            n = pairs + 1;
        } else {
            n = pairs;
        }

        free(jobs);
        free(tids);
    }
}

/* ════════════════════════════════════════════════════════════════════
   Step 4 launcher: 3 fields in 3 parallel threads
   ════════════════════════════════════════════════════════════════════ */
typedef struct {
    const PaillierCtx *he;
    mpz_t             *cts;
    int                n;
} FieldReduceArg;

static void *worker_field_reduce(void *arg)
{
    FieldReduceArg *a = (FieldReduceArg *)arg;
    tree_reduce(a->he, a->cts, a->n, 8);
    return NULL;
}

/* ════════════════════════════════════════════════════════════════════
   Step 9 worker: Vehicle HE Decrypt (3 weights per vehicle)
   ════════════════════════════════════════════════════════════════════ */
typedef struct {
    const PaillierCtx *he;
    mpz_t             *enc_w;   /* 3-element array of encrypted weights */
    int                start, end;
} S9Arg;

static void *worker_step9(void *arg)
{
    S9Arg *a = (S9Arg *)arg;
    mpz_t tmp; mpz_init(tmp);
    for (int i = a->start; i < a->end; i++) {
        paillier_dec(a->he, tmp, a->enc_w[0]);
        paillier_dec(a->he, tmp, a->enc_w[1]);
        paillier_dec(a->he, tmp, a->enc_w[2]);
    }
    mpz_clear(tmp);
    return NULL;
}

/* ════════════════════════════════════════════════════════════════════
   run_hard_v3
   ════════════════════════════════════════════════════════════════════ */
static RunResult run_hard_v3(int N,
                              const PaillierCtx *he,
                              const RSACtx      *rsa_fe)
{
    RunResult res;
    memset(&res, 0, sizeof(res));
    res.N = N; res.K = K_RSUS; res.C = C_CORES;

    Network net = network_default();

    /* vehicle-side RSA for blind signatures */
    RSACtx rsa_v;
    rsa_init(&rsa_v);

    /* allocate vehicles + bundles */
    Vehicle *vehs    = malloc(N * sizeof(Vehicle));
    VBundle *bundles = malloc(N * sizeof(VBundle));
    gen_vehicles(vehs, N, 42);

    /* initialise all mpz fields inside bundles */
    for (int i = 0; i < N; i++) {
        mpz_init(bundles[i].blind_sig);
        mpz_inits(bundles[i].ctP,
                  bundles[i].ctB,
                  bundles[i].ctS,
                  bundles[i].ctI, NULL);
    }

    /* ── Step 1: BlindSig + ZKP + HE Enc  [≤8 threads] ─────────────
       Eq 5, 11, 12                                                   */
    int nt1 = N < 8 ? N : 8;
    S1Arg   *s1args = malloc(nt1 * sizeof(S1Arg));
    pthread_t *tids = malloc(nt1 * sizeof(pthread_t));
    int chunk1 = (N + nt1 - 1) / nt1;

    for (int t = 0; t < nt1; t++) {
        s1args[t].rsa_v   = &rsa_v;
        s1args[t].he      = he;
        s1args[t].vehs    = vehs;
        s1args[t].bundles = bundles;
        s1args[t].start   = t * chunk1;
        s1args[t].end     = (t+1)*chunk1 < N ? (t+1)*chunk1 : N;
    }

    double t0 = now_ms();
    for (int t = 0; t < nt1; t++)
        pthread_create(&tids[t], NULL, worker_step1, &s1args[t]);
    for (int t = 0; t < nt1; t++)
        pthread_join(tids[t], NULL);
    result_add(&res, "Step 1: BlindSig+ZKP+HE Enc (parallel)",
               now_ms() - t0);
    free(s1args);

    /* ── Step 2: 6G OFDMA uplink ─────────────────────────────────── */
    result_add(&res, "Step 2: V->RSU 6G uplink",
               network_delay(&net, net.d_v2r,
                             B_RSA + 3*B_ZKP + 4*B_PAIL));

    /* ── Step 3: RSU ZKP Verify  [K×C threads]
       Eq 16: Verify(πi) = true                                       */
    int n_workers3 = K_RSUS * C_CORES;
    if (n_workers3 > N) n_workers3 = N;
    S3Arg *s3args = malloc(n_workers3 * sizeof(S3Arg));
    int chunk3 = (N + n_workers3 - 1) / n_workers3;

    /* reuse tids — may need more space */
    free(tids);
    tids = malloc(n_workers3 * sizeof(pthread_t));

    t0 = now_ms();
    for (int t = 0; t < n_workers3; t++) {
        s3args[t].bundles = bundles;
        s3args[t].start   = t * chunk3;
        s3args[t].end     = (t+1)*chunk3 < N ? (t+1)*chunk3 : N;
        pthread_create(&tids[t], NULL, worker_step3, &s3args[t]);
    }
    for (int t = 0; t < n_workers3; t++)
        pthread_join(tids[t], NULL);
    result_add(&res, "Step 3: RSU parallel ZKP Verify (KxC)",
               now_ms() - t0);
    free(s3args);
    free(tids);

    /* ── Step 4: RSU HE Aggregation — binary tree-reduction
       Three fields reduced in parallel threads.                      */

    /* Copy pointers into flat arrays (tree_reduce works in-place) */
    mpz_t *arP = malloc(N * sizeof(mpz_t));
    mpz_t *arB = malloc(N * sizeof(mpz_t));
    mpz_t *arS = malloc(N * sizeof(mpz_t));

    for (int i = 0; i < N; i++) {
        mpz_init_set(arP[i], bundles[i].ctP);
        mpz_init_set(arB[i], bundles[i].ctB);
        mpz_init_set(arS[i], bundles[i].ctS);
    }

    FieldReduceArg fra[3] = {
        { he, arP, N },
        { he, arB, N },
        { he, arS, N }
    };
    pthread_t ftids[3];

    t0 = now_ms();
    for (int j = 0; j < 3; j++)
        pthread_create(&ftids[j], NULL, worker_field_reduce, &fra[j]);
    for (int j = 0; j < 3; j++)
        pthread_join(ftids[j], NULL);
    result_add(&res, "Step 4: RSU HE Agg tree-reduction (3 fields)",
               now_ms() - t0);

    /* aggregated ciphertexts are now in arP[0], arB[0], arS[0] */

    /* ── Step 5: RSU MPC — Shamir SSS + 2-round RTT model
       Eq 14: (C_CP_MPC,...) = MPC(C_CP_agg,...)                      */
    t0 = now_ms();

    mpz_t mpc_prime; mpz_init_set_str(mpc_prime, MPC_PRIME_STR, 10);

    mpz_t vP, vB, vS;
    mpz_inits(vP, vB, vS, NULL);
    mpz_mod(vP, arP[0], mpc_prime);
    mpz_mod(vB, arB[0], mpc_prime);
    mpz_mod(vS, arS[0], mpc_prime);

    Share shP[3], shB[3], shS[3];
    mpc_share(shP, vP);
    mpc_share(shB, vB);
    mpc_share(shS, vS);

    mpz_t rec; mpz_init(rec);
    mpc_reconstruct(rec, shP);
    mpc_reconstruct(rec, shB);
    mpc_reconstruct(rec, shS);
    double crypto_ms = now_ms() - t0;

    /* 2 online MPC rounds: each round = RSU↔BS↔CA (full round-trip) */
    double rtt_ms = 0.0;
    for (int r = 0; r < 2; r++) {
        rtt_ms += network_delay(&net, net.d_r2b, B_PAIL*3)
               +  network_delay(&net, net.d_b2c, B_PAIL*3)
               +  network_delay(&net, net.d_b2c, B_PAIL*3)
               +  network_delay(&net, net.d_r2b, B_PAIL*3);
    }
    result_add(&res, "Step 5: RSU MPC (Shamir+2-round RTT)",
               crypto_ms + rtt_ms);

    for (int i = 0; i < 3; i++) {
        share_free(&shP[i]); share_free(&shB[i]); share_free(&shS[i]);
    }

    /* ── Step 6: RSU → BS → CA DTLS + links ─────────────────────── */
    result_add(&res, "Step 6: RSU->BS->CA DTLS+links",
               0.5
               + network_delay(&net, net.d_r2b, B_PAIL*3)
               + network_delay(&net, net.d_b2c, B_PAIL*3));

    /* ── Step 7: CA HE ML Update
       Eq 15: CT_w_new = Eval(PK_HE, f_ML, CT_w_old, C_CP_MPC, ...)
       Proxy: dec aggregates → gradient step → re-enc weights         */
    t0 = now_ms();

    mpz_t dec_P, dec_B, dec_S;
    mpz_inits(dec_P, dec_B, dec_S, NULL);
    paillier_dec(he, dec_P, arP[0]);
    paillier_dec(he, dec_B, arB[0]);
    paillier_dec(he, dec_S, arS[0]);

    /* plaintext aggregates (bounded for gradient) */
    double vPf = (double)(mpz_get_ui(dec_P) % 10000);
    double vBf = (double)(mpz_get_ui(dec_B) % 10);
    double vSf = (double)(mpz_get_ui(dec_S) % 10);

    double w_old[3] = {0.4, 0.3, 0.3};
    double X[3]     = {vPf, vBf, vSf};
    double norm = sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]) + 1e-9;
    for (int j = 0; j < 3; j++) X[j] /= norm;
    double dot = w_old[0]*X[0] + w_old[1]*X[1] + w_old[2]*X[2];
    double w_new[3];
    for (int j = 0; j < 3; j++)
        w_new[j] = w_old[j] - 0.01*(dot - 0.5)*X[j];

    /* re-encrypt updated weights (3 in parallel) */
    mpz_t enc_w[3];
    mpz_inits(enc_w[0], enc_w[1], enc_w[2], NULL);
    {
        mpz_t wm; mpz_init(wm);
        for (int j = 0; j < 3; j++) {
            int64_t wval = (int64_t)(fabs(w_new[j]) * 1e4);
            mpz_set_si(wm, wval);
            paillier_enc(he, enc_w[j], wm);
        }
        mpz_clear(wm);
    }
    result_add(&res, "Step 7: CA HE ML Update (dec+grad+re-enc)",
               now_ms() - t0);

    /* ── Step 8: CA Predict + distribution ──────────────────────── */
    t0 = now_ms();
    for (int i = 0; i < N; i++) {
        double cp = (double)vehs[i].D_CP / 999.0;
        double cb = (double)vehs[i].D_CB;
        double cs = (double)vehs[i].D_CS / 3.0;
        (void)(w_new[0]*cp + w_new[1]*cb + w_new[2]*cs);
    }
    result_add(&res, "Step 8: CA Predict+CA->BS->RSU dist",
               (now_ms() - t0)
               + network_delay(&net, net.d_b2c, B_PAIL)
               + network_delay(&net, net.d_r2b, B_PAIL));

    /* ── Step 9: Vehicle HE Decrypt  [≤8 threads, 3 weights each]
       Y_dec,i = HE.Dec(sk_i, HE_Y_enc,i) = {Strategy_CP, CB, CS}   */
    int nt9 = N < 8 ? N : 8;
    S9Arg     *s9args = malloc(nt9 * sizeof(S9Arg));
    pthread_t *tids9  = malloc(nt9 * sizeof(pthread_t));
    int chunk9 = (N + nt9 - 1) / nt9;

    for (int t = 0; t < nt9; t++) {
        s9args[t].he    = he;
        s9args[t].enc_w = enc_w;
        s9args[t].start = t * chunk9;
        s9args[t].end   = (t+1)*chunk9 < N ? (t+1)*chunk9 : N;
    }
    t0 = now_ms();
    for (int t = 0; t < nt9; t++)
        pthread_create(&tids9[t], NULL, worker_step9, &s9args[t]);
    for (int t = 0; t < nt9; t++)
        pthread_join(tids9[t], NULL);
    result_add(&res, "Step 9: Vehicle HE Decrypt (3 parallel)",
               now_ms() - t0);
    free(s9args); free(tids9);

    /* ── Step 10: ACK ────────────────────────────────────────────── */
    result_add(&res, "Step 10: ACK",
               network_delay(&net, net.d_v2r, 64));

    mpz_clears(enc_w[0], enc_w[1], enc_w[2], NULL);
    mpz_clears(dec_P, dec_B, dec_S, NULL);
    mpz_clears(vP, vB, vS, rec, mpc_prime, NULL);

    for (int i = 0; i < N; i++) {
        mpz_clear(arP[i]);
        mpz_clear(arB[i]);
        mpz_clear(arS[i]);
    }
    free(arP); free(arB); free(arS);

    for (int i = 0; i < N; i++) {
        mpz_clear(bundles[i].blind_sig);
        mpz_clears(bundles[i].ctP, bundles[i].ctB,
                   bundles[i].ctS, bundles[i].ctI, NULL);
        zkp_free(&bundles[i].piP);
        zkp_free(&bundles[i].piB);
        zkp_free(&bundles[i].piS);
    }
    free(bundles);
    free(vehs);
    rsa_free(&rsa_v);

    return res;
}

/* ════════════════════════════════════════════════════════════════════
   JSON output
   ════════════════════════════════════════════════════════════════════ */
static void write_json(const char *path,
                       RunResult *results, int count)
{
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "Cannot open %s\n", path); return; }

    fprintf(f, "[\n");
    for (int i = 0; i < count; i++) {
        RunResult *r = &results[i];
        fprintf(f, "  {\n");
        fprintf(f, "    \"N\": %d,\n", r->N);
        fprintf(f, "    \"K\": %d,\n", r->K);
        fprintf(f, "    \"C\": %d,\n", r->C);
        fprintf(f, "    \"steps\": {\n");
        for (int s = 0; s < r->n_steps; s++)
            fprintf(f, "      \"%s\": %.4f%s\n",
                    r->steps[s].name, r->steps[s].ms,
                    s < r->n_steps-1 ? "," : "");
        fprintf(f, "    },\n");
        fprintf(f, "    \"total_ms\": %.4f\n", r->total_ms);
        fprintf(f, "  }%s\n", i < count-1 ? "," : "");
    }
    fprintf(f, "]\n");
    fclose(f);
}

/* ════════════════════════════════════════════════════════════════════
   main
   ════════════════════════════════════════════════════════════════════ */
int main(void)
{
    printf("==========================================================================\n");
    printf("  HARD PRIVACY v3 — REAL NETWORK SIMULATION\n");
    printf("  Paillier-1024 HE | Schnorr ZKP | Shamir MPC | 6G 3GPP TR38.824\n");
    printf("  K=%d RSUs, C=%d cores/RSU\n", K_RSUS, C_CORES);
    printf("==========================================================================\n");

    printf("  [keygen] Generating Paillier + RSA keys...\n");
    fflush(stdout);

    PaillierCtx he;
    RSACtx      rsa_fe;
    paillier_init(&he);
    rsa_init(&rsa_fe);
    printf("  [keygen] Done.\n");
    fflush(stdout);

    RunResult all_r[N_COUNTS];

    for (int ni = 0; ni < N_COUNTS; ni++) {
        int N = VEHICLE_COUNTS[ni];
        printf("\n-- N=%d --\n", N);
        fflush(stdout);

        all_r[ni] = run_hard_v3(N, &he, &rsa_fe);

        for (int s = 0; s < all_r[ni].n_steps; s++)
            printf("  %-50s %8.2f ms\n",
                   all_r[ni].steps[s].name,
                   all_r[ni].steps[s].ms);
        printf("  %-50s %8.2f ms\n", "Total E2E", all_r[ni].total_ms);
        fflush(stdout);
    }

    write_json("hard_v3_results.json", all_r, N_COUNTS);
    printf("\n[OK] Saved hard_v3_results.json\n");

    paillier_free(&he);
    rsa_free(&rsa_fe);
    return 0;
}
