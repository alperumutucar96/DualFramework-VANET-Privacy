
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

#include "crypto_primitives.h"

/* ── simulation parameters ─────────────────────────────────────────── */
static const int VEHICLE_COUNTS[] = {5, 10, 25, 50, 100, 250, 500};
#define N_COUNTS  (int)(sizeof(VEHICLE_COUNTS)/sizeof(VEHICLE_COUNTS[0]))
#define K_RSUS    4     /* number of RSUs                    */
#define C_CORES   16    /* parallel cores per RSU            */
#define MAX_VEH   500

/* ── timing helper ─────────────────────────────────────────────────── */
static double now_ms(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1e6;
}

/* ── result structures ─────────────────────────────────────────────── */
#define MAX_STEPS 10

typedef struct {
    char   name[64];
    double ms;
} Step;

typedef struct {
    int    N, K, C;
    char   loop[8];       /* "fast" or "slow" */
    Step   steps[MAX_STEPS];
    int    n_steps;
    double total_ms;
} LoopResult;

static void result_add(LoopResult *r, const char *name, double ms)
{
    int i = r->n_steps++;
    snprintf(r->steps[i].name, sizeof(r->steps[i].name), "%s", name);
    r->steps[i].ms = ms;
    r->total_ms   += ms;
}

/* ════════════════════════════════════════════════════════════════════
   Thread-pool worker for Step 1: vehicle BlindSig + ABE Enc
   ════════════════════════════════════════════════════════════════════ */
typedef struct {
    const RSACtx *rsa_auth;
    Vehicle      *vehs;
    int           start, end;   /* slice [start, end) */
    /* outputs */
    mpz_t        *blind_sigs;   /* pre-allocated array */
    AbeBundle    *bundles_cp;
    AbeBundle    *bundles_cb;
    AbeBundle    *bundles_cs;
} VEncArg;

static void *worker_v_enc(void *arg)
{
    VEncArg *a = (VEncArg *)arg;
    char msg[32];
    for (int i = a->start; i < a->end; i++) {
        /* Eq 5: D'_CIDi = BlindSign(sk_auth, D_CIDi) */
        snprintf(msg, sizeof(msg), "%lld", (long long)a->vehs[i].D_CID);
        blind_sign(a->rsa_auth, a->blind_sigs[i],
                   (const uint8_t *)msg, strlen(msg));
        /* Eq 6: D'_CPi = ABE.Enc(PK_ABE, D_CPi, P) */
        abe_enc(a->rsa_auth, &a->bundles_cp[i], a->vehs[i].D_CP);
        abe_enc(a->rsa_auth, &a->bundles_cb[i], a->vehs[i].D_CB);
        abe_enc(a->rsa_auth, &a->bundles_cs[i], a->vehs[i].D_CS);
    }
    return NULL;
}

/* ════════════════════════════════════════════════════════════════════
   Thread-pool worker for Step 3: RSU ABE Decrypt
   ════════════════════════════════════════════════════════════════════ */
typedef struct {
    const RSACtx *rsa_auth;
    AbeBundle    *bundles_cp;
    AbeBundle    *bundles_cb;
    AbeBundle    *bundles_cs;
    int           start, end;
    int64_t      *dec_cp, *dec_cb, *dec_cs;
} VDecArg;

static void *worker_v_dec(void *arg)
{
    VDecArg *a = (VDecArg *)arg;
    for (int i = a->start; i < a->end; i++) {
        a->dec_cp[i] = abe_dec(a->rsa_auth, &a->bundles_cp[i]);
        a->dec_cb[i] = abe_dec(a->rsa_auth, &a->bundles_cb[i]);
        a->dec_cs[i] = abe_dec(a->rsa_auth, &a->bundles_cs[i]);
    }
    return NULL;
}

/* ════════════════════════════════════════════════════════════════════
   Thread-pool worker for Step 7: Vehicle FE Decrypt
   ════════════════════════════════════════════════════════════════════ */
typedef struct {
    const RSACtx *rsa_fe;
    AbeBundle    *enc_strats;
    int           start, end;
} FDecArg;

static void *worker_fe_dec(void *arg)
{
    FDecArg *a = (FDecArg *)arg;
    for (int i = a->start; i < a->end; i++)
        fe_dec(a->rsa_fe, &a->enc_strats[i]);
    return NULL;
}

/* ── generic thread-pool launcher ──────────────────────────────────── */
static void run_parallel(int N, int max_t,
                         void *(*f)(void *), void **args)
{
    int nt = N < max_t ? N : max_t;
    if (nt <= 0) nt = 1;
    pthread_t *tids = malloc(nt * sizeof(pthread_t));
    for (int t = 0; t < nt; t++)
        pthread_create(&tids[t], NULL, f, args[t]);
    for (int t = 0; t < nt; t++)
        pthread_join(tids[t], NULL);
    free(tids);
}

/* ════════════════════════════════════════════════════════════════════
   FAST LOOP
   ════════════════════════════════════════════════════════════════════ */
static LoopResult fast_loop(int N,
                             const RSACtx *rsa_auth,
                             const RSACtx *rsa_fe)
{
    LoopResult res;
    memset(&res, 0, sizeof(res));
    res.N = N; res.K = K_RSUS; res.C = C_CORES;
    strcpy(res.loop, "fast");

    Network net = network_default();

    Vehicle   *vehs       = malloc(N * sizeof(Vehicle));
    mpz_t     *bsigs      = malloc(N * sizeof(mpz_t));
    AbeBundle *bun_cp     = malloc(N * sizeof(AbeBundle));
    AbeBundle *bun_cb     = malloc(N * sizeof(AbeBundle));
    AbeBundle *bun_cs     = malloc(N * sizeof(AbeBundle));
    int64_t   *dec_cp     = malloc(N * sizeof(int64_t));
    int64_t   *dec_cb     = malloc(N * sizeof(int64_t));
    int64_t   *dec_cs     = malloc(N * sizeof(int64_t));
    double    *strats     = malloc(N * sizeof(double));
    AbeBundle *enc_strats = malloc(N * sizeof(AbeBundle));

    gen_vehicles(vehs, N, 42);
    for (int i = 0; i < N; i++) mpz_init(bsigs[i]);

    /* ── Step 1: Vehicle BlindSig + ABE Enc  [N parallel, ≤16 threads] */
    int nt1 = N < 16 ? N : 16;
    VEncArg *eargs = malloc(nt1 * sizeof(VEncArg));
    void   **eptrs = malloc(nt1 * sizeof(void *));
    int chunk1 = (N + nt1 - 1) / nt1;
    for (int t = 0; t < nt1; t++) {
        eargs[t].rsa_auth   = rsa_auth;
        eargs[t].vehs       = vehs;
        eargs[t].start      = t * chunk1;
        eargs[t].end        = (t+1)*chunk1 < N ? (t+1)*chunk1 : N;
        eargs[t].blind_sigs = bsigs;
        eargs[t].bundles_cp = bun_cp;
        eargs[t].bundles_cb = bun_cb;
        eargs[t].bundles_cs = bun_cs;
        eptrs[t] = &eargs[t];
    }
    double t0 = now_ms();
    run_parallel(nt1, nt1, worker_v_enc, eptrs);
    result_add(&res, "Step 1: Vehicle BlindSig+ABE Enc", now_ms() - t0);
    free(eargs); free(eptrs);

    /* ── Step 2: 6G uplink (OFDMA single slot) */
    result_add(&res, "Step 2: V->RSU 6G uplink",
               network_delay(&net, net.d_v2r, B_RSA + 3 * B_ABE));

    /* ── Step 3: RSU parallel ABE Decrypt  [K RSUs × C cores]
       We distribute N vehicles across K×C = 64 worker threads.        */
    int n_workers = K_RSUS * C_CORES;
    if (n_workers > N) n_workers = N;
    VDecArg *dargs = malloc(n_workers * sizeof(VDecArg));
    void   **dptrs = malloc(n_workers * sizeof(void *));
    int chunk3 = (N + n_workers - 1) / n_workers;
    for (int t = 0; t < n_workers; t++) {
        dargs[t].rsa_auth   = rsa_auth;
        dargs[t].bundles_cp = bun_cp;
        dargs[t].bundles_cb = bun_cb;
        dargs[t].bundles_cs = bun_cs;
        dargs[t].start      = t * chunk3;
        dargs[t].end        = (t+1)*chunk3 < N ? (t+1)*chunk3 : N;
        dargs[t].dec_cp     = dec_cp;
        dargs[t].dec_cb     = dec_cb;
        dargs[t].dec_cs     = dec_cs;
        dptrs[t] = &dargs[t];
    }
    t0 = now_ms();
    run_parallel(n_workers, n_workers, worker_v_dec, dptrs);
    result_add(&res, "Step 3: RSU parallel ABE Decrypt", now_ms() - t0);
    free(dargs); free(dptrs);

    /* ── Step 4: RSU aggregation */
    t0 = now_ms();
    int64_t CP_a = 0, CB_a = 0, CS_a = 0;
    for (int i = 0; i < N; i++) {
        CP_a += dec_cp[i]; CB_a += dec_cb[i]; CS_a += dec_cs[i];
    }
    result_add(&res, "Step 4: RSU aggregation", now_ms() - t0);

    /* ── Step 5: RSU cached model inference (no CA)
       w = [0.4, 0.3, 0.3], strat_i = dot(w, [CP%1000, CB, CS])       */
    t0 = now_ms();
    for (int i = 0; i < N; i++) {
        double cp = (double)(dec_cp[i] % 1000);
        double cb = (double) dec_cb[i];
        double cs = (double) dec_cs[i];
        strats[i] = 0.4*cp + 0.3*cb + 0.3*cs;
    }
    result_add(&res, "Step 5: RSU cached model inference", now_ms() - t0);

    /* ── Step 6: FE Enc + 6G downlink
       Eq 7: Y_enc = FE.Enc(pk_vehicles, f_predict(w, D''))            */
    t0 = now_ms();
    for (int i = 0; i < N; i++) {
        int64_t val = (int64_t)(fabs(strats[i]) * 100) % ((int64_t)1 << 60);
        fe_enc(rsa_fe, &enc_strats[i], val);
    }
    double step6 = (now_ms() - t0)
                 + network_delay(&net, net.d_v2r, B_STRAT * N);
    result_add(&res, "Step 6: RSU FE Enc+6G downlink", step6);

    /* ── Step 7: Vehicle FE Decrypt  [N parallel, ≤16 threads]
       Eq 7: Y_dec,i = FE.Dec(sk_Vi, Y_enc) = {Strategy_CP,CB,CS}     */
    int nt7 = N < 16 ? N : 16;
    FDecArg *fargs = malloc(nt7 * sizeof(FDecArg));
    void   **fptrs = malloc(nt7 * sizeof(void *));
    int chunk7 = (N + nt7 - 1) / nt7;
    for (int t = 0; t < nt7; t++) {
        fargs[t].rsa_fe     = rsa_fe;
        fargs[t].enc_strats = enc_strats;
        fargs[t].start      = t * chunk7;
        fargs[t].end        = (t+1)*chunk7 < N ? (t+1)*chunk7 : N;
        fptrs[t] = &fargs[t];
    }
    t0 = now_ms();
    run_parallel(nt7, nt7, worker_fe_dec, fptrs);
    result_add(&res, "Step 7: Vehicle FE Decrypt", now_ms() - t0);
    free(fargs); free(fptrs);

    for (int i = 0; i < N; i++) mpz_clear(bsigs[i]);
    free(vehs); free(bsigs);
    free(bun_cp); free(bun_cb); free(bun_cs);
    free(dec_cp); free(dec_cb); free(dec_cs);
    free(strats); free(enc_strats);

    return res;
}

/* ════════════════════════════════════════════════════════════════════
   SLOW LOOP  (background CA refresh)
   ════════════════════════════════════════════════════════════════════ */
static LoopResult slow_loop(int N,
                             const RSACtx *rsa_auth,
                             const RSACtx *rsa_fe)
{
    LoopResult res;
    memset(&res, 0, sizeof(res));
    res.N = N; res.K = K_RSUS; res.C = C_CORES;
    strcpy(res.loop, "slow");

    Network net = network_default();

    Vehicle *vehs = malloc(N * sizeof(Vehicle));
    gen_vehicles(vehs, N, 42);

    int64_t CP_a = 0, CB_a = 0, CS_a = 0;
    for (int i = 0; i < N; i++) {
        CP_a += vehs[i].D_CP;
        CB_a += vehs[i].D_CB;
        CS_a += vehs[i].D_CS;
    }

    /* ── Step A: RSU PRE transform  Eq 9: CT_CA = ReEncrypt(rk, CT_RSU) */
    double t0 = now_ms();
    AbeBundle bP, bB, bS;
    abe_enc(rsa_auth, &bP, CP_a % ((int64_t)1 << 60));
    abe_enc(rsa_auth, &bB, CB_a % ((int64_t)1 << 60));
    abe_enc(rsa_auth, &bS, CS_a % ((int64_t)1 << 60));

    /* pre_transform operates on the RSA-encapsulated key (k_enc field) */
    mpz_t ke_P, ke_B, ke_S, pre_P, pre_B, pre_S;
    mpz_inits(ke_P, ke_B, ke_S, pre_P, pre_B, pre_S, NULL);

    mpz_import(ke_P, RSA_BYTES, 1, 1, 0, 0, bP.k_enc);
    mpz_import(ke_B, RSA_BYTES, 1, 1, 0, 0, bB.k_enc);
    mpz_import(ke_S, RSA_BYTES, 1, 1, 0, 0, bS.k_enc);

    pre_transform(rsa_auth, pre_P, ke_P);
    pre_transform(rsa_auth, pre_B, ke_B);
    pre_transform(rsa_auth, pre_S, ke_S);
    result_add(&res, "Step A: RSU PRE transform", now_ms() - t0);

    /* ── Step B: RSU → BS 6G link */
    result_add(&res, "Step B: RSU->BS 6G link",
               network_delay(&net, net.d_r2b, 3 * B_ABE));

    /* ── Step C: BS → CA DTLS + backhaul */
    result_add(&res, "Step C: BS->CA DTLS+backhaul",
               0.5 + network_delay(&net, net.d_b2c, 3 * B_ABE));

    /* ── Step D: CA PRE Decrypt */
    t0 = now_ms();
    mpz_t dec_P, dec_B, dec_S;
    mpz_inits(dec_P, dec_B, dec_S, NULL);
    rsa_dec(rsa_auth, dec_P, pre_P);
    rsa_dec(rsa_auth, dec_B, pre_B);
    rsa_dec(rsa_auth, dec_S, pre_S);
    result_add(&res, "Step D: CA PRE Decrypt", now_ms() - t0);

    /* ── Step E: CA ML update
       Eq 10: w_new = f_ML(w_old, CP_agg, CB_agg, CS_agg)
       w_new = w_old - 0.01 * (dot(w_old,X) - 0.5) * X              */
    t0 = now_ms();
    double w[3]  = {0.4, 0.3, 0.3};
    double X[3]  = { (double)(CP_a % 1000),
                     (double)(CB_a % 2),
                     (double)(CS_a % 4) };
    double norm  = sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]) + 1e-9;
    for (int j = 0; j < 3; j++) X[j] /= norm;
    double dot   = w[0]*X[0] + w[1]*X[1] + w[2]*X[2];
    double w_new[3];
    for (int j = 0; j < 3; j++)
        w_new[j] = w[j] - 0.01 * (dot - 0.5) * X[j];
    result_add(&res, "Step E: CA ML update", now_ms() - t0);

    /* ── Step F: FE Enc + model push CA → BS → RSU
       Eq 7 variant: encode each updated weight                        */
    t0 = now_ms();
    AbeBundle wb[3];
    for (int j = 0; j < 3; j++) {
        int64_t wval = (int64_t)(fabs(w_new[j]) * 1e6) % ((int64_t)1 << 60);
        fe_enc(rsa_fe, &wb[j], wval);
    }
    double step_f = (now_ms() - t0)
                  + network_delay(&net, net.d_b2c, 3 * B_STRAT)
                  + network_delay(&net, net.d_r2b, 3 * B_STRAT);
    result_add(&res, "Step F: CA FE Enc+model push", step_f);

    mpz_clears(ke_P, ke_B, ke_S, pre_P, pre_B, pre_S,
               dec_P, dec_B, dec_S, NULL);
    free(vehs);

    return res;
}


static void write_json(const char *path,
                       LoopResult *results, int count)
{
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "Cannot open %s\n", path); return; }

    fprintf(f, "[\n");
    for (int i = 0; i < count; i++) {
        LoopResult *r = &results[i];
        fprintf(f, "  {\n");
        fprintf(f, "    \"N\": %d,\n", r->N);
        fprintf(f, "    \"K\": %d,\n", r->K);
        fprintf(f, "    \"C\": %d,\n", r->C);
        fprintf(f, "    \"loop\": \"%s\",\n", r->loop);
        fprintf(f, "    \"steps\": {\n");
        for (int s = 0; s < r->n_steps; s++) {
            fprintf(f, "      \"%s\": %.4f%s\n",
                    r->steps[s].name, r->steps[s].ms,
                    s < r->n_steps - 1 ? "," : "");
        }
        fprintf(f, "    },\n");
        fprintf(f, "    \"total_ms\": %.4f\n", r->total_ms);
        fprintf(f, "  }%s\n", i < count - 1 ? "," : "");
    }
    fprintf(f, "]\n");
    fclose(f);
}

/* ════════════════════════════════════════════════════════════════════
   main
   ════════════════════════════════════════════════════════════════════ */
int main(void)
{
    printf("========================================================================\n");
    printf("  SOFT PRIVACY v3 — REAL NETWORK SIMULATION\n");
    printf("  CP-ABE: RSA-1024 KEM+AES-256-GCM | PRE: RSA-1024 | FE: RSA-1024+AES\n");
    printf("  6G: 100 Gbps 3GPP TR38.824 | K=%d RSUs  C=%d cores/RSU\n",
           K_RSUS, C_CORES);
    printf("========================================================================\n");

    printf("  [keygen] Generating infrastructure RSA keys...\n");
    fflush(stdout);

    RSACtx rsa_auth, rsa_fe;
    rsa_init(&rsa_auth);
    rsa_init(&rsa_fe);
    printf("  [keygen] Done.\n");
    fflush(stdout);

    LoopResult fast_all[N_COUNTS];
    LoopResult slow_all[N_COUNTS];

    for (int ni = 0; ni < N_COUNTS; ni++) {
        int N = VEHICLE_COUNTS[ni];
        printf("\n-- N=%d --\n", N);
        fflush(stdout);

        /* fast loop */
        fast_all[ni] = fast_loop(N, &rsa_auth, &rsa_fe);
        for (int s = 0; s < fast_all[ni].n_steps; s++)
            printf("  [F] %-43s %8.2f ms\n",
                   fast_all[ni].steps[s].name,
                   fast_all[ni].steps[s].ms);
        printf("  [F] %-43s %8.2f ms\n", "Total E2E",
               fast_all[ni].total_ms);

        /* slow loop */
        slow_all[ni] = slow_loop(N, &rsa_auth, &rsa_fe);
        for (int s = 0; s < slow_all[ni].n_steps; s++)
            printf("  [S] %-43s %8.2f ms\n",
                   slow_all[ni].steps[s].name,
                   slow_all[ni].steps[s].ms);
        printf("  [S] %-43s %8.2f ms\n", "Total slow",
               slow_all[ni].total_ms);

        fflush(stdout);
    }

    /* write JSON results */
    write_json("soft_v3_fast.json", fast_all, N_COUNTS);
    write_json("soft_v3_slow.json", slow_all, N_COUNTS);

    printf("\n[OK] Saved: soft_v3_fast.json  soft_v3_slow.json\n");

    rsa_free(&rsa_auth);
    rsa_free(&rsa_fe);
    return 0;
}
