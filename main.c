/**
 * main.c — demonstration and test harness for crypto_primitives
 *
 * Exercises every primitive and reports timing + correctness.
 * Compile: gcc -O2 -o demo crypto_primitives.c main.c -lgmp -lssl -lcrypto -lm
 * Run:     ./demo
 */
#define _POSIX_C_SOURCE 199309L
#include "crypto_primitives.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

/* ── timing helper ──────────────────────────────────────────────────── */
static double now_ms(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1e6;
}

/* ── pretty separator ───────────────────────────────────────────────── */
static void section(const char *title)
{
    printf("\n══════════════════════════════════════════════════\n");
    printf("  %s\n", title);
    printf("══════════════════════════════════════════════════\n");
}

/* ═══════════════════════════════════════════════════════════════════════
   main
   ═══════════════════════════════════════════════════════════════════════ */
int main(void)
{
    printf("crypto_primitives — 6G Vehicular Network Demo\n");
    printf("(RSA/Paillier 1024-bit simulation keys)\n");

    double t0, t1;

    /* ── 1. RSA-1024 ─────────────────────────────────────────────── */
    section("1. RSA-1024  (proxy for ABE / FE / PRE bilinear pairings)");

    t0 = now_ms();
    RSACtx rsa;
    rsa_init(&rsa);
    t1 = now_ms();
    printf("  Key generation : %.2f ms\n", t1 - t0);

    mpz_t m, c, dec;
    mpz_inits(m, c, dec, NULL);
    mpz_set_ui(m, 123456789UL);

    t0 = now_ms();
    rsa_enc(&rsa, c, m);
    t1 = now_ms();
    printf("  Encrypt (1024b): %.3f ms\n", t1 - t0);

    t0 = now_ms();
    rsa_dec(&rsa, dec, c);
    t1 = now_ms();
    printf("  Decrypt (1024b): %.3f ms\n", t1 - t0);

    int rsa_ok = (mpz_cmp(m, dec) == 0);
    printf("  Correctness    : %s\n", rsa_ok ? "PASS ✓" : "FAIL ✗");
    mpz_clears(m, c, dec, NULL);

    /* ── 2. Paillier-1024 ────────────────────────────────────────── */
    section("2. Paillier-1024  (Eq.11 / Eq.13 HE encryption + eval)");

    t0 = now_ms();
    PaillierCtx paillier;
    paillier_init(&paillier);
    t1 = now_ms();
    printf("  Key generation : %.2f ms\n", t1 - t0);

    mpz_t pa, pb, ca, cb, cagg, da, db, dsum;
    mpz_inits(pa, pb, ca, cb, cagg, da, db, dsum, NULL);
    mpz_set_ui(pa, 300); mpz_set_ui(pb, 450);

    t0 = now_ms();
    paillier_enc(&paillier, ca, pa);
    paillier_enc(&paillier, cb, pb);
    t1 = now_ms();
    printf("  2× Encrypt     : %.3f ms\n", t1 - t0);

    t0 = now_ms();
    paillier_add(&paillier, cagg, ca, cb);
    t1 = now_ms();
    printf("  HE Add (Eq.13) : %.3f ms\n", t1 - t0);

    t0 = now_ms();
    paillier_dec(&paillier, dsum, cagg);
    t1 = now_ms();
    printf("  Decrypt sum    : %.3f ms\n", t1 - t0);

    unsigned long expected = 750, actual = mpz_get_ui(dsum);
    printf("  300+450 = %lu   : %s\n", actual,
           actual == expected ? "PASS ✓" : "FAIL ✗");
    mpz_clears(pa, pb, ca, cb, cagg, da, db, dsum, NULL);

    /* ── 3. CP-ABE / FE (RSA-KEM + AES-256-GCM) ─────────────────── */
    section("3. CP-ABE  (Eq.6 field encryption, RSA-KEM + AES-256-GCM)");

    int64_t plaintext = 987654321LL;
    AbeBundle bundle;

    t0 = now_ms();
    abe_enc(&rsa, &bundle, plaintext);
    t1 = now_ms();
    printf("  ABE Encrypt    : %.3f ms\n", t1 - t0);

    t0 = now_ms();
    int64_t recovered = abe_dec(&rsa, &bundle);
    t1 = now_ms();
    printf("  ABE Decrypt    : %.3f ms\n", t1 - t0);
    printf("  Plaintext=%lld recovered=%lld : %s\n",
           (long long)plaintext, (long long)recovered,
           plaintext == recovered ? "PASS ✓" : "FAIL ✗");

    /* ── 4. PRE ──────────────────────────────────────────────────── */
    section("4. Proxy Re-Encryption  (Eq.9)");

    mpz_t ke, ke2;
    mpz_inits(ke, ke2, NULL);
    mpz_set_ui(ke, 42UL);
    rsa_enc(&rsa, ke, ke);          /* simulate an encapsulated key */

    t0 = now_ms();
    pre_transform(&rsa, ke2, ke);
    t1 = now_ms();
    printf("  PRE transform  : %.3f ms\n", t1 - t0);

    mpz_t orig;
    mpz_init(orig);
    rsa_dec(&rsa, orig, ke2);
    printf("  Original value : %lu\n", mpz_get_ui(orig));
    mpz_clears(ke, ke2, orig, NULL);

    /* ── 5. Blind Signature ──────────────────────────────────────── */
    section("5. Blind Signature  (Eq.5, Chaum 1983)");

    const uint8_t *msg = (const uint8_t *)"vehicle_identity_42";
    mpz_t sig;
    mpz_init(sig);

    t0 = now_ms();
    blind_sign(&rsa, sig, msg, strlen((const char *)msg));
    t1 = now_ms();
    printf("  Blind sign     : %.3f ms\n", t1 - t0);
    gmp_printf("  Signature (hex): %Zx...  (%zu bits)\n",
               sig, mpz_sizeinbase(sig, 2));
    mpz_clear(sig);

    /* ── 6. Schnorr NIZK ─────────────────────────────────────────── */
    section("6. Schnorr NIZK  (Eq.12, proxy for Bulletproof range proof)");

    ZkpProof proof;

    t0 = now_ms();
    zkp_prove(&proof, 42LL);
    t1 = now_ms();
    printf("  Prove          : %.3f ms\n", t1 - t0);

    t0 = now_ms();
    int valid = zkp_verify(&proof);
    t1 = now_ms();
    printf("  Verify         : %.3f ms\n", t1 - t0);
    printf("  Valid proof    : %s\n", valid ? "PASS ✓" : "FAIL ✗");

    /* tamper test */
    mpz_add_ui(proof.s, proof.s, 1);
    printf("  Tampered proof : %s\n", !zkp_verify(&proof) ? "Rejected ✓" : "Accepted ✗");
    zkp_free(&proof);

    /* ── 7. Shamir MPC (2-of-3) ──────────────────────────────────── */
    section("7. Shamir MPC  (Eq.14, 2-of-3 secret sharing)");

    mpz_t secret, recovered2;
    mpz_inits(secret, recovered2, NULL);
    mpz_set_ui(secret, 999999999UL);

    Share shares[3];
    t0 = now_ms();
    mpc_share(shares, secret);
    t1 = now_ms();
    printf("  Share          : %.3f ms\n", t1 - t0);

    t0 = now_ms();
    mpc_reconstruct(recovered2, shares);
    t1 = now_ms();
    printf("  Reconstruct    : %.3f ms\n", t1 - t0);

    printf("  Secret=%lu reconstructed=%lu : %s\n",
           mpz_get_ui(secret), mpz_get_ui(recovered2),
           mpz_cmp(secret, recovered2) == 0 ? "PASS ✓" : "FAIL ✗");

    for (int i = 0; i < 3; i++) share_free(&shares[i]);
    mpz_clears(secret, recovered2, NULL);

    /* ── 8. 6G network latency ────────────────────────────────────── */
    section("8. 6G Network Model  (3GPP TR 38.824)");

    Network net = network_default();
    printf("  V→RSU  (%4.0fm, %3dB ABE  ): %.4f ms\n",
           net.d_v2r, B_ABE,  network_delay(&net, net.d_v2r, B_ABE));
    printf("  RSU→BS (%4.0fm, %3dB Pail ): %.4f ms\n",
           net.d_r2b, B_PAIL, network_delay(&net, net.d_r2b, B_PAIL));
    printf("  BS→CA  (%4.0fm, %3dB ZKP  ): %.4f ms\n",
           net.d_b2c, B_ZKP,  network_delay(&net, net.d_b2c, B_ZKP));
    printf("  V→RSU  (%4.0fm, %3dB Strat): %.4f ms\n",
           net.d_v2r, B_STRAT,network_delay(&net, net.d_v2r, B_STRAT));

    /* ── 9. Vehicle generation + full pipeline (4 vehicles) ─────── */
    section("9. Full Pipeline  — 4 vehicles, all primitives");

    int N = 4;
    Vehicle vehicles[4];
    gen_vehicles(vehicles, N, 42);

    for (int i = 0; i < N; i++) {
        Vehicle *v = &vehicles[i];
        printf("\n  Vehicle %d  CID=%lld  CP=%lld  CB=%lld  CS=%lld\n",
               v->vid,
               (long long)v->D_CID, (long long)v->D_CP,
               (long long)v->D_CB,  (long long)v->D_CS);

        /* Eq.6: ABE encrypt CP field */
        AbeBundle ab;
        abe_enc(&rsa, &ab, v->D_CP);
        int64_t cp_back = abe_dec(&rsa, &ab);
        printf("    ABE(D_CP=%lld) → dec=%lld  %s\n",
               (long long)v->D_CP, (long long)cp_back,
               v->D_CP == cp_back ? "✓" : "✗");

        /* Eq.11: Paillier encrypt CP */
        mpz_t pcp, pcp_enc, pcp_dec;
        mpz_inits(pcp, pcp_enc, pcp_dec, NULL);
        mpz_set_ui(pcp, (unsigned long)v->D_CP);
        paillier_enc(&paillier, pcp_enc, pcp);
        paillier_dec(&paillier, pcp_dec, pcp_enc);
        printf("    Paillier(D_CP=%lu) → dec=%lu  %s\n",
               mpz_get_ui(pcp), mpz_get_ui(pcp_dec),
               mpz_cmp(pcp, pcp_dec) == 0 ? "✓" : "✗");

        /* Eq.12: ZKP for CP */
        ZkpProof zp;
        zkp_prove(&zp, v->D_CP);
        printf("    ZKP proof valid: %s\n", zkp_verify(&zp) ? "✓" : "✗");
        zkp_free(&zp);

        mpz_clears(pcp, pcp_enc, pcp_dec, NULL);
    }

    /* Eq.13: aggregate Paillier ciphertexts */
    printf("\n  Aggregating Paillier ciphertexts (Eq.13)...\n");
    mpz_t agg; mpz_init(agg);
    long long total_cp = 0;

    /* encrypt first vehicle's CP */
    mpz_t tmp_m, tmp_c;
    mpz_inits(tmp_m, tmp_c, NULL);
    mpz_set_ui(tmp_m, (unsigned long)vehicles[0].D_CP);
    paillier_enc(&paillier, agg, tmp_m);
    total_cp += vehicles[0].D_CP;

    for (int i = 1; i < N; i++) {
        mpz_set_ui(tmp_m, (unsigned long)vehicles[i].D_CP);
        paillier_enc(&paillier, tmp_c, tmp_m);
        paillier_add(&paillier, agg, agg, tmp_c);
        total_cp += vehicles[i].D_CP;
    }
    mpz_t agg_dec; mpz_init(agg_dec);
    paillier_dec(&paillier, agg_dec, agg);
    printf("  Sum(CP) expected=%lld got=%lu  %s\n",
           total_cp, mpz_get_ui(agg_dec),
           (long long)mpz_get_ui(agg_dec) == total_cp ? "✓" : "✗");

    /* Eq.14: MPC share the aggregate */
    Share mpc_shares[3];
    mpc_share(mpc_shares, agg_dec);
    mpz_t mpc_out; mpz_init(mpc_out);
    mpc_reconstruct(mpc_out, mpc_shares);
    printf("  MPC reconstruct sum: %lu  %s\n",
           mpz_get_ui(mpc_out),
           mpz_cmp(agg_dec, mpc_out) == 0 ? "✓" : "✗");

    for (int i = 0; i < 3; i++) share_free(&mpc_shares[i]);
    mpz_clears(tmp_m, tmp_c, agg, agg_dec, mpc_out, NULL);

    /* ── cleanup ─────────────────────────────────────────────────── */
    rsa_free(&rsa);
    paillier_free(&paillier);

    printf("\n══════════════════════════════════════════════════\n");
    printf("  Demo complete.\n");
    printf("══════════════════════════════════════════════════\n\n");
    return 0;
}
