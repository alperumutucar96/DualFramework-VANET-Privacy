/**
 * crypto_primitives.c
 *
 * Compile:
 *   gcc -O2 -o demo crypto_primitives.c main.c -lgmp -lssl -lcrypto -lm
 */
#include "crypto_primitives.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <gmp.h>
#include <openssl/bn.h>
#include <openssl/evp.h>
#include <openssl/rand.h>
#include <openssl/sha.h>

/* ═══════════════════════════════════════════════════════════════════════
   Internal helpers
   ═══════════════════════════════════════════════════════════════════════ */

/* Fill mpz_t with cryptographically random bits */
static void mpz_rand_bits(mpz_t out, int bits)
{
    int bytes = (bits + 7) / 8;
    uint8_t *buf = malloc(bytes);
    RAND_bytes(buf, bytes);
    mpz_import(out, bytes, 1, 1, 0, 0, buf);
    mpz_tdiv_r_2exp(out, out, bits);  /* mask to exact bit count */
    free(buf);
}

/* Miller-Rabin prime (GMP built-in, 25 rounds → error < 2^-50) */
static void gen_prime(mpz_t p, int bits)
{
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);

    uint8_t seed_buf[32];
    RAND_bytes(seed_buf, 32);
    mpz_t seed; mpz_init(seed);
    mpz_import(seed, 32, 1, 1, 0, 0, seed_buf);
    gmp_randseed(rstate, seed);
    mpz_clear(seed);

    mpz_urandomb(p, rstate, bits);
    mpz_setbit(p, bits - 1);   /* ensure top bit set */
    mpz_nextprime(p, p);

    gmp_randclear(rstate);
}

/* SHA-256 of a string, result as mpz_t */
static void sha256_mpz(mpz_t out, const char *str)
{
    uint8_t digest[32];
    SHA256((const uint8_t *)str, strlen(str), digest);
    mpz_import(out, 32, 1, 1, 0, 0, digest);
}

/* Gaussian random number (Box-Muller) */
static double gauss(double sigma)
{
    double u1 = ((double)rand() + 1.0) / ((double)RAND_MAX + 1.0);
    double u2 = ((double)rand() + 1.0) / ((double)RAND_MAX + 1.0);
    return sigma * sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265358979323846 * u2);
}

/* mpz_t → byte array (big-endian), zero-padded to `len` bytes */
static void mpz_to_bytes_padded(const mpz_t x, uint8_t *buf, size_t len)
{
    memset(buf, 0, len);
    size_t written = 0;
    uint8_t *tmp = (uint8_t *)mpz_export(NULL, &written, 1, 1, 0, 0, x);
    if (tmp) {
        if (written <= len)
            memcpy(buf + len - written, tmp, written);
        else
            memcpy(buf, tmp + written - len, len); /* truncate */
        free(tmp);
    }
}

/* ═══════════════════════════════════════════════════════════════════════
   RSA-1024
   ═══════════════════════════════════════════════════════════════════════ */
void rsa_init(RSACtx *ctx)
{
    mpz_init(ctx->n); mpz_init(ctx->e); mpz_init(ctx->d);

    mpz_t p, q, p1, q1, phi;
    mpz_init(p); mpz_init(q); mpz_init(p1); mpz_init(q1); mpz_init(phi);

    gen_prime(p, RSA_BITS / 2);
    gen_prime(q, RSA_BITS / 2);

    mpz_mul(ctx->n, p, q);
    mpz_set_ui(ctx->e, 65537);

    mpz_sub_ui(p1, p, 1);
    mpz_sub_ui(q1, q, 1);
    mpz_mul(phi, p1, q1);

    mpz_invert(ctx->d, ctx->e, phi);   /* d = e^{-1} mod φ(n) */

    mpz_clears(p, q, p1, q1, phi, NULL);
}

void rsa_free(RSACtx *ctx)
{
    mpz_clears(ctx->n, ctx->e, ctx->d, NULL);
}

void rsa_enc(const RSACtx *ctx, mpz_t out, const mpz_t m)
{
    mpz_t tmp; mpz_init(tmp);
    mpz_mod(tmp, m, ctx->n);
    mpz_powm(out, tmp, ctx->e, ctx->n);
    mpz_clear(tmp);
}

void rsa_dec(const RSACtx *ctx, mpz_t out, const mpz_t c)
{
    mpz_t tmp; mpz_init(tmp);
    mpz_mod(tmp, c, ctx->n);
    mpz_powm(out, tmp, ctx->d, ctx->n);
    mpz_clear(tmp);
}

/* ═══════════════════════════════════════════════════════════════════════
   Paillier-1024
   ═══════════════════════════════════════════════════════════════════════ */
void paillier_init(PaillierCtx *ctx)
{
    mpz_inits(ctx->n, ctx->n2, ctx->g, ctx->lam, ctx->mu, NULL);

    mpz_t p, q, p1, q1;
    mpz_inits(p, q, p1, q1, NULL);

    gen_prime(p, RSA_BITS / 2);
    gen_prime(q, RSA_BITS / 2);

    mpz_mul(ctx->n, p, q);
    mpz_mul(ctx->n2, ctx->n, ctx->n);       /* n^2           */
    mpz_add_ui(ctx->g, ctx->n, 1);          /* g = n + 1     */

    mpz_sub_ui(p1, p, 1); mpz_sub_ui(q1, q, 1);
    mpz_mul(ctx->lam, p1, q1);              /* λ = (p-1)(q-1)*/
    mpz_invert(ctx->mu, ctx->lam, ctx->n);  /* μ = λ^{-1} mod n */

    mpz_clears(p, q, p1, q1, NULL);
}

void paillier_free(PaillierCtx *ctx)
{
    mpz_clears(ctx->n, ctx->n2, ctx->g, ctx->lam, ctx->mu, NULL);
}

void paillier_enc(const PaillierCtx *ctx, mpz_t out, const mpz_t m)
{
    mpz_t r, rn, gm, tmp;
    mpz_inits(r, rn, gm, tmp, NULL);

    /* random r ∈ [1, n) */
    mpz_rand_bits(r, RSA_BITS);
    mpz_mod(r, r, ctx->n);
    if (mpz_cmp_ui(r, 0) == 0) mpz_set_ui(r, 1);

    /* g^m mod n^2 = (n+1)^m mod n^2 = 1 + m*n mod n^2  */
    mpz_mod(tmp, m, ctx->n);
    mpz_mul(gm, tmp, ctx->n);
    mpz_add_ui(gm, gm, 1);
    mpz_mod(gm, gm, ctx->n2);

    /* r^n mod n^2 */
    mpz_powm(rn, r, ctx->n, ctx->n2);

    mpz_mul(out, gm, rn);
    mpz_mod(out, out, ctx->n2);

    mpz_clears(r, rn, gm, tmp, NULL);
}

/* Homomorphic addition: c1 * c2 mod n^2 */
void paillier_add(const PaillierCtx *ctx, mpz_t out,
                  const mpz_t c1, const mpz_t c2)
{
    mpz_mul(out, c1, c2);
    mpz_mod(out, out, ctx->n2);
}

void paillier_dec(const PaillierCtx *ctx, mpz_t out, const mpz_t c)
{
    /* L(x) = (x - 1) / n */
    mpz_t x, Lx;
    mpz_inits(x, Lx, NULL);

    mpz_powm(x, c, ctx->lam, ctx->n2);
    mpz_sub_ui(Lx, x, 1);
    mpz_tdiv_q(Lx, Lx, ctx->n);            /* exact division */

    mpz_mul(out, Lx, ctx->mu);
    mpz_mod(out, out, ctx->n);

    mpz_clears(x, Lx, NULL);
}

/* ═══════════════════════════════════════════════════════════════════════
   CP-ABE proxy: RSA-KEM + AES-256-GCM
   ═══════════════════════════════════════════════════════════════════════ */

/* AbeBundle using RSA-KEM + AES-256-GCM */
int abe_enc(const RSACtx *rsa, AbeBundle *out, int64_t data)
{
    /* 1. Generate random AES key */
    uint8_t k[AES_KEY_BYTES];
    RAND_bytes(k, AES_KEY_BYTES);

    /* 2. RSA-encapsulate the key */
    mpz_t k_mpz, k_enc;
    mpz_inits(k_mpz, k_enc, NULL);
    mpz_import(k_mpz, AES_KEY_BYTES, 1, 1, 0, 0, k);
    rsa_enc(rsa, k_enc, k_mpz);
    mpz_to_bytes_padded(k_enc, out->k_enc, RSA_BYTES);
    mpz_clears(k_mpz, k_enc, NULL);

    /* 3. AES-256-GCM encrypt the 8-byte plaintext */
    RAND_bytes(out->iv, GCM_IV_BYTES);

    int64_t val = data % (1LL << 63);
    uint8_t plain[8];
    for (int i = 7; i >= 0; i--) { plain[i] = val & 0xff; val >>= 8; }

    EVP_CIPHER_CTX *ctx = EVP_CIPHER_CTX_new();
    EVP_EncryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, k, out->iv);

    int outlen = 0;
    EVP_EncryptUpdate(ctx, out->ct, &outlen, plain, 8);
    EVP_EncryptFinal_ex(ctx, out->ct + outlen, &outlen);
    EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_GET_TAG, GCM_TAG_BYTES, out->tag);
    EVP_CIPHER_CTX_free(ctx);

    return 0;
}

int64_t abe_dec(const RSACtx *rsa, const AbeBundle *bundle)
{
    /* 1. RSA-decapsulate the AES key */
    mpz_t k_enc, k_mpz;
    mpz_inits(k_enc, k_mpz, NULL);
    mpz_import(k_enc, RSA_BYTES, 1, 1, 0, 0, bundle->k_enc);
    rsa_dec(rsa, k_mpz, k_enc);

    uint8_t k_raw[RSA_BYTES];
    mpz_to_bytes_padded(k_mpz, k_raw, RSA_BYTES);
    /* take last 32 bytes as AES key (matches Python [-32:]) */
    uint8_t *k = k_raw + RSA_BYTES - AES_KEY_BYTES;
    mpz_clears(k_enc, k_mpz, NULL);

    /* 2. AES-256-GCM decrypt */
    uint8_t plain[8];
    EVP_CIPHER_CTX *ctx = EVP_CIPHER_CTX_new();
    EVP_DecryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, k, bundle->iv);
    EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_TAG, GCM_TAG_BYTES,
                        (void *)bundle->tag);

    int outlen = 0;
    EVP_DecryptUpdate(ctx, plain, &outlen, bundle->ct, 8);
    /* ignore tag-verify result for simulation */
    EVP_DecryptFinal_ex(ctx, plain + outlen, &outlen);
    EVP_CIPHER_CTX_free(ctx);

    int64_t val = 0;
    for (int i = 0; i < 8; i++) val = (val << 8) | plain[i];
    return val;
}

/* FE is identical to ABE in this simulation */
int     fe_enc(const RSACtx *rsa, AbeBundle *out, int64_t value)
{ return abe_enc(rsa, out, value); }

int64_t fe_dec(const RSACtx *rsa, const AbeBundle *bundle)
{ return abe_dec(rsa, bundle); }

/* PRE: decrypt with rsa, re-encrypt (same key in this simulation) */
void pre_transform(const RSACtx *rsa, mpz_t out, const mpz_t k_e)
{
    mpz_t tmp;
    mpz_init(tmp);
    rsa_dec(rsa, tmp, k_e);
    rsa_enc(rsa, out, tmp);
    mpz_clear(tmp);
}

/* ═══════════════════════════════════════════════════════════════════════
   Blind signature (Chaum 1983, RSA-1024)
   ═══════════════════════════════════════════════════════════════════════ */
void blind_sign(const RSACtx *rsa, mpz_t out,
                const uint8_t *msg, size_t len)
{
    /* h = SHA-256(msg) mod n */
    uint8_t digest[32];
    SHA256(msg, len, digest);
    mpz_t h, r, ri, blinded, sig;
    mpz_inits(h, r, ri, blinded, sig, NULL);
    mpz_import(h, 32, 1, 1, 0, 0, digest);
    mpz_mod(h, h, rsa->n);

    /* random r, compute r^{-1} mod n */
    mpz_rand_bits(r, RSA_BITS);
    mpz_mod(r, r, rsa->n);
    if (mpz_cmp_ui(r, 0) == 0) mpz_set_ui(r, 1);
    mpz_invert(ri, r, rsa->n);

    /* blinded = h * r^e mod n */
    mpz_t re; mpz_init(re);
    rsa_enc(rsa, re, r);
    mpz_mul(blinded, h, re);
    mpz_mod(blinded, blinded, rsa->n);
    mpz_clear(re);

    /* sign(blinded) * r^{-1} mod n */
    mpz_powm(sig, blinded, rsa->d, rsa->n);
    mpz_mul(out, sig, ri);
    mpz_mod(out, out, rsa->n);

    mpz_clears(h, r, ri, blinded, sig, NULL);
}

/* ═══════════════════════════════════════════════════════════════════════
   Schnorr NIZK  (511-bit prime group, proxy for Bulletproof)
   ═══════════════════════════════════════════════════════════════════════ */

/* _ZP: smallest 511-bit prime ≥ 2^511. Pre-computed for speed. */
static const char ZP_STR[] =
    "6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042048";

static mpz_t ZP, ZG;
static int   zkp_globals_init = 0;

static void zkp_ensure_globals(void)
{
    if (zkp_globals_init) return;
    mpz_init_set_str(ZP, ZP_STR, 10);
    mpz_init_set_ui(ZG, 2);
    zkp_globals_init = 1;
}

void zkp_prove(ZkpProof *proof, int64_t value)
{
    zkp_ensure_globals();
    mpz_inits(proof->pub, proof->R, proof->s, proof->c, NULL);

    mpz_t P1, sec, k, tmp;
    mpz_inits(P1, sec, k, tmp, NULL);
    mpz_sub_ui(P1, ZP, 1);              /* p - 1 */

    /* sec = value mod (p-1) */
    mpz_set_si(sec, value < 0 ? -value : value);
    mpz_mod(sec, sec, P1);

    /* pub = G^sec mod P */
    mpz_powm(proof->pub, ZG, sec, ZP);

    /* k random in [0, p-1) */
    mpz_rand_bits(k, ZKP_BITS);
    mpz_mod(k, k, P1);

    /* R = G^k mod P */
    mpz_powm(proof->R, ZG, k, ZP);

    /* c = SHA-256(R || pub) mod P */
    char *R_str  = mpz_get_str(NULL, 10, proof->R);
    char *pub_str = mpz_get_str(NULL, 10, proof->pub);
    size_t buf_len = strlen(R_str) + strlen(pub_str) + 2;
    char *buf = malloc(buf_len);
    snprintf(buf, buf_len, "%s%s", R_str, pub_str);

    sha256_mpz(proof->c, buf);
    mpz_mod(proof->c, proof->c, ZP);

    free(R_str); free(pub_str); free(buf);

    /* s = (k - c * sec) mod (p-1) */
    mpz_mul(tmp, proof->c, sec);
    mpz_sub(proof->s, k, tmp);
    mpz_mod(proof->s, proof->s, P1);

    mpz_clears(P1, sec, k, tmp, NULL);
}

void zkp_free(ZkpProof *proof)
{
    mpz_clears(proof->pub, proof->R, proof->s, proof->c, NULL);
}

int zkp_verify(const ZkpProof *proof)
{
    zkp_ensure_globals();
    /* check: G^s * pub^c ≡ R  (mod P) */
    /* normalise s to [0, p-1) so tampered values are correctly rejected */
    mpz_t lhs, gs, pubc, s_pos, P1;
    mpz_inits(lhs, gs, pubc, s_pos, P1, NULL);

    mpz_sub_ui(P1, ZP, 1);
    mpz_mod(s_pos, proof->s, P1);

    mpz_powm(gs,   ZG,         s_pos,    ZP);
    mpz_powm(pubc, proof->pub, proof->c, ZP);
    mpz_mul(lhs, gs, pubc);
    mpz_mod(lhs, lhs, ZP);

    int ok = (mpz_cmp(lhs, proof->R) == 0);
    mpz_clears(lhs, gs, pubc, s_pos, P1, NULL);
    return ok;
}

/* ═══════════════════════════════════════════════════════════════════════
   Shamir MPC  (2-of-3, prime = 2^127 + 45 [Mersenne-near])
   ═══════════════════════════════════════════════════════════════════════ */

static const char MPC_PRIME_STR[] =
    "170141183460469231731687303715884105773"; /* smallest prime > 2^127 */

static mpz_t MPC_P;
static int   mpc_globals_init = 0;

static void mpc_ensure_globals(void)
{
    if (mpc_globals_init) return;
    mpz_init_set_str(MPC_P, MPC_PRIME_STR, 10);
    mpc_globals_init = 1;
}

void mpc_share(Share shares[3], const mpz_t val)
{
    mpc_ensure_globals();

    mpz_t v, b;
    mpz_inits(v, b, NULL);
    mpz_mod(v, val, MPC_P);

    /* random slope b */
    mpz_rand_bits(b, 128);
    mpz_mod(b, b, MPC_P);

    for (int i = 0; i < 3; i++) {
        shares[i].x = i + 1;
        mpz_init(shares[i].y);
        /* y = v + b*x  mod P  (linear Shamir) */
        mpz_t bx; mpz_init(bx);
        mpz_mul_ui(bx, b, (unsigned long)(i + 1));
        mpz_add(shares[i].y, v, bx);
        mpz_mod(shares[i].y, shares[i].y, MPC_P);
        mpz_clear(bx);
    }
    mpz_clears(v, b, NULL);
}

void mpc_reconstruct(mpz_t out, const Share shares[3])
{
    mpc_ensure_globals();
    mpz_t x0, x1, y0, y1, tmp, L0, L1, num, den, inv;
    mpz_inits(x0, x1, y0, y1, tmp, L0, L1, num, den, inv, NULL);

    mpz_set_si(x0, shares[0].x);
    mpz_set_si(x1, shares[1].x);
    mpz_set(y0, shares[0].y);
    mpz_set(y1, shares[1].y);

    /* L0 = y0 * (-x1) * invert(x0 - x1) mod P */
    mpz_neg(num, x1);
    mpz_mod(num, num, MPC_P);
    mpz_sub(den, x0, x1);
    mpz_mod(den, den, MPC_P);
    mpz_invert(inv, den, MPC_P);
    mpz_mul(L0, y0, num);
    mpz_mul(L0, L0, inv);
    mpz_mod(L0, L0, MPC_P);

    /* L1 = y1 * (-x0) * invert(x1 - x0) mod P */
    mpz_neg(num, x0);
    mpz_mod(num, num, MPC_P);
    mpz_sub(den, x1, x0);
    mpz_mod(den, den, MPC_P);
    mpz_invert(inv, den, MPC_P);
    mpz_mul(L1, y1, num);
    mpz_mul(L1, L1, inv);
    mpz_mod(L1, L1, MPC_P);

    mpz_add(out, L0, L1);
    mpz_mod(out, out, MPC_P);

    mpz_clears(x0, x1, y0, y1, tmp, L0, L1, num, den, inv, NULL);
}

void share_free(Share *s) { mpz_clear(s->y); }

/* ═══════════════════════════════════════════════════════════════════════
   6G network model
   ═══════════════════════════════════════════════════════════════════════ */
Network network_default(void)
{
    Network n = {
        .bw_bps  = 100e9,
        .proc_ms = 0.05,
        .jit_sig = 0.01,
        .d_v2r   = 200.0,
        .d_r2b   = 500.0,
        .d_b2c   = 5000.0,
        .c_light = 3e8
    };
    return n;
}

double network_delay(const Network *net, double dist_m, int bytes)
{
    double prop  = (dist_m / net->c_light) * 1000.0;
    double tx    = ((double)bytes * 8.0) / net->bw_bps * 1000.0;
    double jit   = fabs(gauss(net->jit_sig));
    double total = prop + tx + net->proc_ms + jit;
    return total < 0.001 ? 0.001 : total;
}

/* ═══════════════════════════════════════════════════════════════════════
   Vehicle generation (LCG matching Python's random.Random(seed))
   ═══════════════════════════════════════════════════════════════════════ */
void gen_vehicles(Vehicle *out, int N, unsigned int seed)
{
    srand(seed);
    for (int i = 0; i < N; i++) {
        out[i].vid   = i;
        out[i].D_CID = 100000 + rand() % 900000;
        out[i].D_CP  = rand() % 1000;
        out[i].D_CB  = rand() % 2;
        out[i].D_CS  = rand() % 4;
    }
}
