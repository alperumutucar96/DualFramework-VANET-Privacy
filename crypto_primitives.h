/**
 * crypto_primitives.h
 *
 * Dependencies: OpenSSL (RSA, AES-256-GCM, SHA-256), GMP (big integers)
 * Key sizes: RSA/Paillier 1024-bit (simulation); deploy at 2048-bit.
 */
#ifndef CRYPTO_PRIMITIVES_H
#define CRYPTO_PRIMITIVES_H

#include <stdint.h>
#include <stddef.h>
#include <gmp.h>
#include <openssl/rsa.h>
#include <openssl/bn.h>
#include <openssl/evp.h>

/* ── constants ─────────────────────────────────────────────────────── */
#define RSA_BITS       1024
#define RSA_BYTES      (RSA_BITS / 8)          /* 128 */
#define PAILLIER_BYTES (RSA_BITS * 2 / 8)      /* 256 */
#define AES_KEY_BYTES  32
#define GCM_IV_BYTES   16
#define GCM_TAG_BYTES  16
#define ZKP_BITS       511

/* payload sizes (bytes) */
#define B_RSA    RSA_BYTES
#define B_ABE    (RSA_BYTES + GCM_IV_BYTES + 8 + GCM_TAG_BYTES) /* 168 */
#define B_PAIL   PAILLIER_BYTES
#define B_ZKP    (3 * RSA_BYTES)               /* 384 */
#define B_STRAT  128

/* ── RSA context ───────────────────────────────────────────────────── */
typedef struct {
    mpz_t n, e, d;
} RSACtx;

/* ── Paillier context ──────────────────────────────────────────────── */
typedef struct {
    mpz_t n, n2, g, lam, mu;
} PaillierCtx;

/* ── ABE / FE ciphertext bundle ────────────────────────────────────── */
typedef struct {
    uint8_t k_enc[RSA_BYTES];   /* RSA-encapsulated AES key  */
    uint8_t iv[GCM_IV_BYTES];   /* AES-GCM nonce             */
    uint8_t ct[8];              /* AES-GCM ciphertext        */
    uint8_t tag[GCM_TAG_BYTES]; /* AES-GCM auth tag          */
} AbeBundle;

/* ── Schnorr / NIZK proof ──────────────────────────────────────────── */
typedef struct {
    mpz_t pub, R, s, c;
} ZkpProof;

/* ── Shamir share ──────────────────────────────────────────────────── */
typedef struct { int x; mpz_t y; } Share;

/* ── Vehicle ───────────────────────────────────────────────────────── */
typedef struct {
    int     vid;
    int64_t D_CID, D_CP, D_CB, D_CS;
} Vehicle;

/* ── 6G network model ──────────────────────────────────────────────── */
typedef struct {
    double bw_bps;   
    double proc_ms;  
    double jit_sig; 
    double d_v2r;   
    double d_r2b;   
    double d_b2c;  
    double c_light; 
} Network;

/* ── function prototypes ───────────────────────────────────────────── */

/* RSA */
void   rsa_init  (RSACtx *ctx);
void   rsa_free  (RSACtx *ctx);
void   rsa_enc   (const RSACtx *ctx, mpz_t out, const mpz_t m);
void   rsa_dec   (const RSACtx *ctx, mpz_t out, const mpz_t c);

/* Paillier */
void   paillier_init (PaillierCtx *ctx);
void   paillier_free (PaillierCtx *ctx);
void   paillier_enc  (const PaillierCtx *ctx, mpz_t out, const mpz_t m);
void   paillier_add  (const PaillierCtx *ctx, mpz_t out,
                      const mpz_t c1, const mpz_t c2);
void   paillier_dec  (const PaillierCtx *ctx, mpz_t out, const mpz_t c);

/* ABE / FE */
int    abe_enc   (const RSACtx *rsa, AbeBundle *out, int64_t data);
int64_t abe_dec  (const RSACtx *rsa, const AbeBundle *bundle);
int    fe_enc    (const RSACtx *rsa, AbeBundle *out, int64_t value);
int64_t fe_dec   (const RSACtx *rsa, const AbeBundle *bundle);
void   pre_transform(const RSACtx *rsa, mpz_t out, const mpz_t k_e);

/* Blind signature */
void   blind_sign(const RSACtx *rsa, mpz_t out, const uint8_t *msg, size_t len);

/* ZKP */
void   zkp_prove  (ZkpProof *proof, int64_t value);
void   zkp_free   (ZkpProof *proof);
int    zkp_verify (const ZkpProof *proof);

/* Shamir MPC */
void   mpc_share       (Share shares[3], const mpz_t val);
void   mpc_reconstruct (mpz_t out, const Share shares[3]);
void   share_free      (Share *s);

/* Network */
Network network_default(void);
double  network_delay  (const Network *net, double dist_m, int bytes);

/* Vehicles */
void gen_vehicles(Vehicle *out, int N, unsigned int seed);

#endif /* CRYPTO_PRIMITIVES_H */
