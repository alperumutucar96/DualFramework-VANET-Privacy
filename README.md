This repository is a full C11 port of the Python simulation described in "Privacy-Preserving ML in 6G Vehicular Networks". It implements two privacy frameworks using real cryptographic operations and a parameterised 6G network
latency model.
 
Soft Privacy  (Eq. 5-7, 9-10)
    CP-ABE (RSA-KEM + AES-256-GCM), Proxy Re-Encryption, Functional Encryption. The fast loop bypasses the CA entirely; the slow loop runs  as a background CA model-refresh.
 
Hard Privacy  (Eq. 11-16)
    Paillier-1024 homomorphic encryption, Schnorr NIZK proofs, Shamir 2-of-3 MPC. All aggregation happens under encryption; the CA only sees ciphertexts.
 
  Vehicle counts : N = 5, 10, 25, 50, 100, 250, 500
  RSUs           : K = 4
  Cores per RSU  : C = 16
 
--------------------------------------------------------------------------------
1. FILE STRUCTURE
--------------------------------------------------------------------------------
 
  crypto_primitives.h    Public API — all types, constants, function prototypes
  crypto_primitives.c    Core library:
                           RSA-1024, Paillier-1024, ABE/FE (RSA-KEM+AES-GCM),
                           PRE, Blind Signature (Chaum 1983), Schnorr NIZK,
                           Shamir MPC (2-of-3), 6G network model, vehicle gen
  main.c                 Unit-test harness — correctness checks + timing for
                           every primitive
  soft_privacy_v3.c      Soft Privacy simulation: 7-step fast loop +
                           6-step slow loop, pthreads, JSON output
  hard_privacy_v3.c      Hard Privacy simulation: 10-step pipeline,
                           Paillier tree-reduction, ZKP verify, MPC, JSON output



 
--------------------------------------------------------------------------------
2. BUILD
--------------------------------------------------------------------------------
  
     gcc -O2 -o demo.exe crypto_primitives.c main.c -lgmp -lssl -lcrypto -lm -lpthread
 
     gcc -O2 -o soft_v3.exe soft_privacy_v3.c crypto_primitives.c -lgmp -lssl -lcrypto -lm -lpthread
 
     gcc -O2 -o hard_v3.exe hard_privacy_v3.c crypto_primitives.c -lgmp -lssl -lcrypto -lm -lpthread
 
 
--------------------------------------------------------------------------------
3. RUN
--------------------------------------------------------------------------------
 
  Unit tests:
     ./demo.exe
 
  Soft Privacy simulation:
     ./soft_v3.exe
     Outputs: soft_v3_fast.json   soft_v3_slow.json
 
  Hard Privacy simulation:
     ./hard_v3.exe
     Output:  hard_v3_results.json
 
--------------------------------------------------------------------------------
4. CRYPTOGRAPHIC PRIMITIVES
--------------------------------------------------------------------------------
 
  Primitive              C implementation                  Paper equation
  ---------------------  --------------------------------  ----------------
  RSA-1024               GMP mpz_powm + gen_prime          ABE/FE/PRE proxy
  Paillier-1024 HE       GMP, g=n+1 optimisation           Eq. 11, 13
  CP-ABE                 RSA-KEM + AES-256-GCM (OpenSSL)   Eq. 6
  Functional Encryption  Same as ABE proxy                 Eq. 7
  Proxy Re-Encryption    RSA dec then re-enc               Eq. 9
  Blind Signature        Chaum 1983, RSA-1024              Eq. 5
  Schnorr NIZK           511-bit group, SHA-256 FS         Eq. 12, 16
  Shamir MPC (2-of-3)    Lagrange interp, 127-bit prime    Eq. 14
  6G Network Model       Prop + TX delay + Gaussian jitter 3GPP TR 38.824
 
--------------------------------------------------------------------------------
5. SIMULATION ARCHITECTURE
--------------------------------------------------------------------------------
 
  Soft Privacy — Fast Loop (7 steps, Eq. 5-7)
  --------------------------------------------
  Step 1  N vehicles in parallel (<=16 threads):
            BlindSign(sk_auth, D_CID)          
            ABE.Enc(PK_ABE, D_CP / D_CB / D_CS)
  Step 2  6G OFDMA uplink delay (V -> RSU, 200 m)
  Step 3  RSU parallel ABE decrypt, K x C = 64 threads
  Step 4  Plaintext aggregation: CP_agg, CB_agg, CS_agg
  Step 5  Cached model inference: strat = w . [CP, CB, CS]
            w = [0.4, 0.3, 0.3]
  Step 6  FE.Enc each strategy + 6G downlink delay     
  Step 7  N vehicles FE decrypt (<=16 threads)       
 
  Soft Privacy — Slow Loop (6 steps, Eq. 9-10)
  ---------------------------------------------
  Runs in the background. N-independent — processes 3 aggregates only.
 
  Step A  RSU PRE transform on 3 aggregate bundles   
  Step B  RSU -> BS 6G link delay
  Step C  BS -> CA backhaul + DTLS overhead (0.5 ms)
  Step D  CA PRE decrypt
  Step E  CA gradient step:                       
            w_new = w_old - 0.01 * (w.X - 0.5) * X
  Step F  FE.Enc 3 weights + model push CA -> BS -> RSU
 
  Hard Privacy (10 steps, Eq. 11-16)
  -----------------------------------
  Step 1  Per vehicle (<=8 threads):
            BlindSign + ZKP.Prove x3 + Paillier.Enc x4
  Step 2  6G uplink (B_RSA + 3*B_ZKP + 4*B_PAIL per vehicle)
  Step 3  ZKP.Verify all proofs, K x C = 64 threads  
  Step 4  Binary tree reduction of N Paillier CTs,
            3 fields in parallel (log2(N) levels)     
  Step 5  Shamir MPC share + 2-round RTT model         
  Step 6  RSU -> BS -> CA DTLS + link delays
  Step 7  CA: dec aggregates -> gradient -> re-enc     
  Step 8  CA prediction + CA -> BS -> RSU distribution
  Step 9  N vehicles each decrypt 3 weights (<=8 threads)
  Step 10 ACK link delay (V -> RSU, 64 bytes)
 
--------------------------------------------------------------------------------
6. PAYLOAD SIZE CONSTANTS
--------------------------------------------------------------------------------
 
  Constant   Value    Description
  ---------  -------  -----------------------------------------------
  B_RSA      128 B    RSA-1024 ciphertext (blind signature)
  B_ABE      168 B    ABE bundle: KEM(128) + IV(16) + data(8) + tag(16)
  B_PAIL     256 B    Paillier-1024 ciphertext (2048-bit number)
  B_ZKP      384 B    Schnorr proof: pub + R + s  (3 x 128 B)
  B_STRAT    128 B    Driving strategy packet
 
--------------------------------------------------------------------------------
7. KEY RESULTS  (Windows build, MSYS2 MinGW-64)
--------------------------------------------------------------------------------
 
  N      Soft Fast (ms)   Soft Slow (ms)   Hard (ms)   Hard/Soft ratio
  -----  ---------------  ---------------  ----------  ---------------
    5         10.31             2.59          38.47          3.7x
   10          5.81             2.70          60.34         10.4x
   25         13.28             2.58          48.36          3.6x
   50         21.78             2.61          94.66          4.3x
  100         34.15             2.67         184.01          5.4x
  250         65.52             2.66         434.07          6.6x
  500        108.74             2.62         736.09          6.8x
 

 
