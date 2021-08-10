#include <assert.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

#include <gmp.h>
#include <openssl/sha.h>

#include "octet_string.h"
#include "finite_field.h"

#include "BLS12_381.h"

#define DST_MAX_LEN 255

#define B_IN_BYTES (256 / 8)
#define R_IN_BYTES 64

/**
 * Parameter L as defined here:
 *      https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-11#section-5.3
 */
#define HTF_PARAM_L 64

#define CEIL(x,y) 1 + ((x - 1) / y)

#define HASH(ctx,dst,src)                     \
    SHA256_Init(&ctx);                        \
    SHA256_Update(&ctx, src->data, src->len); \
    SHA256_Final(dst, &ctx);

octet_string *I2OSP(octet_string *o, int x, int x_len)
{
    octet_string *tmp;

    octet_string_alloc(&tmp, x_len);

    for (int i = 1; i <= x_len && x != 0; i++) {
        tmp->data[x_len - i] = x % 256;
        x /= 256;
    }
    tmp->len = x_len;

    octet_string_reset(o);
    octet_strcpy(o, tmp);
    octet_string_free(tmp);

    return o;
}

void OS2IP(mpz_t x, const octet_string *o)
{
    char *str;

    // TODO(nr): This works fine, but is probably way
    // less efficient than it could be.
    str = octet_string_to_str(o);
    mpz_set_str(x, str, 16);
    free(str);
}

int expand_message_xmd(octet_string *bytes, const char *msg, const char *DST, uint32_t len_in_bytes)
{
    uint32_t ell, DST_len;
    octet_string *DST_prime, *msg_prime, *uniform_bytes,
                 *tmp, *arg;

    SHA256_CTX sha_ctx;
    uint8_t b0[SHA256_DIGEST_LENGTH], bi[SHA256_DIGEST_LENGTH];

    ell = CEIL(len_in_bytes, B_IN_BYTES);
    if (ell > 255) {
        return 1;
    }

    DST_len = strnlen(DST, DST_MAX_LEN);

    octet_string_alloc(&uniform_bytes, SHA256_DIGEST_LENGTH * ell);
    octet_string_alloc(&DST_prime, DST_len + 2);
    octet_string_alloc(&msg_prime, R_IN_BYTES + DST_prime->len + strlen(msg) + 2);
    octet_string_alloc(&tmp, SHA256_DIGEST_LENGTH + DST_len + 1);
    octet_string_alloc(&arg, 0);

    octet_string_appendn(DST_prime, (uint8_t*)DST, DST_len);
    octet_strcat(DST_prime, I2OSP(arg, DST_len, 1));

    octet_strcat(msg_prime, I2OSP(arg, 0, R_IN_BYTES));
    octet_string_appendn(msg_prime, (uint8_t*)msg, strlen(msg));
    octet_strcat(msg_prime, I2OSP(arg, len_in_bytes, 2));
    octet_string_append(msg_prime, 0x0);
    octet_strcat(msg_prime, DST_prime);

    HASH(sha_ctx, b0, msg_prime);

    octet_string_appendn(tmp, b0, SHA256_DIGEST_LENGTH);
    octet_string_append(tmp, 0x1);
    octet_strcat(tmp, DST_prime);

    HASH(sha_ctx, bi, tmp);

    octet_string_appendn(uniform_bytes, bi, SHA256_DIGEST_LENGTH);

    if (ell < 2) goto out;

    for (int i = 2; i <= ell; i++) {
        octet_string_reset(tmp);

        for (int j = 0; j < SHA256_DIGEST_LENGTH; j++) {
            octet_string_append(tmp, b0[j] ^ bi[j]);
        }

        octet_strcat(tmp, I2OSP(arg, i, 1));
        octet_strcat(tmp, DST_prime);

        HASH(sha_ctx, bi, tmp);

        octet_string_appendn(uniform_bytes, bi, SHA256_DIGEST_LENGTH);
    }

out:
    octet_strncpy(bytes, uniform_bytes, len_in_bytes);

    octet_string_free(DST_prime);
    octet_string_free(msg_prime);
    octet_string_free(uniform_bytes);
    octet_string_free(tmp);
    octet_string_free(arg);

    return 0;
}

mpz_t *hash_to_field_fp(const char *msg, const char *DST, uint32_t count)
{
    mpz_t *elems;
    mpz_t p;
    uint32_t len_in_bytes;
    octet_string *uniform_bytes, *tv;

    elems = calloc(count, sizeof(mpz_t));

    mpz_init_set_str(p, BLS12_381_P, 0);

    len_in_bytes = count * HTF_PARAM_L;
    octet_string_alloc(&uniform_bytes, len_in_bytes);
    octet_string_alloc(&tv, HTF_PARAM_L);

    assert(expand_message_xmd(uniform_bytes, msg, DST, len_in_bytes) == 0);

    for (int i = 0; i < count; i++) {
        octet_substr(tv, uniform_bytes, HTF_PARAM_L * i, HTF_PARAM_L);

        mpz_init(elems[i]);
        OS2IP(elems[i], tv);
        mpz_mod(elems[i], elems[i], p);
    }

    octet_string_free(uniform_bytes);
    octet_string_free(tv);
    mpz_clear(p);

    return elems;
}

fp2_elem **hash_to_field_fp2(const char *msg, const char *DST, uint32_t count)
{
    fp2_elem **elems;
    mpz_t p;
    uint32_t len_in_bytes;
    octet_string *uniform_bytes, *tv;

    elems = calloc(count, sizeof(fp2_elem*));

    mpz_init_set_str(p, BLS12_381_P, 0);

    len_in_bytes = count * 2 * HTF_PARAM_L;
    octet_string_alloc(&uniform_bytes, len_in_bytes);
    octet_string_alloc(&tv, HTF_PARAM_L);

    assert(expand_message_xmd(uniform_bytes, msg, DST, len_in_bytes) == 0);

    for (int i = 0; i < count; i++) {
        elems[i] = calloc(1, sizeof(fp2_elem));
        fp2_elem_init(elems[i]);

        octet_substr(tv, uniform_bytes, HTF_PARAM_L * (2 * i), HTF_PARAM_L);
        OS2IP(elems[i]->a, tv);
        mpz_mod(elems[i]->a, elems[i]->a, p);

        octet_substr(tv, uniform_bytes, HTF_PARAM_L * (1 + (2 * i)), HTF_PARAM_L);
        OS2IP(elems[i]->b, tv);
        mpz_mod(elems[i]->b, elems[i]->b, p);
    }

    octet_string_free(uniform_bytes);
    octet_string_free(tv);
    mpz_clear(p);

    return elems;
}

static int sgn0_fp(const mpz_t x)
{
    /**
     * We want to return x (mod 2). This returns non-zero if x is
     * even, and zero otherwise. Thus, we return the negation of
     * the result.
     */
    return !mpz_even_p(x);
}

static int sgn0_fp2(const fp2_elem *x)
{
    int sign_0, zero_0, sign_1;

    sign_0 = sgn0_fp(x->a);
    zero_0 = mpz_sgn(x->a) == 0;
    sign_1 = sgn0_fp(x->b);

    return sign_0 || (zero_0 && sign_1);
}

#define CMOV(x,a,b,c)  \
    if (!c) {          \
        mpz_set(x, a); \
    } else {           \
        mpz_set(x, b); \
    }

#define CMOV2(x,a,b,c)      \
    if (!c) {               \
        fp2_elem_set(x, a); \
    } else {                \
        fp2_elem_set(x, b); \
    }

/**
 * Implementation of the simplified SWU method, specifically for curves
 * where order of the finite field is equivalent to 3 mod 4. This is the
 * case for the curve E'/F_p which is isogenoes to E: y^2 = x^3 + 4. This
 * is required for G1 hashing.
 *
 * https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-11#appendix-G.2.1
 */
#define SWU_3MOD4_C1 "0x680447a8e5ff9a692c6e9ed90d2eb35d91dd2e13ce144afd9cc34a83dac3d8907aaffffac54ffffee7fbfffffffeaaa"
#define SWU_3MOD4_C2 "0x3d689d1e0e762cef9f2bec6130316806b4c80eda6fc10ce77ae83eab1ea8b8b8a407c9c6db195e06f2dbeabc2baeff5"
void map_to_curve_simple_swu_3mod4(mpz_t xn, mpz_t xd, mpz_t y, const mpz_t u)
{
    mpz_t c1, c2, A, B, p,
          tv1, tv2, tv3, tv4,
          x1n, x2n, gxd, gx1,
          y1, y2, tmp;
    int e1, e2, e3;

    mpz_init_set_str(c1, SWU_3MOD4_C1, 0);
    mpz_init_set_str(c2, SWU_3MOD4_C2, 0);
    mpz_init_set_str(p, BLS12_381_P, 0);
    mpz_init_set_str(A, BLS12_381_ISOGENY_A, 0);
    mpz_init_set_str(B, BLS12_381_ISOGENY_B, 0);

    mpz_inits(tv1, tv2, tv3, tv4,
              x1n, x2n, gxd, gx1,
              y1, y2, tmp, NULL);

    fp_mul(tv1, u, u);             // 1.  tv1 = u^2
    fp_mul_ui(tv3, tv1, 11);       // 2.  tv3 = Z * tv1
    fp_mul(tv2, tv3, tv3);         // 3.  tv2 = tv3^2
    fp_add(xd, tv2, tv3);          // 4.   xd = tv2 + tv3
    fp_add_ui(x1n, xd, 1);         // 5.  x1n = xd + 1
    fp_mul(x1n, x1n, B);           // 6.  x1n = x1n * B
    fp_mul(xd, xd, A);             // 7.   xd = -A * xd
    fp_negate(xd, xd);
    e1 = mpz_cmp_ui(xd, 0) == 0;   // 8.   e1 = xd == 0
    fp_mul_ui(tmp, A, 11);         // 9.   xd = CMOV(xd, Z * A, e1)  # If xd == 0, set xd = Z * A
    CMOV(xd, xd, tmp, e1);
    fp_mul(tv2, xd, xd);           // 10. tv2 = xd^2
    fp_mul(gxd, tv2, xd);          // 11. gxd = tv2 * xd             # gxd == xd^3
    fp_mul(tv2, A, tv2);           // 12. tv2 = A * tv2
    fp_mul(gx1, x1n, x1n);         // 13. gx1 = x1n^2
    fp_add(gx1, gx1, tv2);         // 14. gx1 = gx1 + tv2            # x1n^2 + A * xd^2
    fp_mul(gx1, gx1, x1n);         // 15. gx1 = gx1 * x1n            # x1n^3 + A * x1n * xd^2
    fp_mul(tv2, B, gxd);           // 16. tv2 = B * gxd
    fp_add(gx1, gx1, tv2);         // 17. gx1 = gx1 + tv2            # x1n^3 + A * x1n * xd^2 + B * xd^3
    fp_mul(tv4, gxd, gxd);         // 18. tv4 = gxd^2
    fp_mul(tv2, gx1, gxd);         // 19. tv2 = gx1 * gxd
    fp_mul(tv4, tv4, tv2);         // 20. tv4 = tv4 * tv2            # gx1 * gxd^3
    mpz_powm(y1, tv4, c1, p);      // 21.  y1 = tv4^c1               # (gx1 * gxd^3)^((q - 3) / 4)
    fp_mul(y1, y1, tv2);           // 22.  y1 = y1 * tv2             # gx1 * gxd * (gx1 * gxd^3)^((q - 3) / 4)
    fp_mul(x2n, tv3, x1n);         // 23. x2n = tv3 * x1n            # x2 = x2n / xd = Z * u^2 * x1n / xd
    fp_mul(y2, y1, c2);            // 24.  y2 = y1 * c2              # y2 = y1 * sqrt(-Z^3)
    fp_mul(y2, y2, tv1);           // 25.  y2 = y2 * tv1
    fp_mul(y2, y2, u);             // 26.  y2 = y2 * u
    fp_mul(tv2, y1, y1);           // 27. tv2 = y1^2
    fp_mul(tv2, tv2, gxd);         // 28. tv2 = tv2 * gxd
    e2 = mpz_cmp(tv2, gx1) == 0;   // 29.  e2 = tv2 == gx1
    CMOV(xn, x2n, x1n, e2);        // 30.  xn = CMOV(x2n, x1n, e2)   # If e2, x = x1, else x = x2
    CMOV(y, y2, y1, e2);           // 31.   y = CMOV(y2, y1, e2)     # If e2, y = y1, else y = y2
    e3 = sgn0_fp(u) == sgn0_fp(y); // 32.  e3 = sgn0(u) == sgn0(y)   # Fix sign of y
    fp_negate(tmp, y);             // 33.   y = CMOV(-y, y, e3)
    CMOV(y, tmp, y, e3);
                                   // 34. return (xn, xd, y, 1)

    mpz_clears(c1, c2, A, B, p,
               tv1, tv2, tv3, tv4,
               x1n, x2n, gxd, gx1,
               y1, y2, tmp, NULL);
}

/**
 * Implementation of the simplified SWU method, specifically for curves
 * where order of the finite field is equivalent to 9 mod 16. This is the
 * case for the curve E'/F_p^2 which is isogenoes to E: y^2 = x^3 + 4(u + 1).
 * This is required for G2 hashing.
 *
 * https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-11#appendix-G.2.2
 */
#define SWU_9MOD16_C1 "0x2a437a4b8c35fc74bd278eaa22f25e9e2dc90e50e7046b466e59e4"\
                      "9349e8bd050a62cfd16ddca6ef53149330978ef011d68619c86185c7"\
                      "b292e85a87091a04966bf91ed3e71b743162c338362113cfd7ced6b1"\
                      "d76382eab26aa00001c718e3"

#define SWU_9MOD16_C2_A "0"
#define SWU_9MOD16_C2_B "1"

#define SWU_9MOD16_C3_A "0x135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60ef396"\
                        "489f61eb45e304466cf3e67fa0af1ee7b04121bdea2"
#define SWU_9MOD16_C3_B "0x06af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76"\
                        "e17009241c5ee67992f72ec05f4c81084fbede3cc09"

#define SWU_9MOD16_C4_A "0x0699be3b8c6870965e5bf892ad5d2cc7b0e85a117402dfd83b7f4"\
                        "a947e02d978498255a2aaec0ac627b5afbdf1bf1c90"
#define SWU_9MOD16_C4_B "0x08157cd83046453f5dd0972b6e3949e4288020b5b8a9cc99ca07e"\
                        "27089a2ce2436d965026adad3ef7baba37f2183e9b5"

#define SWU_9MOD16_C5_A "0x0f5d0d63d2797471e6d39f306cc0dc0ab85de3bd9f39ce46f3649"\
                        "ac0de9e844417cc8de88716c1fd323fa68040801aea"
#define SWU_9MOD16_C5_B "0x0ab1c2ffdd6c253ca155231eb3e71ba044fd562f6f72bc5bad5ec"\
                        "46a0b7a3b0247cf08ce6c6317f40edbc653a72dee17"

void map_to_curve_simple_swu_9mod16(fp2_elem *xn, fp2_elem *xd, fp2_elem *y, const fp2_elem *u)
{
    fp2_elem tv1, tv2, tv3, tv4, tv5,
             x1n, gxd, gx1, gx2,
             c2, c3, c4, c5,
             Z, A, B, tmp;
    mpz_t c1;
    int e;

    fp2_elem_init(&tv1);
    fp2_elem_init(&tv2);
    fp2_elem_init(&tv3);
    fp2_elem_init(&tv4);
    fp2_elem_init(&tv5);
    fp2_elem_init(&x1n);
    fp2_elem_init(&gxd);
    fp2_elem_init(&gx1);
    fp2_elem_init(&gx2);
    fp2_elem_init(&Z);
    fp2_elem_init(&tmp);

    // Z  = -(2 + I)
    fp2_elem_set_si(&Z, -2, -1);

    mpz_init_set_str(c1, SWU_9MOD16_C1, 0);
    fp2_elem_from_str(&c2, SWU_9MOD16_C2_A, SWU_9MOD16_C2_B);
    fp2_elem_from_str(&c3, SWU_9MOD16_C3_A, SWU_9MOD16_C3_B);
    fp2_elem_from_str(&c4, SWU_9MOD16_C4_A, SWU_9MOD16_C4_B);
    fp2_elem_from_str(&c5, SWU_9MOD16_C5_A, SWU_9MOD16_C5_B);

    fp2_elem_from_str(&A,
            BLS12_381_TWIST_ISOGENY_A_X,
            BLS12_381_TWIST_ISOGENY_A_Y
    );

    fp2_elem_from_str(&B,
            BLS12_381_TWIST_ISOGENY_B_X,
            BLS12_381_TWIST_ISOGENY_B_Y
    );

    fp2_square(&tv1, u);             // 1.  tv1 = u^2
    fp2_mul(&tv3, &Z, &tv1);         // 2.  tv3 = Z * tv1
    fp2_square(&tv5, &tv3);          // 3.  tv5 = tv3^2
    fp2_add(xd, &tv5, &tv3);         // 4.   xd = tv5 + tv3
    fp2_elem_set_si(&tmp, 1, 0);     // 5.  x1n = xd + 1
    fp2_add(&x1n, xd, &tmp);
    fp2_mul(&x1n, &x1n, &B);         // 6.  x1n = x1n * B
    fp2_mul(xd, &A, xd);             // 7.   xd = -A * xd
    fp2_negate(xd, xd);
    fp2_elem_set_si(&tmp, 0, 0);     // 8.   e1 = xd == 0
    e = fp2_equal(xd, &tmp);
    fp2_mul(&tmp, &Z, &A);           // 9.   xd = CMOV(xd, Z * A, e1)   # If xd == 0, set xd = Z * A
    CMOV2(xd, xd, &tmp, e);
    fp2_square(&tv2, xd);            // 10. tv2 = xd^2
    fp2_mul(&gxd, &tv2, xd);         // 11. gxd = tv2 * xd              # gxd == xd^3
    fp2_mul(&tv2, &A, &tv2);         // 12. tv2 = A * tv2
    fp2_square(&gx1, &x1n);          // 13. gx1 = x1n^2
    fp2_add(&gx1, &gx1, &tv2);       // 14. gx1 = gx1 + tv2             # x1n^2 + A * xd^2
    fp2_mul(&gx1, &gx1, &x1n);       // 15. gx1 = gx1 * x1n             # x1n^3 + A * x1n * xd^2
    fp2_mul(&tv2, &B, &gxd);         // 16. tv2 = B * gxd
    fp2_add(&gx1, &gx1, &tv2);       // 17. gx1 = gx1 + tv2             # x1n^3 + A * x1n * xd^2 + B * xd^3
    fp2_square(&tv4, &gxd);          // 18. tv4 = gxd^2
    fp2_mul(&tv2, &tv4, &gxd);       // 19. tv2 = tv4 * gxd             # gxd^3
    fp2_square(&tv4, &tv4);          // 20. tv4 = tv4^2                 # gxd^4
    fp2_mul(&tv2, &tv2, &tv4);       // 21. tv2 = tv2 * tv4             # gxd^7
    fp2_mul(&tv2, &tv2, &gx1);       // 22. tv2 = tv2 * gx1             # gx1 * gxd^7
    fp2_square(&tv4, &tv4);          // 23. tv4 = tv4^2                 # gxd^8
    fp2_mul(&tv4, &tv2, &tv4);       // 24. tv4 = tv2 * tv4             # gx1 * gxd^15
    fp2_pow(y, &tv4, c1);            // 25.   y = tv4^c1                # (gx1 * gxd^15)^((q - 9) / 16)
    fp2_mul(y, y, &tv2);             // 26.   y = y * tv2               # This is almost sqrt(gx1)
    fp2_mul(&tv4, y, &c2);           // 27. tv4 = y * c2                # check the four possible sqrts
    fp2_square(&tv2, &tv4);          // 28. tv2 = tv4^2
    fp2_mul(&tv2, &tv2, &gxd);       // 29. tv2 = tv2 * gxd
    e = fp2_equal(&tv2, &gx1);       // 30.  e2 = tv2 == gx1
    CMOV2(y, y, &tv4, e);            // 31.   y = CMOV(y, tv4, e2)
    fp2_mul(&tv4, y, &c3);           // 32. tv4 = y * c3
    fp2_square(&tv2, &tv4);          // 33. tv2 = tv4^2
    fp2_mul(&tv2, &tv2, &gxd);       // 34. tv2 = tv2 * gxd
    e = fp2_equal(&tv2, &gx1);       // 35.  e3 = tv2 == gx1
    CMOV2(y, y, &tv4, e);            // 36.   y = CMOV(y, tv4, e3)
    fp2_mul(&tv4, &tv4, &c2);        // 37. tv4 = tv4 * c2
    fp2_square(&tv2, &tv4);          // 38. tv2 = tv4^2
    fp2_mul(&tv2, &tv2, &gxd);       // 39. tv2 = tv2 * gxd
    e = fp2_equal(&tv2, &gx1);       // 40.  e4 = tv2 == gx1
    CMOV2(y, y, &tv4, e);            // 41.   y = CMOV(y, tv4, e4)      # if x1 is square, this is its sqrt
    fp2_mul(&gx2, &gx1, &tv5);       // 42. gx2 = gx1 * tv5
    fp2_mul(&gx2, &gx2, &tv3);       // 43. gx2 = gx2 * tv3             # gx2 = gx1 * Z^3 * u^6
    fp2_mul(&tv5, y, &tv1);          // 44. tv5 = y * tv1
    fp2_mul(&tv5, &tv5, u);          // 45. tv5 = tv5 * u               # This is almost sqrt(gx2)
    fp2_mul(&tv1, &tv5, &c4);        // 46. tv1 = tv5 * c4              # check the four possible sqrts
    fp2_mul(&tv4, &tv1, &c2);        // 47. tv4 = tv1 * c2
    fp2_square(&tv2, &tv4);          // 48. tv2 = tv4^2
    fp2_mul(&tv2, &tv2, &gxd);       // 49. tv2 = tv2 * gxd
    e = fp2_equal(&tv2, &gx2);       // 50.  e5 = tv2 == gx2
    CMOV2(&tv1, &tv1, &tv4, e);      // 51. tv1 = CMOV(tv1, tv4, e5)
    fp2_mul(&tv4, &tv5, &c5);        // 52. tv4 = tv5 * c5
    fp2_square(&tv2, &tv4);          // 53. tv2 = tv4^2
    fp2_mul(&tv2, &tv2, &gxd);       // 54. tv2 = tv2 * gxd
    e = fp2_equal(&tv2, &gx2);       // 55.  e6 = tv2 == gx2
    CMOV2(&tv1, &tv1, &tv4, e);      // 56. tv1 = CMOV(tv1, tv4, e6)
    fp2_mul(&tv4, &tv4, &c2);        // 57. tv4 = tv4 * c2
    fp2_square(&tv2, &tv4);          // 58. tv2 = tv4^2
    fp2_mul(&tv2, &tv2, &gxd);       // 59. tv2 = tv2 * gxd
    e = fp2_equal(&tv2, &gx2);       // 60.  e7 = tv2 == gx2
    CMOV2(&tv1, &tv1, &tv4, e);      // 61. tv1 = CMOV(tv1, tv4, e7)
    fp2_square(&tv2, y);             // 62. tv2 = y^2
    fp2_mul(&tv2, &tv2, &gxd);       // 63. tv2 = tv2 * gxd
    e = fp2_equal(&tv2, &gx1);       // 64.  e8 = tv2 == gx1
    CMOV2(y, &tv1, y, e);            // 65.   y = CMOV(tv1, y, e8)      # choose correct y-coordinate
    fp2_mul(&tv2, &tv3, &x1n);       // 66. tv2 = tv3 * x1n             # x2n = x2n / xd = Z * u^2 * x1n / xd
    CMOV2(xn, &tv2, &x1n, e);        // 67.  xn = CMOV(tv2, x1n, e8)    # choose correct x-coordinate
    e = sgn0_fp2(u) == sgn0_fp2(y);  // 68.  e9 = sgn0(u) == sgn0(y)    # Fix sign of y
    fp2_negate(&tmp, y);             // 69.   y = CMOV(-y, y, e9)
    CMOV2(y, &tmp, y, e);
                                     // 70. return (xn, xd, y, 1)

    fp2_elem_free(&tv1);
    fp2_elem_free(&tv2);
    fp2_elem_free(&tv3);
    fp2_elem_free(&tv4);
    fp2_elem_free(&tv5);
    fp2_elem_free(&x1n);
    fp2_elem_free(&gxd);
    fp2_elem_free(&gx1);
    fp2_elem_free(&gx2);
    fp2_elem_free(&Z);
    fp2_elem_free(&tmp);
    fp2_elem_free(&c2);
    fp2_elem_free(&c3);
    fp2_elem_free(&c4);
    fp2_elem_free(&c5);
    fp2_elem_free(&A);
    fp2_elem_free(&B);
    mpz_clear(c1);
}

/**
 * The constants used in evaluating the 11-isogeny from E'/F to E/F.
 */
const char *iso_G1_k1[12] = {
    "0x11a05f2b1e833340b809101dd99815856b303e88a2d7005ff2627b56cdb4e2c85610c2d5f2e62d6eaeac1662734649b7",
    "0x17294ed3e943ab2f0588bab22147a81c7c17e75b2f6a8417f565e33c70d1e86b4838f2a6f318c356e834eef1b3cb83bb",
    "0xd54005db97678ec1d1048c5d10a9a1bce032473295983e56878e501ec68e25c958c3e3d2a09729fe0179f9dac9edcb0",
    "0x1778e7166fcc6db74e0609d307e55412d7f5e4656a8dbf25f1b33289f1b330835336e25ce3107193c5b388641d9b6861",
    "0xe99726a3199f4436642b4b3e4118e5499db995a1257fb3f086eeb65982fac18985a286f301e77c451154ce9ac8895d9",
    "0x1630c3250d7313ff01d1201bf7a74ab5db3cb17dd952799b9ed3ab9097e68f90a0870d2dcae73d19cd13c1c66f652983",
    "0xd6ed6553fe44d296a3726c38ae652bfb11586264f0f8ce19008e218f9c86b2a8da25128c1052ecaddd7f225a139ed84",
    "0x17b81e7701abdbe2e8743884d1117e53356de5ab275b4db1a682c62ef0f2753339b7c8f8c8f475af9ccb5618e3f0c88e",
    "0x80d3cf1f9a78fc47b90b33563be990dc43b756ce79f5574a2c596c928c5d1de4fa295f296b74e956d71986a8497e317",
    "0x169b1f8e1bcfa7c42e0c37515d138f22dd2ecb803a0c5c99676314baf4bb1b7fa3190b2edc0327797f241067be390c9e",
    "0x10321da079ce07e272d8ec09d2565b0dfa7dccdde6787f96d50af36003b14866f69b771f8c285decca67df3f1605fb7b",
    "0x6e08c248e260e70bd1e962381edee3d31d79d7e22c837bc23c0bf1bc24c6b68c24b1b80b64d391fa9c8ba2e8ba2d229"
};

const char *iso_G1_k2[10] = {
    "0x8ca8d548cff19ae18b2e62f4bd3fa6f01d5ef4ba35b48ba9c9588617fc8ac62b558d681be343df8993cf9fa40d21b1c",
    "0x12561a5deb559c4348b4711298e536367041e8ca0cf0800c0126c2588c48bf5713daa8846cb026e9e5c8276ec82b3bff",
    "0xb2962fe57a3225e8137e629bff2991f6f89416f5a718cd1fca64e00b11aceacd6a3d0967c94fedcfcc239ba5cb83e19",
    "0x3425581a58ae2fec83aafef7c40eb545b08243f16b1655154cca8abc28d6fd04976d5243eecf5c4130de8938dc62cd8",
    "0x13a8e162022914a80a6f1d5f43e7a07dffdfc759a12062bb8d6b44e833b306da9bd29ba81f35781d539d395b3532a21e",
    "0xe7355f8e4e667b955390f7f0506c6e9395735e9ce9cad4d0a43bcef24b8982f7400d24bc4228f11c02df9a29f6304a5",
    "0x772caacf16936190f3e0c63e0596721570f5799af53a1894e2e073062aede9cea73b3538f0de06cec2574496ee84a3a",
    "0x14a7ac2a9d64a8b230b3f5b074cf01996e7f63c21bca68a81996e1cdf9822c580fa5b9489d11e2d311f7d99bbdcc5a5e",
    "0xa10ecf6ada54f825e920b3dafc7a3cce07f8d1d7161366b74100da67f39883503826692abba43704776ec3a79a1d641",
    "0x95fc13ab9e92ad4476d6e3eb3a56680f682b4ee96f7d03776df533978f31c1593174e4b4b7865002d6384d168ecdd0a"
};

const char *iso_G1_k3[16] = {
    "0x90d97c81ba24ee0259d1f094980dcfa11ad138e48a869522b52af6c956543d3cd0c7aee9b3ba3c2be9845719707bb33",
    "0x134996a104ee5811d51036d776fb46831223e96c254f383d0f906343eb67ad34d6c56711962fa8bfe097e75a2e41c696",
    "0xcc786baa966e66f4a384c86a3b49942552e2d658a31ce2c344be4b91400da7d26d521628b00523b8dfe240c72de1f6",
    "0x1f86376e8981c217898751ad8746757d42aa7b90eeb791c09e4a3ec03251cf9de405aba9ec61deca6355c77b0e5f4cb",
    "0x8cc03fdefe0ff135caf4fe2a21529c4195536fbe3ce50b879833fd221351adc2ee7f8dc099040a841b6daecf2e8fedb",
    "0x16603fca40634b6a2211e11db8f0a6a074a7d0d4afadb7bd76505c3d3ad5544e203f6326c95a807299b23ab13633a5f0",
    "0x4ab0b9bcfac1bbcb2c977d027796b3ce75bb8ca2be184cb5231413c4d634f3747a87ac2460f415ec961f8855fe9d6f2",
    "0x987c8d5333ab86fde9926bd2ca6c674170a05bfe3bdd81ffd038da6c26c842642f64550fedfe935a15e4ca31870fb29",
    "0x9fc4018bd96684be88c9e221e4da1bb8f3abd16679dc26c1e8b6e6a1f20cabe69d65201c78607a360370e577bdba587",
    "0xe1bba7a1186bdb5223abde7ada14a23c42a0ca7915af6fe06985e7ed1e4d43b9b3f7055dd4eba6f2bafaaebca731c30",
    "0x19713e47937cd1be0dfd0b8f1d43fb93cd2fcbcb6caf493fd1183e416389e61031bf3a5cce3fbafce813711ad011c132",
    "0x18b46a908f36f6deb918c143fed2edcc523559b8aaf0c2462e6bfe7f911f643249d9cdf41b44d606ce07c8a4d0074d8e",
    "0xb182cac101b9399d155096004f53f447aa7b12a3426b08ec02710e807b4633f06c851c1919211f20d4c04f00b971ef8",
    "0x245a394ad1eca9b72fc00ae7be315dc757b3b080d4c158013e6632d3c40659cc6cf90ad1c232a6442d9d3f5db980133",
    "0x5c129645e44cf1102a159f748c4a3fc5e673d81d7e86568d9ab0f5d396a7ce46ba1049b6579afb7866b1e715475224b",
    "0x15e6be4e990f03ce4ea50b3b42df2eb5cb181d8f84965a3957add4fa95af01b2b665027efec01c7704b456be69c8b604"
};

const char *iso_G1_k4[15] = {
    "0x16112c4c3a9c98b252181140fad0eae9601a6de578980be6eec3232b5be72e7a07f3688ef60c206d01479253b03663c1",
    "0x1962d75c2381201e1a0cbd6c43c348b885c84ff731c4d59ca4a10356f453e01f78a4260763529e3532f6102c2e49a03d",
    "0x58df3306640da276faaae7d6e8eb15778c4855551ae7f310c35a5dd279cd2eca6757cd636f96f891e2538b53dbf67f2",
    "0x16b7d288798e5395f20d23bf89edb4d1d115c5dbddbcd30e123da489e726af41727364f2c28297ada8d26d98445f5416",
    "0xbe0e079545f43e4b00cc912f8228ddcc6d19c9f0f69bbb0542eda0fc9dec916a20b15dc0fd2ededda39142311a5001d",
    "0x8d9e5297186db2d9fb266eaac783182b70152c65550d881c5ecd87b6f0f5a6449f38db9dfa9cce202c6477faaf9b7ac",
    "0x166007c08a99db2fc3ba8734ace9824b5eecfdfa8d0cf8ef5dd365bc400a0051d5fa9c01a58b1fb93d1a1399126a775c",
    "0x16a3ef08be3ea7ea03bcddfabba6ff6ee5a4375efa1f4fd7feb34fd206357132b920f5b00801dee460ee415a15812ed9",
    "0x1866c8ed336c61231a1be54fd1d74cc4f9fb0ce4c6af5920abc5750c4bf39b4852cfe2f7bb9248836b233d9d55535d4a",
    "0x167a55cda70a6e1cea820597d94a84903216f763e13d87bb5308592e7ea7d4fbc7385ea3d529b35e346ef48bb8913f55",
    "0x4d2f259eea405bd48f010a01ad2911d9c6dd039bb61a6290e591b36e636a5c871a5c29f4f83060400f8b49cba8f6aa8",
    "0xaccbb67481d033ff5852c1e48c50c477f94ff8aefce42d28c0f9a88cea7913516f968986f7ebbea9684b529e2561092",
    "0xad6b9514c767fe3c3613144b45f1496543346d98adf02267d5ceef9a00d9b8693000763e3b90ac11e99b138573345cc",
    "0x2660400eb2e4f3b628bdd0d53cd76f2bf565b94e72927c1cb748df27942480e420517bd8714cc80d1fadc1326ed06f7",
    "0xe0fa1d816ddc03e6b24255e0d7819c171c40f65e273b853324efcd6356caa205ca2f570f13497804415473a1d634b8f"
};

void iso_map_G1(mpz_t x, mpz_t y, const mpz_t x_prime, const mpz_t y_prime)
{
    mpz_t xn, xd, yn, yd,
          p, tmp, k;

    mpz_inits(xn, xd, yn, yd, tmp, k, NULL);
    mpz_init_set_str(p, BLS12_381_P, 0);

    /**
     * TODO(nr): As a future optimization, Horner's rule can
     *           be used here. This literal approach works
     *           fine, however.
     */

    // x_num
    for (int i = 0; i < 12; i++) {
        mpz_set_str(k, iso_G1_k1[i], 0);

        mpz_powm_ui(tmp, x_prime, i, p);
        fp_mul(tmp, tmp, k);

        fp_add(xn, xn, tmp);
    }

    // x_den
    mpz_powm_ui(xd, x_prime, 10, p);
    for (int i = 0; i < 10; i++) {
        mpz_set_str(k, iso_G1_k2[i], 0);

        mpz_powm_ui(tmp, x_prime, i, p);
        fp_mul(tmp, tmp, k);

        fp_add(xd, xd, tmp);
    }

    // y_num
    for (int i = 0; i < 16; i++) {
        mpz_set_str(k, iso_G1_k3[i], 0);

        mpz_powm_ui(tmp, x_prime, i, p);
        fp_mul(tmp, tmp, k);

        fp_add(yn, yn, tmp);
    }

    // y_den
    mpz_powm_ui(yd, x_prime, 15, p);
    for (int i = 0; i < 15; i++) {
        mpz_set_str(k, iso_G1_k4[i], 0);

        mpz_powm_ui(tmp, x_prime, i, p);
        fp_mul(tmp, tmp, k);

        fp_add(yd, yd, tmp);
    }

    fp_inv(x, xd);
    fp_mul(x, x, xn);

    fp_inv(y, yd);
    fp_mul(yn, y_prime, yn);
    fp_mul(y, yn, y);

    mpz_clears(xn, xd, yn, yd, tmp, p, k, NULL);
}

/**
 * TODO(nr): It may be more convenient to leave things in
 *           projective coordinates, but let's get things
 *           working for now.
 */
void map_to_curve_G1(G1_elem_affine *P, const mpz_t u)
{
    mpz_t x, xn, xd, y;

    mpz_inits(x, xn, xd, y, NULL);

    map_to_curve_simple_swu_3mod4(xn, xd, y, u);

    fp_inv(x, xd);
    fp_mul(x, x, xn);

    iso_map_G1(P->x, P->y, x, y);
    P->infinity = false;

    mpz_clears(x, xn, xd, y, NULL);
}

void clear_cofactor_G1(G1_elem_affine *P)
{
    mpz_t h_eff;

    mpz_init_set_str(h_eff, BLS12_381_G1_H_EFF, 0);
    G1_mul_scalar(P, P, h_eff);
    mpz_clear(h_eff);
}

void BLS12_381_hash_to_G1(G1_elem_affine *P, const uint8_t *bytes, const uint8_t *DST)
{
    mpz_t *u;
    G1_elem_affine Q0, Q1;
    G1_elem_proj R, Q0_proj, Q1_proj;

    G1_identity_init_affine(&Q0);
    G1_identity_init_affine(&Q1);
    G1_identity_init_proj(&R);
    G1_identity_init_proj(&Q0_proj);
    G1_identity_init_proj(&Q1_proj);

    u = hash_to_field_fp(bytes, DST, 2);

    map_to_curve_G1(&Q0, u[0]);
    map_to_curve_G1(&Q1, u[1]);

    G1_affine2proj(&Q0_proj, &Q0);
    G1_affine2proj(&Q1_proj, &Q1);

    G1_add_proj(&R, &Q0_proj, &Q1_proj);

    G1_proj2affine(P, &R);
    clear_cofactor_G1(P);

    G1_elem_free_affine(&Q0);
    G1_elem_free_affine(&Q1);
    G1_elem_free_proj(&R);
    G1_elem_free_proj(&Q0_proj);
    G1_elem_free_proj(&Q1_proj);
}
