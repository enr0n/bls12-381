#include <assert.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

#include <stdio.h>

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
