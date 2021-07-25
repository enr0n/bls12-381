/**
 * Finite field arithmetic for BLS12-381.
 *
 * The p<k> in a symbol indicates that type/function is assocatiated with
 * the finite field F_p^k. For k > 1, these fields are defined as extensions
 * of the base field F_p as such:
 *
 *      GF(p^2)  = GF(p)[u] / (u^2 + 1)
 *      GF(p^6)  = GF(p^2)[v] / (v^3 - u - 1)
 *      GF(p^12) = GF(p^6)[w] / (w^2 - v)
 *
 * When a comment refers "Algorithm X", it is referring to that algorithm
 * as defined in "Guide to Pairing-based Cryptography."
 */
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>

#include "finite_field.h"

#include "BLS12_381.h"

/**
 * The prime P which is the characteristic of the finite fields.
 *
 * This variable must be initialized with fp_params_init(), and
 * free'd with fp_params_free(). Each of these functions must
 * only be called once.
 *
 * If we could define this as static const, we would. This variable
 * is never modified.
 */
static mpz_t g_BLS12_381_P;

void fp_params_init()
{
    mpz_set_str(g_BLS12_381_P, BLS12_381_P, 0);
}

void fp_params_free()
{
    mpz_clear(g_BLS12_381_P);
}

bool fp_equiv(const mpz_t x, const mpz_t y)
{
    return mpz_congruent_p(x, y, g_BLS12_381_P);
}

void fp_add(mpz_t x, const mpz_t y, const mpz_t z)
{
    mpz_add(x, y, z);
    mpz_mod(x, x, g_BLS12_381_P);
}

void fp_sub(mpz_t x, const mpz_t y, const mpz_t z)
{
    mpz_sub(x, y, z);
    mpz_mod(x, x, g_BLS12_381_P);
}

void fp_mul(mpz_t x, const mpz_t y, const mpz_t z)
{
    mpz_mul(x, y, z);
    mpz_mod(x, x, g_BLS12_381_P);
}

void fp_inv(mpz_t x, const mpz_t y)
{
    // We are in a finite field, the inverse MUST exist.
    assert(mpz_invert(x, y, g_BLS12_381_P) != 0);
}

void fp2_elem_init(fp2_elem *e)
{
    mpz_inits(e->a, e->b, NULL);
}

void fp2_elem_set(fp2_elem *e1, const fp2_elem *e2)
{
    mpz_set(e1->a, e2->a);
    mpz_set(e1->b, e2->b);
}

void fp2_elem_set_si(fp2_elem *e, signed long int a, signed long int b)
{
    mpz_set_si(e->a, a);
    mpz_set_si(e->b, b);
}

void fp2_elem_from_str(fp2_elem *e, const char *a, const char *b)
{
    mpz_init_set_str(e->a, a, 0);
    mpz_init_set_str(e->b, b, 0);
}

void fp2_elem_free(fp2_elem *e)
{
    mpz_clears(e->a, e->b, NULL);
}

char *fp2_elem_get_str(const fp2_elem *e)
{
    char *str;

    str = calloc(1, sizeof(char));
    sprintf(str, "%s + %s * u",
            mpz_get_str(NULL, 16, e->a),
            mpz_get_str(NULL, 16, e->b));

    return str;
}

void fp2_add(fp2_elem *x, const fp2_elem *y, const fp2_elem *z)
{
    /* Algorithm 5.15 */
    fp_add(x->a, y->a, z->a);
    fp_add(x->b, y->b, z->b);
}

void fp2_sub(fp2_elem *x, const fp2_elem *y, const fp2_elem *z)
{
    fp_sub(x->a, y->a, z->a);
    fp_sub(x->b, y->b, z->b);
}

void fp2_mul(fp2_elem *x, const fp2_elem *y, const fp2_elem *z)
{
    /* Algorithm 5.16 */
    mpz_t v0, v1, tmp, ra, rb;

    mpz_inits(v0, v1, tmp, ra, rb, NULL);

    fp_mul(v0, y->a, z->a);
    fp_mul(v1, y->b, z->b);

    fp_sub(ra, v0, v1);

    fp_add(rb, y->a, y->b);
    fp_add(tmp, z->a, z->b);
    fp_mul(rb, rb, tmp);
    fp_sub(rb, rb, v0);
    fp_sub(rb, rb, v1);

    mpz_set(x->a, ra);
    mpz_set(x->b, rb);

    mpz_clears(v0, v1, tmp, ra, rb, NULL);
}

void fp2_square(fp2_elem *x, const fp2_elem *y)
{
    /* Algorithm 5.17 */
    mpz_t u, v, w, ra, rb;
    mpz_inits(u, v, w, ra, rb, NULL);

    fp_sub(u, y->a, y->b);

    mpz_neg(w, y->b);
    fp_sub(w, y->a, w);
    fp_mul(v, y->a, y->b);
    fp_mul(u, u, w);
    fp_add(u, u, v);

    fp_add(rb, v, v);

    mpz_neg(v, v);
    fp_add(ra, u, v);

    mpz_set(x->a, ra);
    mpz_set(x->b, rb);

    mpz_clears(u, v, w, ra, rb, NULL);
}

void fp2_inv(fp2_elem *x, const fp2_elem *y)
{
    /* Algorithm 5.19 */
    mpz_t v, w, exp, ra, rb;
    mpz_inits(v, w, exp, ra, rb, NULL);

    mpz_set_ui(exp, 2);

    mpz_powm(v, y->a, exp, g_BLS12_381_P);
    mpz_powm(w, y->b, exp, g_BLS12_381_P);
    mpz_neg(w, w);
    fp_sub(v, v, w);
    fp_inv(w, v);

    fp_mul(ra, y->a, w);

    mpz_neg(rb, y->b);
    fp_mul(rb, rb, w);

    mpz_set(x->a, ra);
    mpz_set(x->b, rb);

    mpz_clears(v, w, exp, ra, rb, NULL);
}

void fp2_negate(fp2_elem *x, const fp2_elem *y)
{
    mpz_neg(x->a, y->a);
    mpz_mod(x->a, x->a, g_BLS12_381_P);
    mpz_neg(x->b, y->b);
    mpz_mod(x->b, x->b, g_BLS12_381_P);
}

void fp2_conjugate(fp2_elem *x, const fp2_elem *y)
{
    mpz_set(x->a, y->a);
    mpz_neg(x->b, y->b);
    mpz_mod(x->b, x->b, g_BLS12_381_P);
}

void fp2_mul_nonresidue(fp2_elem *x, const fp2_elem *y)
{
    /**
     * Need to multiply an element in F_p^2 by quadratic/cubic non-residue
     * of this field that is used to construct F_p^6. In this case that is
     * u + 1.
     *
     * This amounts to (a + bu)(u + 1) = au + a + bu^2 + bu
     *                                 = au + a - b + bu      // since u^2 = beta = -1
     *                                 = (a - b) + (a + b)u
     */
    fp2_elem r;

    fp2_elem_init(&r);

    fp_sub(r.a, y->a, y->b);
    fp_add(r.b, y->a, y->b);

    fp2_elem_set(x, &r);

    fp2_elem_free(&r);
}

void fp2_mul_scalar(fp2_elem *x, const fp2_elem *y, const mpz_t m)
{
   fp_mul(x->a, y->a, m);
   fp_mul(x->b, y->b, m);
}

bool fp2_equal(const fp2_elem *e1, const fp2_elem *e2)
{
    return (mpz_cmp(e1->a, e2->a) == 0) && (mpz_cmp(e1->b, e2->b) == 0);
}

void fp6_elem_init(fp6_elem *e)
{
    e->a = calloc(1, sizeof(fp2_elem));
    e->b = calloc(1, sizeof(fp2_elem));
    e->c = calloc(1, sizeof(fp2_elem));

    fp2_elem_init(e->a);
    fp2_elem_init(e->b);
    fp2_elem_init(e->c);
}

void fp6_elem_set(fp6_elem *e1, const fp6_elem *e2)
{
    fp2_elem_set(e1->a, e2->a);
    fp2_elem_set(e1->b, e2->b);
    fp2_elem_set(e1->c, e2->c);
}

void fp6_elem_from_str(fp6_elem *e,
                       const char *a1, const char *a2,
                       const char *b1, const char *b2,
                       const char *c1, const char *c2)
{
    e->a = calloc(1, sizeof(fp2_elem));
    e->b = calloc(1, sizeof(fp2_elem));
    e->c = calloc(1, sizeof(fp2_elem));

    fp2_elem_from_str(e->a, a1, a2);
    fp2_elem_from_str(e->b, b1, b2);
    fp2_elem_from_str(e->c, c1, c2);
}

void fp6_elem_free(fp6_elem *e)
{
    fp2_elem_free(e->a);
    fp2_elem_free(e->b);
    fp2_elem_free(e->c);
}

char *fp6_elem_get_str(const fp6_elem *e)
{
    char *str;

    str = calloc(1, sizeof(char));
    sprintf(str, "(%s) + (%s) * v + (%s) * v^2",
            fp2_elem_get_str(e->a),
            fp2_elem_get_str(e->b),
            fp2_elem_get_str(e->c));

    return str;
}

void fp6_add(fp6_elem *x, const fp6_elem *y, const fp6_elem *z)
{
    /* Algorithm 5.20 */
    fp2_add(x->a, y->a, z->a);
    fp2_add(x->b, y->b, z->b);
    fp2_add(x->c, y->c, z->c);
}

void fp6_sub(fp6_elem *x, const fp6_elem *y, const fp6_elem *z)
{
    fp2_sub(x->a, y->a, z->a);
    fp2_sub(x->b, y->b, z->b);
    fp2_sub(x->c, y->c, z->c);
}

void fp6_mul(fp6_elem *x, const fp6_elem *y, const fp6_elem *z)
{
    /* Algorithm 5.21 */
    fp2_elem v0, v1, v2, tmp,
             ra, rb, rc;

    fp2_elem_init(&v0);
    fp2_elem_init(&v1);
    fp2_elem_init(&v2);
    fp2_elem_init(&tmp);
    fp2_elem_init(&ra);
    fp2_elem_init(&rb);
    fp2_elem_init(&rc);

    fp2_mul(&v0, y->a, z->a);
    fp2_mul(&v1, y->b, z->b);
    fp2_mul(&v2, y->c, z->c);

    /* Line 4, assignment to c0 (x->a in our case). */
    fp2_add(&ra, y->b, y->c);
    fp2_add(&tmp, z->b, z->c);
    fp2_mul(&ra, &ra, &tmp);
    fp2_sub(&ra, &ra, &v1);
    fp2_sub(&ra, &ra, &v2);
    fp2_mul_nonresidue(&tmp, &ra);
    fp2_add(&ra, &tmp, &v0);

    /* Line 5, assignment to c1 (x->b in our case). */
    fp2_add(&rb, y->a, y->b);
    fp2_add(&tmp, z->a, z->b);
    fp2_mul(&rb, &rb, &tmp);
    fp2_sub(&rb, &rb, &v0);
    fp2_sub(&rb, &rb, &v1);
    fp2_mul_nonresidue(&tmp, &v2);
    fp2_add(&rb, &rb, &tmp);

    /* Line 6, assignement to c2 (x->c in our case). */
    fp2_add(&rc, y->a, y->c);
    fp2_add(&tmp, z->a, z->c);
    fp2_mul(&rc, &rc, &tmp);
    fp2_sub(&rc, &rc, &v0);
    fp2_sub(&rc, &rc, &v2);
    fp2_add(&rc, &rc, &v1);

    fp2_elem_set(x->a, &ra);
    fp2_elem_set(x->b, &rb);
    fp2_elem_set(x->c, &rc);
    fp2_elem_free(&ra);
    fp2_elem_free(&rb);
    fp2_elem_free(&rc);

    fp2_elem_free(&v0);
    fp2_elem_free(&v1);
    fp2_elem_free(&v2);
    fp2_elem_free(&tmp);
}

void fp6_square(fp6_elem *x, const fp6_elem *y)
{
    /* Algorithm 5.22 */
    fp2_elem v2, v3, v4, v5,
             ra, rb, rc;

    fp2_elem_init(&v2);
    fp2_elem_init(&v3);
    fp2_elem_init(&v4);
    fp2_elem_init(&v5);
    fp2_elem_init(&ra);
    fp2_elem_init(&rb);
    fp2_elem_init(&rc);

    fp2_mul(&v4, y->a, y->b);
    fp2_add(&v4, &v4, &v4);
    fp2_square(&v5, y->c);
    fp2_mul_nonresidue(&rb, &v5);
    fp2_add(&rb, &rb, &v4);

    fp2_sub(&v2, &v4, &v5);
    fp2_square(&v3, y->a);
    fp2_sub(&v4, y->a, y->b);
    fp2_add(&v4, &v4, y->c);
    fp2_square(&v4, &v4);
    fp2_mul(&v5, y->b, y->c);
    fp2_add(&v5, &v5, &v5);

    fp2_mul_nonresidue(&ra, &v5);
    fp2_add(&ra, &ra, &v3);

    fp2_add(&rc, &v2, &v4);
    fp2_add(&rc, &rc, &v5);
    fp2_sub(&rc, &rc, &v3);

    fp2_elem_set(x->a, &ra);
    fp2_elem_set(x->b, &rb);
    fp2_elem_set(x->c, &rc);
    fp2_elem_free(&ra);
    fp2_elem_free(&rb);
    fp2_elem_free(&rc);

    fp2_elem_free(&v2);
    fp2_elem_free(&v3);
    fp2_elem_free(&v4);
    fp2_elem_free(&v5);
}

void fp6_inv(fp6_elem *x, const fp6_elem *y)
{
    /* Algorithm 5.23 */
    fp2_elem v0, v1, v2, v3, v4, v5, v6,
             A, B, C, F,
             ra, rb, rc,
             tmp;

    fp2_elem_init(&v0);
    fp2_elem_init(&v1);
    fp2_elem_init(&v2);
    fp2_elem_init(&v3);
    fp2_elem_init(&v4);
    fp2_elem_init(&v5);
    fp2_elem_init(&v6);
    fp2_elem_init(&A);
    fp2_elem_init(&B);
    fp2_elem_init(&C);
    fp2_elem_init(&F);
    fp2_elem_init(&tmp);
    fp2_elem_init(&ra);
    fp2_elem_init(&rb);
    fp2_elem_init(&rc);

    fp2_square(&v0, y->a);
    fp2_square(&v1, y->b);
    fp2_square(&v2, y->c);

    fp2_mul(&v3, y->a, y->b);
    fp2_mul(&v4, y->a, y->c);
    fp2_mul(&v5, y->b, y->c);

    fp2_mul_nonresidue(&tmp, &v5);
    fp2_sub(&A, &v0, &tmp);

    fp2_mul_nonresidue(&tmp, &v2);
    fp2_sub(&B, &tmp, &v3);

    fp2_sub(&C, &v1, &v4);

    fp2_mul(&v6, y->a, &A);

    fp2_mul_nonresidue(&tmp, y->c);
    fp2_mul(&tmp, &tmp, &B);
    fp2_add(&v6, &v6, &tmp);
    fp2_mul_nonresidue(&tmp, y->b);
    fp2_mul(&tmp, &tmp, &C);
    fp2_add(&v6, &v6, &tmp);

    fp2_inv(&F, &v6);

    fp2_mul(&ra, &A, &F);
    fp2_mul(&rb, &B, &F);
    fp2_mul(&rc, &C, &F);

    fp2_elem_set(x->a, &ra);
    fp2_elem_set(x->b, &rb);
    fp2_elem_set(x->c, &rc);
    fp2_elem_free(&ra);
    fp2_elem_free(&rb);
    fp2_elem_free(&rc);

    fp2_elem_free(&v0);
    fp2_elem_free(&v1);
    fp2_elem_free(&v2);
    fp2_elem_free(&v3);
    fp2_elem_free(&v4);
    fp2_elem_free(&v5);
    fp2_elem_free(&v6);
    fp2_elem_free(&A);
    fp2_elem_free(&B);
    fp2_elem_free(&C);
    fp2_elem_free(&F);
    fp2_elem_free(&tmp);
}

void fp6_negate(fp6_elem *x, const fp6_elem *y)
{
    fp2_negate(x->a, y->a);
    fp2_negate(x->b, y->b);
    fp2_negate(x->c, y->c);
}

bool fp6_equal(const fp6_elem *e1, const fp6_elem *e2)
{
    return fp2_equal(e1->a, e2->a) &&
           fp2_equal(e1->b, e2->b) &&
           fp2_equal(e1->c, e2->c);
}

void fp12_elem_init(fp12_elem *e)
{
    e->a = calloc(1, sizeof(fp6_elem));
    e->b = calloc(1, sizeof(fp6_elem));

    fp6_elem_init(e->a);
    fp6_elem_init(e->b);
}

void fp12_elem_set(fp12_elem *e1, const fp12_elem *e2)
{
    fp6_elem_set(e1->a, e2->a);
    fp6_elem_set(e1->b, e2->b);
}

void fp12_elem_from_str(fp12_elem *e, const char *a[6], const char *b[6])
{
    e->a = calloc(1, sizeof(fp6_elem));
    e->b = calloc(1, sizeof(fp6_elem));

    fp6_elem_from_str(e->a, a[0], a[1], a[2], a[3], a[4], a[5]);
    fp6_elem_from_str(e->b, b[0], b[1], b[2], b[3], b[4], b[5]);
}

void fp12_elem_free(fp12_elem *e)
{
    fp6_elem_free(e->a);
    fp6_elem_free(e->b);
}

void fp12_add(fp12_elem *x, const fp12_elem *y, const fp12_elem *z)
{
    fp6_add(x->a, y->a, z->a);
    fp6_add(x->b, y->b, z->b);
}

void fp12_sub(fp12_elem *x, const fp12_elem *y, const fp12_elem *z)
{
    fp6_sub(x->a, y->a, z->a);
    fp6_sub(x->b, y->b, z->b);
}

static inline void fp6_mul_nonresidue(fp6_elem *x, const fp6_elem *y)
{
    /**
     * Need to multiply an element in F_p^6 by the quadratic non-resiude in this
     * field which is used to construct F_p^12.
     *
     * This amounts to (a + bv + bv^2)v = av + bv^2 + cv^3
     *                                  = av + bv^2 + c(u+1)  // since v^3 = u + 1
     *                                  = c(u+1) + av + bv^2  // re-arranged to F_p^6 format
     */
    fp2_elem tmp;

    fp2_elem_init(&tmp);

    fp2_mul_nonresidue(&tmp, y->c);

    fp2_elem_set(x->c, y->b);
    fp2_elem_set(x->b, y->a);
    fp2_elem_set(x->a, &tmp);

    fp2_elem_free(&tmp);
}

void fp12_mul(fp12_elem *x, const fp12_elem *y, const fp12_elem *z)
{
    /* Algorithm 5.16 */
    fp6_elem v0, v1, tmp,
             ra, rb;

    fp6_elem_init(&v0);
    fp6_elem_init(&v1);
    fp6_elem_init(&tmp);
    fp6_elem_init(&ra);
    fp6_elem_init(&rb);

    fp6_mul(&v0, y->a, z->a);
    fp6_mul(&v1, y->b, z->b);

    fp6_mul_nonresidue(&tmp, &v1);
    fp6_add(&ra, &v0, &tmp);

    fp6_add(&rb, y->a, y->b);
    fp6_add(&tmp, z->a, z->b);
    fp6_mul(&rb, &rb, &tmp);
    fp6_sub(&rb, &rb, &v0);
    fp6_sub(&rb, &rb, &v1);

    fp6_elem_set(x->a, &ra);
    fp6_elem_set(x->b, &rb);

    fp6_elem_free(&v0);
    fp6_elem_free(&v1);
    fp6_elem_free(&tmp);
    fp6_elem_free(&ra);
    fp6_elem_free(&rb);
}

void fp12_square(fp12_elem *x, const fp12_elem *y)
{
    /* Algorithm 5.17 */
    fp6_elem v0, v1, v2, ra, rb;

    fp6_elem_init(&v0);
    fp6_elem_init(&v1);
    fp6_elem_init(&v2);
    fp6_elem_init(&ra);
    fp6_elem_init(&rb);

    fp6_sub(&v0, y->a, y->b);

    fp6_mul_nonresidue(&v2, y->b);
    fp6_sub(&v2, y->a, &v2);

    fp6_mul(&v1, y->a, y->b);
    fp6_mul(&v0, &v0, &v2);
    fp6_add(&v0, &v0, &v1);

    fp6_add(&rb, &v1, &v1);

    fp6_mul_nonresidue(&v1, &v1);
    fp6_add(&ra, &v0, &v1);

    fp6_elem_set(x->a, &ra);
    fp6_elem_set(x->b, &rb);

    fp6_elem_free(&v0);
    fp6_elem_free(&v1);
    fp6_elem_free(&v2);
    fp6_elem_free(&ra);
    fp6_elem_free(&rb);
}

void fp12_inv(fp12_elem *x, const fp12_elem *y)
{
    /* Algorithm 5.19 */
    fp6_elem v0, v1, ra, rb, tmp;

    fp6_elem_init(&v0);
    fp6_elem_init(&v1);
    fp6_elem_init(&ra);
    fp6_elem_init(&rb);
    fp6_elem_init(&tmp);

    fp6_square(&v0, y->a);
    fp6_square(&v1, y->b);

    fp6_mul_nonresidue(&v1, &v1);
    fp6_sub(&v0, &v0, &v1);
    fp6_inv(&v1, &v0);

    fp6_mul(&ra, y->a, &v1);

    // TODO(nr): I should probably have fp<k>_negate functions.
    fp6_sub(&tmp, &tmp, y->b);
    fp6_mul(&rb, &tmp, &v1);

    fp6_elem_set(x->a, &ra);
    fp6_elem_set(x->b, &rb);

    fp6_elem_free(&v0);
    fp6_elem_free(&v1);
    fp6_elem_free(&ra);
    fp6_elem_free(&rb);
    fp6_elem_free(&tmp);
}

void fp12_pow(fp12_elem *x, const fp12_elem *y, const mpz_t exp)
{
    fp12_elem r;
    mpz_t e;
    mp_bitcnt_t L;
    int sign;

    fp12_elem_init(&r);
    mpz_init_set(e, exp);

    fp12_elem_set(&r, y);

    /**
     * We always need a positive integer representation
     * to use the GMP bit fiddling functions. If the exp
     * is negative, that just means the binary representation
     * looks like:
     *     exp = -2^n0 - 2^n2 ... -2^nL-1
     * which means that at the end of the loop, we just take
     * the inverse of the result.
     */
    if ((sign = mpz_sgn(e)) < 0) {
        mpz_neg(e, e);
    }

    L = (mp_bitcnt_t)mpz_sizeinbase(e, 2) - 2;
    for (;;) {
        fp12_square(&r, &r);

        if (mpz_tstbit(e, L)) {
            fp12_mul(&r, &r, y);
        }

        if (!L) break;
        L--;
    }

    if (sign < 0) {
        fp12_inv(&r, &r);
    }

    fp12_elem_set(x, &r);

    fp12_elem_free(&r);
    mpz_clear(e);
}

void fp12_conjugate(fp12_elem *x, const fp12_elem *y)
{
    fp6_elem_set(x->a, y->a);
    fp6_negate(x->b, y->b);
}

bool fp12_equal(const fp12_elem *e1, const fp12_elem *e2)
{
    return fp6_equal(e1->a, e2->a) && fp6_equal(e1->b, e2->b);
}

/**
 * From https://eprint.iacr.org/2010/354.pdf
 *
 * We will use algorithm 28 to compute the frobenius in F_p^12. This involves
 * using some coefficients, denoted as "gammas" in the paper, which are constant
 * in terms of the field characteristic. So, we can pre-compute the values and
 * store them here.
 *
 * The coefficients are elements of F_p^2, so are of the form a + bu. For convenience,
 * use parallel arrays to store the a's and b's.
 */
const char* const g_fp12_frob_coeff_a[5] = {
    "0x1904d3bf02bb0667c231beb4202c0d1f0fd603fd3cbd5f4f7b2443d784bab9c4f67ea53d63e7813d8d0775ed92235fb8",
    "0x0",
    "0x06af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09",
    "0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad",
    "0x05b2cfd9013a5fd8df47fa6b48b1e045f39816240c0b8fee8beadf4d8e9c0566c63a3e6e257f87329b18fae980078116"
};

const char* const g_fp12_frob_coeff_b[5] = {
    "0xfc3e2b36c4e03288e9e902231f9fb854a14787b6c7b36fec0c8ec971f63c5f282d5ac14d6c7ec22cf78a126ddc4af3",
    "0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac",
    "0x06af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09",
    "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
    "0x144e4211384586c16bd3ad4afa99cc9170df3560e77982d0db45f3536814f0bd5871c1908bd478cd1ee605167ff82995"
};

void fp12_frobenius(fp12_elem *x, const fp12_elem *y)
{
    fp2_elem gamma1, gamma2, gamma3, gamma4, gamma5,
             t1, t2, t3, t4, t5, t6;

    fp2_elem_from_str(&gamma1, g_fp12_frob_coeff_a[0], g_fp12_frob_coeff_b[0]);
    fp2_elem_from_str(&gamma2, g_fp12_frob_coeff_a[1], g_fp12_frob_coeff_b[1]);
    fp2_elem_from_str(&gamma3, g_fp12_frob_coeff_a[2], g_fp12_frob_coeff_b[2]);
    fp2_elem_from_str(&gamma4, g_fp12_frob_coeff_a[3], g_fp12_frob_coeff_b[3]);
    fp2_elem_from_str(&gamma5, g_fp12_frob_coeff_a[4], g_fp12_frob_coeff_b[4]);

    fp2_elem_init(&t1);
    fp2_elem_init(&t2);
    fp2_elem_init(&t3);
    fp2_elem_init(&t4);
    fp2_elem_init(&t5);
    fp2_elem_init(&t6);

    fp2_conjugate(&t1, y->a->a);
    fp2_conjugate(&t2, y->b->a);
    fp2_conjugate(&t3, y->a->b);
    fp2_conjugate(&t4, y->b->b);
    fp2_conjugate(&t5, y->a->c);
    fp2_conjugate(&t6, y->b->c);

    fp2_mul(&t2, &t2, &gamma1);
    fp2_mul(&t3, &t3, &gamma2);
    fp2_mul(&t4, &t4, &gamma3);
    fp2_mul(&t5, &t5, &gamma4);
    fp2_mul(&t6, &t6, &gamma5);

    fp2_elem_set(x->a->a, &t1);
    fp2_elem_set(x->b->a, &t2);
    fp2_elem_set(x->a->b, &t3);
    fp2_elem_set(x->b->b, &t4);
    fp2_elem_set(x->a->c, &t5);
    fp2_elem_set(x->b->c, &t6);

    fp2_elem_free(&t1);
    fp2_elem_free(&t2);
    fp2_elem_free(&t3);
    fp2_elem_free(&t4);
    fp2_elem_free(&t5);
    fp2_elem_free(&t6);
    fp2_elem_free(&gamma1);
    fp2_elem_free(&gamma2);
    fp2_elem_free(&gamma3);
    fp2_elem_free(&gamma4);
    fp2_elem_free(&gamma5);
}
