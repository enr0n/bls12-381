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
#include "params.h"

/**
 * The prime P which is the characteristic of the finite fields.
 *
 * This variable must be initialized with fp_params_init(), and
 * free'd iwth fp_params_free(). Each of these functions must
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

void fp2_elem_from_str(fp2_elem *e, const char *a, const char *b)
{
    mpz_init_set_str(e->a, a, 0);
    mpz_init_set_str(e->b, b, 0);
}

void fp2_elem_clear(fp2_elem *e)
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

int fp2_equal(const fp2_elem *e1, const fp2_elem *e2)
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

void fp6_elem_clear(fp6_elem *e)
{
    fp2_elem_clear(e->a);
    fp2_elem_clear(e->b);
    fp2_elem_clear(e->c);
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

static inline void fp2_mul_by_fp6_alpha(fp2_elem *x, const fp2_elem *y)
{
    /**
     * Need to multiply an element in F_p^2 by the "alpha" parameter of F_p^6,
     * which in this case is u + 1.
     *
     * This amounts to (a + bu)(u + 1) = au + a + bu^2 + bu
     *                                 = au + a - b + bu      // since u^2 = beta = -1
     *                                 = (a - b) + (a + bu)
     */
    assert(x != y);
    fp_sub(x->a, y->a, y->b);
    fp_add(x->b, y->a, y->b);
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
    fp2_mul_by_fp6_alpha(&tmp, &ra);
    fp2_add(&ra, &tmp, &v0);

    /* Line 5, assignment to c1 (x->b in our case). */
    fp2_add(&rb, y->a, y->b);
    fp2_add(&tmp, z->a, z->b);
    fp2_mul(&rb, &rb, &tmp);
    fp2_sub(&rb, &rb, &v0);
    fp2_sub(&rb, &rb, &v1);
    fp2_mul_by_fp6_alpha(&tmp, &v2);
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
    fp2_elem_clear(&ra);
    fp2_elem_clear(&rb);
    fp2_elem_clear(&rc);

    fp2_elem_clear(&v0);
    fp2_elem_clear(&v1);
    fp2_elem_clear(&v2);
    fp2_elem_clear(&tmp);
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
    fp2_mul_by_fp6_alpha(&rb, &v5);
    fp2_add(&rb, &rb, &v4);

    fp2_sub(&v2, &v4, &v5);
    fp2_square(&v3, y->a);
    fp2_sub(&v4, y->a, y->b);
    fp2_add(&v4, &v4, y->c);
    fp2_square(&v4, &v4);
    fp2_mul(&v5, y->b, y->c);
    fp2_add(&v5, &v5, &v5);

    fp2_mul_by_fp6_alpha(&ra, &v5);
    fp2_add(&ra, &ra, &v3);

    fp2_add(&rc, &v2, &v4);
    fp2_add(&rc, &rc, &v5);
    fp2_sub(&rc, &rc, &v3);

    fp2_elem_set(x->a, &ra);
    fp2_elem_set(x->b, &rb);
    fp2_elem_set(x->c, &rc);
    fp2_elem_clear(&ra);
    fp2_elem_clear(&rb);
    fp2_elem_clear(&rc);

    fp2_elem_clear(&v2);
    fp2_elem_clear(&v3);
    fp2_elem_clear(&v4);
    fp2_elem_clear(&v5);
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

    fp2_mul_by_fp6_alpha(&tmp, &v5);
    fp2_sub(&A, &v0, &tmp);

    fp2_mul_by_fp6_alpha(&tmp, &v2);
    fp2_sub(&B, &tmp, &v3);

    fp2_sub(&C, &v1, &v4);

    fp2_mul(&v6, y->a, &A);

    fp2_mul_by_fp6_alpha(&tmp, y->c);
    fp2_mul(&tmp, &tmp, &B);
    fp2_add(&v6, &v6, &tmp);
    fp2_mul_by_fp6_alpha(&tmp, y->b);
    fp2_mul(&tmp, &tmp, &C);
    fp2_add(&v6, &v6, &tmp);

    fp2_inv(&F, &v6);

    fp2_mul(&ra, &A, &F);
    fp2_mul(&rb, &B, &F);
    fp2_mul(&rc, &C, &F);

    fp2_elem_set(x->a, &ra);
    fp2_elem_set(x->b, &rb);
    fp2_elem_set(x->c, &rc);
    fp2_elem_clear(&ra);
    fp2_elem_clear(&rb);
    fp2_elem_clear(&rc);

    fp2_elem_clear(&v0);
    fp2_elem_clear(&v1);
    fp2_elem_clear(&v2);
    fp2_elem_clear(&v3);
    fp2_elem_clear(&v4);
    fp2_elem_clear(&v5);
    fp2_elem_clear(&v6);
    fp2_elem_clear(&A);
    fp2_elem_clear(&B);
    fp2_elem_clear(&C);
    fp2_elem_clear(&F);
    fp2_elem_clear(&tmp);
}

int fp6_equal(const fp6_elem *e1, const fp6_elem *e2)
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

void fp12_elem_from_str(fp12_elem *e, const char *a[6], const char *b[6])
{
    e->a = calloc(1, sizeof(fp6_elem));
    e->b = calloc(1, sizeof(fp6_elem));

    fp6_elem_from_str(e->a, a[0], a[1], a[2], a[3], a[4], a[5]);
    fp6_elem_from_str(e->b, b[0], b[1], b[2], b[3], b[4], b[5]);
}

void fp12_elem_clear(fp12_elem *e)
{
    fp6_elem_clear(e->a);
    fp6_elem_clear(e->b);
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

static inline void fp6_mul_by_fp12_alpha(fp6_elem *x, const fp6_elem *y)
{
    /**
     * Need to multiply an element in F_p^6 by the "alpha" parameter of F_p^12,
     * which in this case is v.
     *
     * This amounts to (a + bv + bv^2)v = av + bv^2 + cv^3
     *                                  = av + bv^2 + c(u+1)  // since v^3 = u + 1
     *                                  = c(u+1) + av + bv^2  // re-arranged to F_p^6 format
     */
    fp2_elem tmp;

    fp2_elem_init(&tmp);

    fp2_mul_by_fp6_alpha(&tmp, y->c);

    fp2_elem_set(x->c, y->b);
    fp2_elem_set(x->b, y->a);
    fp2_elem_set(x->a, &tmp);

    fp2_elem_clear(&tmp);
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

    fp6_mul_by_fp12_alpha(&tmp, &v1);
    fp6_add(&ra, &v0, &tmp);

    fp6_add(&rb, y->a, y->b);
    fp6_add(&tmp, z->a, z->b);
    fp6_mul(&rb, &rb, &tmp);
    fp6_sub(&rb, &rb, &v0);
    fp6_sub(&rb, &rb, &v1);

    fp6_elem_set(x->a, &ra);
    fp6_elem_set(x->b, &rb);

    fp6_elem_clear(&v0);
    fp6_elem_clear(&v1);
    fp6_elem_clear(&tmp);
    fp6_elem_clear(&ra);
    fp6_elem_clear(&rb);
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

    fp6_mul_by_fp12_alpha(&v2, y->b);
    fp6_sub(&v2, y->a, &v2);

    fp6_mul(&v1, y->a, y->b);
    fp6_mul(&v0, &v0, &v2);
    fp6_add(&v0, &v0, &v1);

    fp6_add(&rb, &v1, &v1);

    fp6_mul_by_fp12_alpha(&v1, &v1);
    fp6_add(&ra, &v0, &v1);

    fp6_elem_set(x->a, &ra);
    fp6_elem_set(x->b, &rb);

    fp6_elem_clear(&v0);
    fp6_elem_clear(&v1);
    fp6_elem_clear(&v2);
    fp6_elem_clear(&ra);
    fp6_elem_clear(&rb);
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

    fp6_mul_by_fp12_alpha(&v1, &v1);
    fp6_sub(&v0, &v0, &v1);
    fp6_inv(&v1, &v0);

    fp6_mul(&ra, y->a, &v1);

    // TODO(nr): I should probably have fp<k>_negate functions.
    fp6_sub(&tmp, &tmp, y->b);
    fp6_mul(&rb, &tmp, &v1);

    fp6_elem_set(x->a, &ra);
    fp6_elem_set(x->b, &rb);

    fp6_elem_clear(&v0);
    fp6_elem_clear(&v1);
    fp6_elem_clear(&ra);
    fp6_elem_clear(&rb);
    fp6_elem_clear(&tmp);
}

int fp12_equal(const fp12_elem *e1, const fp12_elem *e2)
{
    return fp6_equal(e1->a, e2->a) && fp6_equal(e1->b, e2->b);
}
