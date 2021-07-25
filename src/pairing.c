#include <assert.h>
#include <stdio.h>
#include <gmp.h>

#include "finite_field.h"
#include "pairing.h"

#include "BLS12_381.h"

/**
 * The curve parameter t which parameterizes p(t) and r(t).
 *
 * This variable is initailized once with BLS12_381_init(), and free'd
 * with BLS12_381_free(). Each of these must only be called once, and
 * only by external callers.
 *
 * This variable is never modified.
 */
static mpz_t g_BLS12_381_t;

void BLS12_381_init()
{
    fp_params_init();
    mpz_init_set_str(g_BLS12_381_t, BLS12_381_t, 0);
}

void BLS12_381_free()
{
    fp_params_free();
    mpz_clear(g_BLS12_381_t);
}

void final_exponentiation(fp12_elem *r, const fp12_elem *f)
{
    /**
     * Follow the theorem outlined in https://eprint.iacr.org/2020/875.pdf
     * for the hard part of the exponentiation.
     */
    fp12_elem fr, t0, t1, t2, t3;

    fp12_elem_init(&fr);
    fp12_elem_init(&t0);
    fp12_elem_init(&t1);
    fp12_elem_init(&t2);
    fp12_elem_init(&t3);

    // f = f^(p^6 - 1)
    fp12_conjugate(&t0, f);
    fp12_inv(&t1, f);
    fp12_mul(&fr, &t0, &t1);

    // f = f * f^(p^2 + 1)
    //   = f^(p^6 - 1)(p^2 + 1)
    //
    // ...which concludes the "easy" part of final exponetiation...
    fp12_frobenius(&t0, &fr);
    fp12_frobenius(&t0, &t0);
    fp12_mul(&fr, &t0, &fr);

    // And now we do the "hard" part...
    fp12_square(&t0, &fr);
    fp12_mul(&t0, &t0, &fr);

    // (x-1)^2
    fp12_pow(&t1, &fr, g_BLS12_381_t);
    fp12_inv(&t2, &fr);
    fp12_mul(&fr, &t1, &t2);
    fp12_pow(&t1, &fr, g_BLS12_381_t);
    fp12_inv(&t2, &fr);
    fp12_mul(&fr, &t1, &t2);

    // (x + p)
    fp12_pow(&t1, &fr, g_BLS12_381_t);
    fp12_frobenius(&t2, &fr);
    fp12_mul(&fr, &t1, &t2);

    // (x^2 + p^2 - 1)
    fp12_pow(&t1, &fr, g_BLS12_381_t);
    fp12_pow(&t1, &t1, g_BLS12_381_t);
    fp12_frobenius(&t2, &fr);
    fp12_frobenius(&t2, &t2);
    fp12_inv(&t3, &fr);
    fp12_mul(&fr, &t1, &t2);
    fp12_mul(&fr, &fr, &t3);

    fp12_mul(&fr, &fr, &t0);

    fp12_elem_set(r, &fr);

    fp12_elem_free(&fr);
    fp12_elem_free(&t0);
    fp12_elem_free(&t1);
    fp12_elem_free(&t2);
    fp12_elem_free(&t3);
}

static void line_function_double(fp12_elem *l, const G2_elem_proj *A, const G1_elem_affine *P)
{
    // The N corresponds to the degree of the term.
    fp2_elem r0, r2, r3, tmp;
    fp12_elem r;

    fp2_elem_init(&r0);
    fp2_elem_init(&r2);
    fp2_elem_init(&r3);
    fp2_elem_init(&tmp);
    fp12_elem_init(&r);

    // r3 = -2YZ*y_P
    fp2_mul(&r3, A->y, A->z);
    fp2_negate(&r3, &r3);
    fp2_add(&r3, &r3, &r3);
    fp2_mul_scalar(&r3, &r3, P->y);

    // r2 = 3X^2*x_P
    fp2_square(&tmp, A->x);
    fp2_add(&r2, &tmp, &tmp);
    fp2_add(&r2, &r2, &tmp);
    fp2_mul_scalar(&r2, &r2, P->x);

    // r0 = 3bZ^2(u+1) - Y^2
    fp2_square(&r0, A->z);
    fp2_elem_set_si(&tmp, 12, 0);
    fp2_mul_nonresidue(&tmp, &tmp);
    fp2_mul(&r0, &r0, &tmp);
    fp2_square(&tmp, A->y);
    fp2_sub(&r0, &r0, &tmp);

    fp2_elem_set(r.a->a, &r0);
    fp2_elem_set(r.a->b, &r2);
    fp2_elem_set(r.b->b, &r3);

    fp12_elem_set(l, &r);

    fp2_elem_free(&r0);
    fp2_elem_free(&r2);
    fp2_elem_free(&r3);
    fp2_elem_free(&tmp);
    fp12_elem_free(&r);
}

static void line_function_add(fp12_elem *l, const G2_elem_proj *A, const G2_elem_proj *B, const G1_elem_affine *P)
{
    // The N corresponds to the degree of the term.
    fp2_elem r0, r2, r3,
             t1, t2, tmp;
    fp12_elem r;

    fp2_elem_init(&r0);
    fp2_elem_init(&r2);
    fp2_elem_init(&r3);
    fp2_elem_init(&t1);
    fp2_elem_init(&t2);
    fp2_elem_init(&tmp);
    fp12_elem_init(&r);

    // t1 = (X1 - X2*Z1)
    fp2_mul(&t1, B->x, A->z);
    fp2_sub(&t1, A->x, &t1);

    // t2 = (Y1 - Y2*Z1)
    fp2_mul(&t2, B->y, A->z);
    fp2_sub(&t2, A->y, &t2);

    // r3 = (X1 - X2*Z1)*y_p
    fp2_mul_scalar(&r3, &t1, P->y);

    // r2 = -(X1 - Y2*Z1)*x_p
    fp2_mul_scalar(&r2, &t2, P->x);
    fp2_negate(&r2, &r2);

    // r0 = (Y1 - Y2*Z1)*X2 - (X1 - X2*Z1)*Y2
    fp2_mul(&r0, &t2, B->x);
    fp2_mul(&tmp, &t1, B->y);
    fp2_sub(&r0, &r0, &tmp);

    fp2_elem_set(r.a->a, &r0);
    fp2_elem_set(r.a->b, &r2);
    fp2_elem_set(r.b->b, &r3);

    fp12_elem_set(l, &r);

    fp2_elem_free(&r0);
    fp2_elem_free(&r2);
    fp2_elem_free(&r3);
    fp2_elem_free(&t1);
    fp2_elem_free(&t2);
    fp2_elem_free(&tmp);
    fp12_elem_free(&r);
}

/**
 * Follow the line equations for projective cooridnates laid out
 * here:
 *     https://ieeexplore.ieee.org/document/7350136
 *
 * N.B. that because we are dealing with an M-type sextic twist,
 * we twist the point P from G_1 onto the twist, rather than untwisting
 * the point in G_2. As noted in the paper, this is simply a less expensive
 * mapping:
 *
 *     (x, y) -> (x*w^2, y*w^3), where w^6 = (u + 1)
 *
 * So, this amounts to assigning particular coefficients to contruct our
 * sparse element in F_p^12.
 */
void line_function(fp12_elem *l, G2_elem_proj *T,
                   const G2_elem_proj *A, const G2_elem_proj *B,
                   const G1_elem_affine *P)
{
    // This is intentionally a pointer comparison.
    if (A == B) {
        // T = 2A
        line_function_double(l, A, P);
        G2_double_proj(T, A);
    } else {
        // T = A + B
        line_function_add(l, A, B, P);
        G2_add_proj(T, A, B);
    }
}

/**
 * During unit tests, we want to inspect the final value of T; the projective
 * point accumalator. This should help failing tests distinguish beween line function
 * errors, and miller loop errors.
 */
void miller_loop(fp12_elem *r, G2_elem_affine *R, const G1_elem_affine *P, const G2_elem_affine *Q)
{
    fp12_elem fr, l, tmp;
    G2_elem_proj T, Q_proj, Q_neg;
    mp_bitcnt_t c, nbits;
    int sign;
    mpz_t e;

    fp12_elem_init(&fr);
    fp12_elem_init(&l);
    fp12_elem_init(&tmp);

    G2_identity_init_proj(&T);
    G2_identity_init_proj(&Q_proj);
    G2_identity_init_proj(&Q_neg);

    mpz_init_set(e, g_BLS12_381_t);

    if ((sign = mpz_sgn(e)) < 0) {
        mpz_neg(e, e);
    }

    G2_affine2proj(&T, Q);
    G2_affine2proj(&Q_proj, Q);
    G2_affine2proj(&Q_neg, Q); G2_negate_proj(&Q_neg, &Q_neg);

    mpz_set_si(fr.a->a->a, 1);

    c = nbits = (mp_bitcnt_t)mpz_sizeinbase(e, 2) - 1;
    for (;;) {
        fp12_square(&fr, &fr);
        line_function(&l, &T, &T, &T, P);

        /* No addition on the first iteration! */
        if (mpz_tstbit(e, c) && c != nbits) {
            /**
             * If the curve parameter is negative, do a subtraction at the first
             * addition step.
             */
            if (sign < 0 && c == nbits - 1) {
                line_function(&tmp, &T, &T, &Q_neg, P);
            } else {
                line_function(&tmp, &T, &T, &Q_proj, P);
            }

            fp12_mul(&l, &l, &tmp);
        }

        fp12_mul(&fr, &fr, &l);

        if (!c) break;
        c--;
    }

    if (sign < 0) {
        fp12_conjugate(&fr, &fr);
    }

    fp12_elem_set(r, &fr);
    G2_proj2affine(R, &T);

    fp12_elem_free(&fr);
    fp12_elem_free(&l);
    fp12_elem_free(&tmp);
    G2_elem_free_proj(&T);
    G2_elem_free_proj(&Q_proj);
    G2_elem_free_proj(&Q_neg);
    mpz_clear(e);
}

void ate(fp12_elem *r, const G2_elem_affine *Q, const G1_elem_affine *P)
{
    G2_elem_affine R;
    G2_identity_init_affine(&R);

    miller_loop(r, &R, P, Q);
    final_exponentiation(r, r);

    G2_elem_free_affine(&R);
}
