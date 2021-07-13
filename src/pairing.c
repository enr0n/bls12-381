#include <assert.h>
#include <stdio.h>
#include <gmp.h>

#include "params.h"
#include "finite_field.h"

#include "pairing.h"

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

    fp12_elem_clear(&fr);
    fp12_elem_clear(&t0);
    fp12_elem_clear(&t1);
    fp12_elem_clear(&t2);
    fp12_elem_clear(&t3);
}
