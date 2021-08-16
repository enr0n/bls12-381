/**
 * These are functions for finite field arithmetic that are likely
 * only needed internally. External users should not normally need
 * to perform arithmetic outside of G_1, G_2, and F_p^12.
 */
#ifndef _FINITE_FIELD_H_
#define _FINITE_FIELD_H_

#include <stdbool.h>
#include <gmp.h>

#include "BLS12_381.h"

void fp_params_init();
void fp_params_free();

/* Arithmetic operations in the base field F_p. */
bool fp_equiv(const mpz_t x, const mpz_t y);
void fp_add(mpz_t x, const mpz_t y, const mpz_t z);
void fp_sub(mpz_t x, const mpz_t y, const mpz_t z);
void fp_mul(mpz_t x, const mpz_t y, const mpz_t z);
void fp_inv(mpz_t x, const mpz_t y);
void fp_negate(mpz_t x, const mpz_t y);
void fp_mul_ui(mpz_t, const mpz_t y, unsigned long int z);
void fp_add_ui(mpz_t, const mpz_t y, unsigned long int z);

void fp2_elem_init(fp2_elem *e);
void fp2_elem_set(fp2_elem *e1, const fp2_elem *e2);
void fp2_elem_set_si(fp2_elem *e, signed long int a, signed long int b);
void fp2_elem_from_str(fp2_elem *e, const char *a, const char *b);
void fp2_elem_set_str(fp2_elem *e, const char *a, const char *b);
void fp2_elem_free(fp2_elem *e);

char *fp2_elem_get_str(const fp2_elem *e);

/* Arithmetic operations in the field extension F_p^2. */
void fp2_add(fp2_elem *x, const fp2_elem *y, const fp2_elem *z);
void fp2_sub(fp2_elem *x, const fp2_elem *y, const fp2_elem *z);
void fp2_mul(fp2_elem *x, const fp2_elem *y, const fp2_elem *z);
void fp2_square(fp2_elem *x, const fp2_elem *y);
void fp2_inv(fp2_elem *x, const fp2_elem *y);
void fp2_pow(fp2_elem *x, const fp2_elem *y, const mpz_t exp);
void fp2_negate(fp2_elem *x, const fp2_elem *y);
void fp2_conjugate(fp2_elem *x, const fp2_elem *y);
void fp2_mul_nonresidue(fp2_elem *x, const fp2_elem *y);
void fp2_mul_scalar(fp2_elem *x, const fp2_elem *y, const mpz_t m);
bool fp2_equal(const fp2_elem *e1, const fp2_elem *e2);

void fp6_elem_init(fp6_elem *e);
void fp6_elem_set(fp6_elem *e1, const fp6_elem *e2);
void fp6_elem_from_str(fp6_elem *e,
                       const char *a1, const char *a2,
                       const char *b1, const char *b2,
                       const char *c1, const char *c2);
void fp6_elem_free(fp6_elem *e);

char *fp6_elem_get_str(const fp6_elem *e);

/**
 * Arithmetic operations in the field extension F_p^6. This field is implemented
 * as a cubic extension of F_p^2.
 */
void fp6_add(fp6_elem *x, const fp6_elem *y, const fp6_elem *z);
void fp6_sub(fp6_elem *x, const fp6_elem *y, const fp6_elem *z);
void fp6_mul(fp6_elem *x, const fp6_elem *y, const fp6_elem *z);
void fp6_square(fp6_elem *x, const fp6_elem *y);
void fp6_inv(fp6_elem *x, const fp6_elem *y);
void fp6_negate(fp6_elem *x, const fp6_elem *y);
bool fp6_equal(const fp6_elem *e1, const fp6_elem *e2);

#endif // _FINITE_FIELD_H_
