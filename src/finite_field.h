#ifndef _FINITE_FIELD_H_
#define _FINITE_FIELD_H_

#include <gmp.h>

void fp_params_init();
void fp_params_free();

/* Arithmetic operations in the base field F_p. */
void fp_add(mpz_t x, const mpz_t y, const mpz_t z);
void fp_sub(mpz_t x, const mpz_t y, const mpz_t z);
void fp_mul(mpz_t x, const mpz_t y, const mpz_t z);
void fp_inv(mpz_t x, const mpz_t y);

/* An element of the field F_p^2, a + bu, where a,b belong to F_p. */
typedef struct fp2_elem {
    mpz_t a;
    mpz_t b;
} fp2_elem;

void fp2_elem_init(fp2_elem *e);
void fp2_elem_set(fp2_elem *e1, const fp2_elem *e2);
void fp2_elem_from_str(fp2_elem *e, const char *a, const char *b);
void fp2_elem_clear(fp2_elem *e);

char *fp2_elem_get_str(const fp2_elem *e);

/* Arithmetic operations in the field extension F_p^2. */
void fp2_add(fp2_elem *x, const fp2_elem *y, const fp2_elem *z);
void fp2_sub(fp2_elem *x, const fp2_elem *y, const fp2_elem *z);
void fp2_mul(fp2_elem *x, const fp2_elem *y, const fp2_elem *z);
void fp2_square(fp2_elem *x, const fp2_elem *y);
void fp2_inv(fp2_elem *x, const fp2_elem *y);

int fp2_equal(const fp2_elem *e1, const fp2_elem *e2);

/* An element of the field F_p^6. */
typedef struct fp6_elem {
    fp2_elem *a;
    fp2_elem *b;
    fp2_elem *c;
} fp6_elem;

void fp6_elem_init(fp6_elem *e);
void fp6_elem_set(fp6_elem *e1, const fp6_elem *e2);
void fp6_elem_from_str(fp6_elem *e,
                       const char *a1, const char *a2,
                       const char *b1, const char *b2,
                       const char *c1, const char *c2);
void fp6_elem_clear(fp6_elem *e);

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

int fp6_equal(const fp6_elem *e1, const fp6_elem *e2);

/* An element of the field F_p^12. */
typedef struct fp12_elem {
    fp6_elem *a;
    fp6_elem *b;
} fp12_elem;

void fp12_elem_init(fp12_elem *e);
void fp12_elem_set(fp12_elem *e1, const fp12_elem *e2);
void fp12_elem_from_str(fp12_elem *e, const char *a[6], const char *b[6]);
void fp12_elem_clear(fp12_elem *e);

/**
 * Arithmetic operations in the field extension F_p^12. This field is implemented
 * as a quadratic extension of F_p^6.
 */
void fp12_add(fp12_elem *x, const fp12_elem *y, const fp12_elem *z);
void fp12_sub(fp12_elem *x, const fp12_elem *y, const fp12_elem *z);
void fp12_mul(fp12_elem *x, const fp12_elem *y, const fp12_elem *z);
void fp12_square(fp12_elem *x, const fp12_elem *y);
void fp12_inv(fp12_elem *x, const fp12_elem *y);

int fp12_equal(const fp12_elem *e1, const fp12_elem *e2);

#endif // _FINITE_FIELD_H_
