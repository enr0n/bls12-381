#ifndef BLS12_381_H
#define BLS12_381_H

#include <stdbool.h>
#include <gmp.h>

/**
 * BLS12-381 is a pairing-friendly curve from the BLS (Barreto, Lynn, Scott)
 * family of curves with embedding degree k=12. BLS12 curves are parameterized by
 * the polynomials p(t) = (t - 1)^2 * (t^4 - t^2 + 1) / 3 + t, and
 * r(t) = t^4 - t^2 + 1. The curve equation is E: y^2 = x^3 + 4, which admits a
 * a sextic twist E': y^2 = x^3 + 4(u + 1) where u is neither a quadratic- nor
 * cubic- resiude mod p, and p is the prime characteristic of finite field F_p.
 *
 * Our pairing is defined e: G_2 x G_1 -> G_T, where G_1, G_2 and G_T each have
 * prime order r.  G_1 is defined over the finite field F_p, G_2 is defined
 * over the quadratic extension F_p^2, and G_T is a subgroup of the full
 * extension field F_p^12. The fields are constructed using the following
 * tower:
 *
 *      GF(p^2)  = GF(p)[u] / (u^2 + 1)
 *      GF(p^6)  = GF(p^2)[v] / (v^3 - u - 1)
 *      GF(p^12) = GF(p^6)[w] / (w^2 - v)
 *
 * Thus, we have G_1 is a subgroup of E(F_p), G_2 is a subgroup of E'(F_p^2),
 * and G_T is the r roots of unity from the multiplicative group of F_p^12.
 *
 * BLS12-381's sextic twist is M-type. This means that the twisting isomorphism
 * from E -> E' is more efficient than the untwisting isomorphism from E' -> E.
 * For this reason, rather than untwisting points on G_2 during the pairing
 * computation, we twist points in G_1, thus computing the whole pairing on the
 * curve. This is called the twisted ate pairing.
 *
 * Below are the parmaters of BLS12-381.
 */

/**
 * The parameter such that p(t), r(t) give the prime field characteristic and group orders.
 */
#define BLS12_381_t "-0xd201000000010000"

/**
 * The prime field characteristic.
 */
#define BLS12_381_P "0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab"

/**
 * The x coordinate for a generator of G_1.
 */
#define BLS12_381_G1_GEN_X "0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb"

/**
 * The y coordinate for a generator of G_1.
 */
#define BLS12_381_G1_GEN_Y "0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1"

/**
 * The first term of the x coordinate for a generator of G_2.
 */
#define BLS12_381_G2_GEN_X0 "0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8"

/**
 * The second term of the x coordinate for a generator of G_2.
 */
#define BLS12_381_G2_GEN_X1 "0x13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e"

/**
 * The first term of the y coordinate for a generator of G_2.
 */
#define BLS12_381_G2_GEN_Y0 "0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801"

/**
 * The second term of the y coordinate for a generator of G_2.
 */
#define BLS12_381_G2_GEN_Y1 "0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be"

/**
 * The curve parameters must be initialized before calling any other BLS12_381
 * functions. It is the caller's responsibility to call BLS12_381_init() once,
 * and call BLS12_381_free() when needed.
 */
void BLS12_381_init();
void BLS12_381_free();

/* An element of the field F_p^2, a + bu, where a,b belong to F_p. */
typedef struct fp2_elem {
    mpz_t a;
    mpz_t b;
} fp2_elem;

/* An element of the field F_p^6. */
typedef struct fp6_elem {
    fp2_elem *a;
    fp2_elem *b;
    fp2_elem *c;
} fp6_elem;

/**
 * An element of the field F_p^12. The result of the pairing is an element in
 * this field.
 */
typedef struct fp12_elem {
    fp6_elem *a;
    fp6_elem *b;
} fp12_elem;

void fp12_elem_init(fp12_elem *e);
void fp12_elem_set(fp12_elem *e1, const fp12_elem *e2);
void fp12_elem_from_str(fp12_elem *e, const char *a[6], const char *b[6]);
void fp12_elem_free(fp12_elem *e);

/**
 * Arithmetic operations in the field extension F_p^12. This field is implemented
 * as a quadratic extension of F_p^6.
 *
 * By convention, the result of the operation is stored in the first argument. It
 * is safe to call fp12_ functions where the result and one or more of the operands
 * point to the same memory.
 */
void fp12_add(fp12_elem *x, const fp12_elem *y, const fp12_elem *z);
void fp12_sub(fp12_elem *x, const fp12_elem *y, const fp12_elem *z);
void fp12_mul(fp12_elem *x, const fp12_elem *y, const fp12_elem *z);
void fp12_square(fp12_elem *x, const fp12_elem *y);
void fp12_pow(fp12_elem *x, const fp12_elem *y, const mpz_t exp);
void fp12_inv(fp12_elem *x, const fp12_elem *y);
void fp12_conjugate(fp12_elem *x, const fp12_elem *y);
void fp12_frobenius(fp12_elem *x, const fp12_elem *y);
bool fp12_equal(const fp12_elem *e1, const fp12_elem *e2);

/**
 * G_1 is the group of r-torsion points of the curve E: y^2 = x^3 + 4 which are
 * contained in E(F_p). In other words, G_1 consists of the r-torsion which are
 * fixed by the p-power Frobenius.
 *
 * Below are operations and utilities for working with elements in BLS12-381's
 * G_1 group.  Generally, each type and function will have an _affine and _proj
 * variant, with the exception of the group operations themselves which are
 * implemented with algorithms using projective coordinates.
 *
 * Generally, one should use G1_elem_affine to store points in G_1 because affine
 * coorindates require less space. The G1_affine2proj and G1_proj2affine helpers
 * will convert to and from each coordinate space when needed.
 */
typedef struct {
    mpz_t x;
    mpz_t y;

    /* Indicate this is a point at infinity. */
    bool infinity;
} G1_elem_affine;

typedef struct {
    mpz_t x;
    mpz_t y;
    mpz_t z;
} G1_elem_proj;

/* Covert affine to projective, and vice versa. */
void G1_affine2proj(G1_elem_proj *proj, const G1_elem_affine *affn);
void G1_proj2affine(G1_elem_affine *affn, const G1_elem_proj *proj);

/**
 * Initialize a group element from the provided strings. Using these
 * functions does not guarantee that the resulting point will be on
 * the curve E. These should only be used when the input is known to
 * produce a group element, such as initializing generators from well
 * known parameters.
 */
void G1_elem_affine_from_str(G1_elem_affine *P, const char *x, const char *y);
void G1_elem_proj_from_str(G1_elem_proj *P, const char *x, const char *y, const char *z);

/**
 * Initialize a group element as the identity, i.e. the point at infinity.
 */
void G1_identity_init_affine(G1_elem_affine *e);
void G1_identity_init_proj(G1_elem_proj *e);

/**
 * Initialize a group element as a particular generator, namely the base point
 * using the coordinates defined above.
 *
 * Note that because G_1 is prime order, all points are generators.
 */
void G1_generator_init_affine(G1_elem_affine *g);
void G1_generator_init_proj(G1_elem_proj *g);

/**
 * Free the memory held by the group element data structure.
 */
void G1_elem_free_affine(G1_elem_affine *P);
void G1_elem_free_proj(G1_elem_proj *P);

/**
 * Make sure the point lies on the curve E. This should always be true
 * when using points generated by this API. This may NOT be true for
 * points initialized by G1_elem_affine_from_str or G1_elem_proj_from_str.
 */
bool G1_is_on_curve_affine(const G1_elem_affine *P);
bool G1_is_on_curve_proj(const G1_elem_proj *P);

/**
 * Check if the element P is the identity element.
 */
bool G1_is_identity_affine(const G1_elem_affine *P);
bool G1_is_identity_proj(const G1_elem_proj *P);

/**
 * Check if the elements P and Q are in the same equivalency class.
 */
bool G1_equiv_affine(const G1_elem_affine *P, const G1_elem_affine *Q);
bool G1_equiv_proj(const G1_elem_proj *P, const G1_elem_proj *Q);

/**
 * Negate an element in G_1.
 */
void G1_negate_affine(G1_elem_affine *r, const G1_elem_affine *P);
void G1_negate_proj(G1_elem_proj *r, const G1_elem_proj *P);

/**
 * Group operations.
 *
 * The algorithms are all implemented using projective coordinates, but
 * G1_mul_scalar uses affine parameters for convenience to the caller.
 */
void G1_double_proj(G1_elem_proj *r, const G1_elem_proj *P);
void G1_add_proj(G1_elem_proj *r, const G1_elem_proj *P, const G1_elem_proj *Q);
void G1_add_mixed(G1_elem_proj *r, const G1_elem_proj *P, const G1_elem_affine *Q);
void G1_mul_scalar(G1_elem_affine *r, const G1_elem_affine *P, const mpz_t m);

/**
 * G_2 is an order-r group of points on the twist E': y^2 = x^3 + 4(u + 1),
 * which is isomorphic to the so-called "trace zero subgroup" of the r-torison
 * of E. The full r-torsion is isomorphic to G_1 x G_2.
 *
 * Below are operations and utilities for working with elements in BLS12-381's
 * G_2 group.  G_2 is an order r subgroup of the curve E': y^2 = x^3 + 4(u+1)
 * defined over the field F_p^2. Generally, each type and function will have an
 * _affine and _proj variant, with the exception of the group operations
 * themselves which are implemented with algorithms using projective
 * coordinates.
 *
 * Generally, one should use G2_elem_affine to store points in G_1 because affine
 * coorindates require less space. The G2_affine2proj and G2_proj2affine helpers
 * will convert to and from each coordinate space when needed.
 */
typedef struct {
    fp2_elem *x;
    fp2_elem *y;

    /* Indicate this is a point at infinity. */
    bool infinity;
} G2_elem_affine;

typedef struct {
    fp2_elem *x;
    fp2_elem *y;
    fp2_elem *z;
} G2_elem_proj;

/* Covert affine to projective, and vice versa. */
void G2_affine2proj(G2_elem_proj *proj, const G2_elem_affine *affn);
void G2_proj2affine(G2_elem_affine *affn, const G2_elem_proj *proj);

/**
 * Initialize a group element from the provided strings. Using these
 * functions does not guarantee that the resulting point will lie on
 * the curve E. These should only be used when the input is known to
 * produce a group element, such as initializing generators from well
 * known parameters.
 */
void G2_elem_affine_from_str(G2_elem_affine *P,
                             const char *x0, const char *x1,
                             const char *y0, const char *y1);
void G2_elem_proj_from_str(G2_elem_proj *P,
                           const char *x0, const char *x1,
                           const char *y0, const char *y1,
                           const char *z0, const char *z1);

/**
 * Initialize a group element as the identity, i.e. the point at infinity.
 */
void G2_identity_init_affine(G2_elem_affine *e);
void G2_identity_init_proj(G2_elem_proj *e);

/**
 * Initialize a group element as a particular generator, namely the base point
 * using the coordinates defined above.
 *
 * Note that because G_2 is prime order, all points are generators.
 */
void G2_generator_init_affine(G2_elem_affine *g);
void G2_generator_init_proj(G2_elem_proj *g);

/**
 * Free the memory held by the group element data structure.
 */
void G2_elem_free_affine(G2_elem_affine *P);
void G2_elem_free_proj(G2_elem_proj *P);

/**
 * Make sure the point lies on the curve E'. This should always be true
 * when using points generated by this API. This may NOT be true for
 * points initialized by G2_elem_affine_from_str or G2_elem_proj_from_str.
 */
bool G2_is_on_curve_affine(const G2_elem_affine *P);
bool G2_is_on_curve_proj(const G2_elem_proj *P);

/**
 * Check if the element P is the identity element.
 */
bool G2_is_identity_affine(const G2_elem_affine *P);
bool G2_is_identity_proj(const G2_elem_proj *P);

/**
 * Check if the elements P and Q are in the same equivalency class.
 */
bool G2_equiv_affine(const G2_elem_affine *P, const G2_elem_affine *Q);
bool G2_equiv_proj(const G2_elem_proj *P, const G2_elem_proj *Q);

/**
 * Negate an element in G_2.
 */
void G2_negate_affine(G2_elem_affine *r, const G2_elem_affine *P);
void G2_negate_proj(G2_elem_proj *r, const G2_elem_proj *P);

/**
 * Group operations.
 *
 * The algorithms are all implemented using projective coordinates, but
 * G2_mul_scalar uses affine parameters for convenience to the caller.
 */
void G2_double_proj(G2_elem_proj *r, const G2_elem_proj *P);
void G2_add_proj(G2_elem_proj *r, const G2_elem_proj *P, const G2_elem_proj *Q);
void G2_add_mixed(G2_elem_proj *r, const G2_elem_proj *P, const G2_elem_affine *Q);
void G2_mul_scalar(G2_elem_affine *r, const G2_elem_affine *P, const mpz_t m);

/**
 * The pairing computation, given by the twisted ate pairing e: G_2 x G_1 -> G_T. This
 * computation includes the final exponentiation.
 */
void BLS12_381_pairing(fp12_elem *r, const G2_elem_affine *Q, const G1_elem_affine *P);

#endif /* BLS12_381_H */
