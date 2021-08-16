#include <stdbool.h>
#include <stdlib.h>
#include <gmp.h>

#include "finite_field.h"

#include "BLS12_381.h"

void G2_elem_set_affine(G2_elem_affine *P, const G2_elem_affine *Q)
{
    fp2_elem_set(P->x, Q->x);
    fp2_elem_set(P->y, Q->y);
    P->infinity = Q->infinity;
}

void G2_elem_set_proj(G2_elem_proj *P, const G2_elem_proj *Q)
{
    fp2_elem_set(P->x, Q->x);
    fp2_elem_set(P->y, Q->y);
    fp2_elem_set(P->z, Q->z);
}

void G2_affine2proj(G2_elem_proj *proj, const G2_elem_affine *affn)
{
    if (affn->infinity) {
        fp2_elem_set_si(proj->x, 0, 0);
        fp2_elem_set_si(proj->y, 1, 0);
        fp2_elem_set_si(proj->z, 0, 0);

        return;
    }

    fp2_elem_set(proj->x, affn->x);
    fp2_elem_set(proj->y, affn->y);
    fp2_elem_set_si(proj->z, 1, 0);
}

void G2_proj2affine(G2_elem_affine *affn, const G2_elem_proj *proj)
{
    if (G2_is_identity_proj(proj)) {
        fp2_elem_set_si(affn->x, 0, 0);
        fp2_elem_set_si(affn->y, 1, 0);
        affn->infinity = true;

        return;
    }

    affn->infinity = false;

    fp2_elem zinv;

    fp2_elem_init(&zinv);

    fp2_inv(&zinv, proj->z);
    fp2_mul(affn->x, proj->x, &zinv);
    fp2_mul(affn->y, proj->y, &zinv);

    fp2_elem_free(&zinv);
}

void G2_elem_affine_from_str(G2_elem_affine *P,
                             const char *x0, const char *x1,
                             const char *y0, const char *y1)
{
    P->x = calloc(1, sizeof(fp2_elem));
    P->y = calloc(1, sizeof(fp2_elem));

    fp2_elem_from_str(P->x, x0, x1);
    fp2_elem_from_str(P->y, y0, y1);
    P->infinity = false;
}

void G2_elem_proj_from_str(G2_elem_proj *P,
                           const char *x0, const char *x1,
                           const char *y0, const char *y1,
                           const char *z0, const char *z1)
{
    P->x = calloc(1, sizeof(fp2_elem));
    P->y = calloc(1, sizeof(fp2_elem));
    P->z = calloc(1, sizeof(fp2_elem));

    fp2_elem_from_str(P->x, x0, x1);
    fp2_elem_from_str(P->y, y0, y1);
    fp2_elem_from_str(P->z, z0, z1);
}

void G2_identity_init_affine(G2_elem_affine *e)
{
    e->x = calloc(1, sizeof(fp2_elem));
    e->y = calloc(1, sizeof(fp2_elem));

    fp2_elem_init(e->x);
    fp2_elem_init(e->y);

    fp2_elem_set_si(e->x, 0, 0); // 0 + 0x
    fp2_elem_set_si(e->y, 1, 0); // 1 + 0x
    e->infinity = true;
}

void G2_identity_init_proj(G2_elem_proj *e)
{
    e->x = calloc(1, sizeof(fp2_elem));
    e->y = calloc(1, sizeof(fp2_elem));
    e->z = calloc(1, sizeof(fp2_elem));

    fp2_elem_init(e->x);
    fp2_elem_init(e->y);
    fp2_elem_init(e->z);

    fp2_elem_set_si(e->x, 0, 0); // 0 + 0x
    fp2_elem_set_si(e->y, 1, 0); // 1 + 0x
    fp2_elem_set_si(e->z, 0, 0); // 0 + 0x
}

void G2_generator_init_affine(G2_elem_affine *g)
{
    g->x = calloc(1, sizeof(fp2_elem));
    g->y = calloc(1, sizeof(fp2_elem));

    fp2_elem_from_str(g->x, BLS12_381_G2_GEN_X0, BLS12_381_G2_GEN_X1);
    fp2_elem_from_str(g->y, BLS12_381_G2_GEN_Y0, BLS12_381_G2_GEN_Y1);
    g->infinity = false;
}

void G2_generator_init_proj(G2_elem_proj *g)
{
    g->x = calloc(1, sizeof(fp2_elem));
    g->y = calloc(1, sizeof(fp2_elem));
    g->z = calloc(1, sizeof(fp2_elem));

    fp2_elem_from_str(g->x, BLS12_381_G2_GEN_X0, BLS12_381_G2_GEN_X1);
    fp2_elem_from_str(g->y, BLS12_381_G2_GEN_Y0, BLS12_381_G2_GEN_Y1);

    fp2_elem_init(g->z);
    fp2_elem_set_si(g->z, 1, 0);
}

void G2_elem_free_affine(G2_elem_affine *P)
{
    fp2_elem_free(P->x);
    fp2_elem_free(P->y);
}

void G2_elem_free_proj(G2_elem_proj *P)
{
    fp2_elem_free(P->x);
    fp2_elem_free(P->y);
    fp2_elem_free(P->z);
}

bool G2_is_on_curve_affine(const G2_elem_affine *P)
{
    /* The identity is on the curve, return early. */
    if (P->infinity) {
        return true;
    }

    /* Check y^2 == x^3 + b' */
    bool ret;
    fp2_elem lhs, rhs, tmp;

    fp2_elem_init(&lhs);
    fp2_elem_init(&rhs);
    fp2_elem_init(&tmp);

    fp2_square(&lhs, P->y);

    // x^3
    fp2_square(&rhs, P->x);
    fp2_mul(&rhs, &rhs, P->x);

    // 4(u+1) = 4 + 4u
    fp2_elem_set_si(&tmp, 4, 4);
    fp2_add(&rhs, &rhs, &tmp);

    ret = fp2_equal(&lhs, &rhs);

    fp2_elem_free(&lhs);
    fp2_elem_free(&rhs);
    fp2_elem_free(&tmp);

    return ret;
}

bool G2_is_on_curve_proj(const G2_elem_proj *P)
{
    /* The identity is on the curve, return early. */
    if (G2_is_identity_proj(P)) {
        return true;
    }

    /* Check y^2 * z == x^3 + b'z^3 */
    bool ret;
    fp2_elem lhs, rhs, tmp1, tmp2;

    fp2_elem_init(&lhs);
    fp2_elem_init(&rhs);
    fp2_elem_init(&tmp1);
    fp2_elem_init(&tmp2);

    // y^2 * z
    fp2_square(&lhs, P->y);
    fp2_mul(&lhs, &lhs, P->z);

    // x^3
    fp2_square(&rhs, P->x);
    fp2_mul(&rhs, &rhs, P->x);

    // z^3
    fp2_square(&tmp1, P->z);
    fp2_mul(&tmp1, &tmp1, P->z);

    // b' = 4(u+1) = 4 + 4u
    fp2_elem_set_si(&tmp2, 4, 4);

    // b'z^3
    fp2_mul(&tmp1, &tmp1, &tmp2);

    fp2_add(&rhs, &rhs, &tmp1);

    ret = fp2_equal(&lhs, &rhs);

    fp2_elem_free(&lhs);
    fp2_elem_free(&rhs);
    fp2_elem_free(&tmp1);
    fp2_elem_free(&tmp2);

    return ret;
}

bool G2_is_identity_affine(const G2_elem_affine *P)
{
    return P->infinity;
}

bool G2_is_identity_proj(const G2_elem_proj *P)
{
    /**
     * In projective coordinates, a point is at infinity if P ~ (0 : 1 : 0).
     * On this curve, setting z = 0 yields x^3 = 0. Check if z == 0 to determine
     * if this point the identity (point at infinity).
     */
    return (mpz_cmp_si(P->z->a, 0) == 0) && (mpz_cmp_si(P->z->b, 0) == 0);
}

bool G2_equiv_affine(const G2_elem_affine *P, const G2_elem_affine *Q)
{
    if (P->infinity && Q->infinity) {
        return true;
    }

    if (P->infinity || Q->infinity) {
        /* We know from above this must be an exclusive or case. */
        return false;
    }

    return fp2_equal(P->x, Q->x) && fp2_equal(P->y, Q->y);
}

bool G2_equiv_proj(const G2_elem_proj *P, const G2_elem_proj *Q)
{
    /**
     * Check if (x1 : y1 : z1) ~ (x2 : y2 : z2). Do this by checking
     * z1*x2 == z2*x1 && z1*y2 == z2*y1. Or, if they are both the
     * identity.
     */
    if (G2_is_identity_proj(P) && G2_is_identity_proj(Q)) {
        return true;
    }

    if (G2_is_identity_proj(P) || G2_is_identity_proj(Q)) {
        /* We know from above this must be an exclusive or case. */
        return false;
    }

    bool ret;
    fp2_elem x_lhs, x_rhs, y_lhs, y_rhs;

    fp2_elem_init(&x_lhs);
    fp2_elem_init(&x_rhs);
    fp2_elem_init(&y_lhs);
    fp2_elem_init(&y_rhs);

    fp2_mul(&x_lhs, P->z, Q->x);
    fp2_mul(&x_rhs, Q->z, P->x);
    fp2_mul(&y_lhs, P->z, Q->y);
    fp2_mul(&y_rhs, Q->z, P->y);

    ret = (fp2_equal(&x_lhs, &x_rhs) && fp2_equal(&y_lhs, &y_rhs));

    fp2_elem_free(&x_lhs);
    fp2_elem_free(&x_rhs);
    fp2_elem_free(&y_lhs);
    fp2_elem_free(&y_rhs);

    return ret;
}

void G2_negate_affine(G2_elem_affine *r, const G2_elem_affine *P)
{
    fp2_elem_set(r->x, P->x);

    if (P->infinity) {
        fp2_elem_set(r->y, P->y);
        r->infinity = true;
    } else {
        mpz_neg(r->y->a, P->y->a);
        mpz_neg(r->y->b, P->y->b);
        r->infinity = false;
    }
}

void G2_negate_proj(G2_elem_proj *r, const G2_elem_proj *P)
{
    fp2_elem_set(r->x, P->x);
    mpz_neg(r->y->a, P->y->a);
    mpz_neg(r->y->b, P->y->b);
    fp2_elem_set(r->z, P->z);
}

void G2_add_proj(G2_elem_proj *r, const G2_elem_proj *P, const G2_elem_proj *Q)
{
    /* Algorithm 7: https://eprint.iacr.org/2015/1060.pdf */
    fp2_elem t0, t1, t2, t3, t4,
             rx, ry, rz,
             b3;

    /* If P or Q are the identity, we can return early. */
    if (G2_is_identity_proj(P)) {
        fp2_elem_set(r->x, Q->x);
        fp2_elem_set(r->y, Q->y);
        fp2_elem_set(r->z, Q->z);

        return;
    }

    if (G2_is_identity_proj(Q)) {
        fp2_elem_set(r->x, P->x);
        fp2_elem_set(r->y, P->y);
        fp2_elem_set(r->z, P->z);

        return;
    }

    fp2_elem_init(&t0);
    fp2_elem_init(&t1);
    fp2_elem_init(&t2);
    fp2_elem_init(&t3);
    fp2_elem_init(&t4);
    fp2_elem_init(&rx);
    fp2_elem_init(&ry);
    fp2_elem_init(&rz);
    fp2_elem_init(&b3);

    /**
     * b3 = 3*b', where b is the constant from the curve equation. It is 4(u+1)
     * in this case.
     */
    fp2_elem_set_si(&b3, 12, 12);

    fp2_mul(&t0, P->x, Q->x); // 1. t0 ← X1·X2
    fp2_mul(&t1, P->y, Q->y); // 2. t1 ← Y1·Y2
    fp2_mul(&t2, P->z, Q->z); // 3. t2 ← Z1·Z2
    fp2_add(&t3, P->x, P->y); // 4. t3 ← X1+Y1
    fp2_add(&t4, Q->x, Q->y); // 5. t4 ← X2+Y2
    fp2_mul(&t3, &t3, &t4);   // 6. t3 ← t3·t4
    fp2_add(&t4, &t0, &t1);   // 7. t4 ← t0+t1
    fp2_sub(&t3, &t3, &t4);   // 8. t3 ← t3−t4
    fp2_add(&t4, P->y, P->z); // 9. t4 ← Y1+Z1
    fp2_add(&rx, Q->y, Q->z); // 10. X3 ← Y2+Z2
    fp2_mul(&t4, &t4, &rx);   // 11. t4 ← t4·X3
    fp2_add(&rx, &t1, &t2);   // 12. X3 ← t1+t2
    fp2_sub(&t4, &t4, &rx);   // 13. t4 ← t4−X3
    fp2_add(&rx, P->x, P->z); // 14. X3 ← X1+Z1
    fp2_add(&ry, Q->x, Q->z); // 15. Y3 ← X2+Z2
    fp2_mul(&rx, &rx, &ry);   // 16. X3 ← X3·Y3
    fp2_add(&ry, &t0, &t2);   // 17. Y3 ← t0+t2
    fp2_sub(&ry, &rx, &ry);   // 18. Y3 ← X3−Y3
    fp2_add(&rx, &t0, &t0);   // 19. X3 ← t0+t0
    fp2_add(&t0, &rx, &t0);   // 20. t0 ← X3+t0
    fp2_mul(&t2, &b3, &t2);   // 21. t2 ← b3·t2
    fp2_add(&rz, &t1, &t2);   // 22. Z3 ← t1+t2
    fp2_sub(&t1, &t1, &t2);   // 23. t1 ← t1−t2
    fp2_mul(&ry, &b3, &ry);   // 24. Y3 ← b3·Y3
    fp2_mul(&rx, &t4, &ry);   // 25. X3 ← t4·Y3
    fp2_mul(&t2, &t3, &t1);   // 26. t2 ← t3·t1
    fp2_sub(&rx, &t2, &rx);   // 27. X3 ← t2−X3
    fp2_mul(&ry, &ry, &t0);   // 28. Y3 ← Y3·t0
    fp2_mul(&t1, &t1, &rz);   // 29. t1 ← t1·Z3
    fp2_add(&ry, &t1, &ry);   // 30. Y3 ← t1+Y3
    fp2_mul(&t0, &t0, &t3);   // 31. t0 ← t0·t3
    fp2_mul(&rz, &rz, &t4);   // 32. Z3 ← Z3·t4
    fp2_add(&rz, &rz, &t0);   // 33. Z3 ← Z3+t0

    fp2_elem_set(r->x, &rx);
    fp2_elem_set(r->y, &ry);
    fp2_elem_set(r->z, &rz);

    fp2_elem_free(&t0);
    fp2_elem_free(&t1);
    fp2_elem_free(&t2);
    fp2_elem_free(&t3);
    fp2_elem_free(&t4);
    fp2_elem_free(&rx);
    fp2_elem_free(&ry);
    fp2_elem_free(&rz);
    fp2_elem_free(&b3);
}

void G2_double_proj(G2_elem_proj *r, const G2_elem_proj *P)
{
    /* Algorithm 9: https://eprint.iacr.org/2015/1060.pdf */
    fp2_elem t0, t1, t2,
             rx, ry, rz,
             b3;

    /* If P is the identity, we can return early. */
    if (G2_is_identity_proj(P)) {
        fp2_elem_set(r->x, P->x);
        fp2_elem_set(r->y, P->y);
        fp2_elem_set(r->z, P->z);

        return;
    }

    fp2_elem_init(&t0);
    fp2_elem_init(&t1);
    fp2_elem_init(&t2);
    fp2_elem_init(&rx);
    fp2_elem_init(&ry);
    fp2_elem_init(&rz);
    fp2_elem_init(&b3);

    /**
     * b3 = 3*b', where b is the constant from the curve equation. It is 4(u+1)
     * in this case.
     */
    fp2_elem_set_si(&b3, 12, 12);

    fp2_mul(&t0, P->y, P->y); // 1. t0 ← Y·Y
    fp2_add(&rz, &t0, &t0);   // 2. Z3 ← t0+t0
    fp2_add(&rz, &rz, &rz);   // 3. Z3 ← Z3+Z3
    fp2_add(&rz, &rz, &rz);   // 4. Z3 ← Z3+Z3
    fp2_mul(&t1, P->y, P->z); // 5. t1 ← Y·Z
    fp2_mul(&t2, P->z, P->z); // 6. t2 ← Z·Z
    fp2_mul(&t2, &b3, &t2);   // 7. t2 ← b3·t2
    fp2_mul(&rx, &t2, &rz);   // 8. X3 ← t2·Z3
    fp2_add(&ry, &t0, &t2);   // 9. Y3 ← t0+t2
    fp2_mul(&rz, &t1, &rz);   // 10. Z3 ← t1·Z3
    fp2_add(&t1, &t2, &t2);   // 11. t1 ← t2+t2
    fp2_add(&t2, &t1, &t2);   // 12. t2 ← t1+t2
    fp2_sub(&t0, &t0, &t2);   // 13. t0 ← t0−t2
    fp2_mul(&ry, &t0, &ry);   // 14. Y3 ← t0·Y3
    fp2_add(&ry, &rx, &ry);   // 15. Y3 ← X3+Y3
    fp2_mul(&t1, P->x, P->y); // 16. t1 ← X·Y
    fp2_mul(&rx, &t0, &t1);   // 17. X3 ← t0·t1
    fp2_add(&rx, &rx, &rx);   // 18. X3 ← X3+X3

    fp2_elem_set(r->x, &rx);
    fp2_elem_set(r->y, &ry);
    fp2_elem_set(r->z, &rz);

    fp2_elem_free(&t0);
    fp2_elem_free(&t1);
    fp2_elem_free(&t2);
    fp2_elem_free(&rx);
    fp2_elem_free(&ry);
    fp2_elem_free(&rz);
    fp2_elem_free(&b3);
}

void G2_add_mixed(G2_elem_proj *r, const G2_elem_proj *P, const G2_elem_affine *Q)
{
    /* Algorithm 8: https://eprint.iacr.org/2015/1060.pdf */
    fp2_elem t0, t1, t2, t3, t4,
             rx, ry, rz,
             b3;

    /* If P or Q are the identity, we can return early. */
    if (G2_is_identity_proj(P)) {
        G2_affine2proj(r, Q);

        return;
    }

    if (G2_is_identity_affine(Q)) {
        fp2_elem_set(r->x, P->x);
        fp2_elem_set(r->y, P->y);
        fp2_elem_set(r->z, P->z);

        return;
    }

    fp2_elem_init(&t0);
    fp2_elem_init(&t1);
    fp2_elem_init(&t2);
    fp2_elem_init(&t3);
    fp2_elem_init(&t4);
    fp2_elem_init(&rx);
    fp2_elem_init(&ry);
    fp2_elem_init(&rz);
    fp2_elem_init(&b3);

    /**
     * b3 = 3*b', where b is the constant from the curve equation. It is 4(u+1)
     * in this case.
     */
    fp2_elem_set_si(&b3, 12, 12);

    fp2_mul(&t0, P->x, Q->x); // 1. t0 ← X1·X2
    fp2_mul(&t1, P->y, Q->y); // 2. t1 ← Y1·Y2
    fp2_add(&t3, Q->x, Q->y); // 3. t3 ← X2+Y2
    fp2_add(&t4, P->x, P->y); // 4. t4 ← X1+Y1
    fp2_mul(&t3, &t3, &t4);   // 5. t3 ← t3·t4
    fp2_add(&t4, &t0, &t1);   // 6. t4 ← t0+t1
    fp2_sub(&t3, &t3, &t4);   // 7. t3 ← t3−t4
    fp2_mul(&t4, Q->y, P->z); // 8. t4 ← Y2·Z1
    fp2_add(&t4, &t4, P->y);  // 9. t4 ← t4+Y1
    fp2_mul(&ry, Q->x, P->z); // 10. Y3 ← X2·Z1
    fp2_add(&ry, &ry, P->x);  // 11. Y3 ← Y3+X1
    fp2_add(&rx, &t0, &t0);   // 12. X3 ← t0+t0
    fp2_add(&t0, &rx, &t0);   // 13. t0 ← X3+t0
    fp2_mul(&t2, &b3, P->z);  // 14. t2 ← b3·Z1
    fp2_add(&rz, &t1, &t2);   // 15. Z3 ← t1+t2
    fp2_sub(&t1, &t1, &t2);   // 16. t1 ← t1−t2
    fp2_mul(&ry, &b3, &ry);   // 17. Y3 ← b3·Y3
    fp2_mul(&rx, &t4, &ry);   // 18. X3 ← t4·Y3
    fp2_mul(&t2, &t3, &t1);   // 19. t2 ← t3·t1
    fp2_sub(&rx, &t2, &rx);   // 20. X3 ← t2−X3
    fp2_mul(&ry, &ry, &t0);   // 21. Y3 ← Y3·t0
    fp2_mul(&t1, &t1, &rz);   // 22. t1 ← t1·Z3
    fp2_add(&ry, &t1, &ry);   // 23. Y3 ← t1+Y3
    fp2_mul(&t0, &t0, &t3);   // 24. t0 ← t0·t3
    fp2_mul(&rz, &rz, &t4);   // 25. Z3 ← Z3·t4
    fp2_add(&rz, &rz, &t0);   // 26. Z3 ← Z3+t0

    fp2_elem_set(r->x, &rx);
    fp2_elem_set(r->y, &ry);
    fp2_elem_set(r->z, &rz);

    fp2_elem_free(&t0);
    fp2_elem_free(&t1);
    fp2_elem_free(&t2);
    fp2_elem_free(&t3);
    fp2_elem_free(&t4);
    fp2_elem_free(&rx);
    fp2_elem_free(&ry);
    fp2_elem_free(&rz);
    fp2_elem_free(&b3);
}

void G2_mul_scalar(G2_elem_affine *r, const G2_elem_affine *P, const mpz_t m)
{
    G2_elem_proj r_proj, P_proj;
    mp_bitcnt_t c;

    G2_identity_init_proj(&P_proj);
    G2_identity_init_proj(&r_proj);

    G2_affine2proj(&P_proj, P);
    G2_affine2proj(&r_proj, P);

    c = (mp_bitcnt_t)mpz_sizeinbase(m, 2) - 2;
    for (;;) {
        G2_double_proj(&r_proj, &r_proj);

        if (mpz_tstbit(m, c)) {
            G2_add_proj(&r_proj, &r_proj, &P_proj);
        }

        if (!c) break;
        c--;
    }

    G2_proj2affine(r, &r_proj);

    G2_elem_free_proj(&P_proj);
    G2_elem_free_proj(&r_proj);
}
