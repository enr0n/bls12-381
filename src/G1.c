#include <stdbool.h>
#include <gmp.h>

#include "finite_field.h"

#include "BLS12_381.h"

void G1_affine2proj(G1_elem_proj *proj, const G1_elem_affine *affn)
{
    if (affn->infinity) {
        mpz_set_ui(proj->x, 0);
        mpz_set_ui(proj->y, 1);
        mpz_set_ui(proj->z, 0);

        return;
    }

    mpz_set(proj->x, affn->x);
    mpz_set(proj->y, affn->y);
    mpz_set_si(proj->z, 1);
}

void G1_proj2affine(G1_elem_affine *affn, const G1_elem_proj *proj)
{
    if (G1_is_identity_proj(proj)) {
        mpz_set_si(affn->x, 0);
        mpz_set_si(affn->y, 1);
        affn->infinity = true;

        return;
    }

    affn->infinity = false;

    mpz_t zinv;

    mpz_init(zinv);

    fp_inv(zinv, proj->z);
    fp_mul(affn->x, proj->x, zinv);
    fp_mul(affn->y, proj->y, zinv);

    mpz_clear(zinv);
}

void G1_elem_affine_from_str(G1_elem_affine *P, const char *x, const char *y)
{
    mpz_init_set_str(P->x, x, 0);
    mpz_init_set_str(P->y, y, 0);
    P->infinity = false;
}

void G1_elem_proj_from_str(G1_elem_proj *P, const char *x, const char *y, const char *z)
{
    mpz_init_set_str(P->x, x, 0);
    mpz_init_set_str(P->y, y, 0);
    mpz_init_set_str(P->z, z, 0);
}

void G1_identity_init_affine(G1_elem_affine *e)
{
    mpz_init_set_si(e->x, 0);
    mpz_init_set_si(e->y, 1);
    e->infinity = true;
}

void G1_identity_init_proj(G1_elem_proj *e)
{
    mpz_init_set_si(e->x, 0);
    mpz_init_set_si(e->y, 1);
    mpz_init_set_si(e->z, 0);
}

void G1_generator_init_affine(G1_elem_affine *g)
{
    mpz_init_set_str(g->x, BLS12_381_G1_GEN_X, 0);
    mpz_init_set_str(g->y, BLS12_381_G1_GEN_Y, 0);
    g->infinity = false;
}

void G1_generator_init_proj(G1_elem_proj *g)
{
    mpz_init_set_str(g->x, BLS12_381_G1_GEN_X, 0);
    mpz_init_set_str(g->y, BLS12_381_G1_GEN_Y, 0);
    mpz_init_set_si(g->z, 1);
}

void G1_elem_free_affine(G1_elem_affine *P)
{
    mpz_clears(P->x, P->y, NULL);
}

void G1_elem_free_proj(G1_elem_proj *P)
{
    mpz_clears(P->x, P->y, P->z, NULL);
}

bool G1_is_on_curve_affine(const G1_elem_affine *P)
{
    /* The identity is on the curve, return early. */
    if (P->infinity) {
        return true;
    }

    /* Check y^2 == x^3 + b */
    bool ret;
    mpz_t lhs, rhs;

    mpz_inits(lhs, rhs, NULL);

    mpz_pow_ui(lhs, P->y, 2);

    mpz_pow_ui(rhs, P->x, 3);
    mpz_add_ui(rhs, rhs, 4);

    ret = fp_equiv(lhs, rhs);

    mpz_clears(lhs, rhs, NULL);

    return ret;
}

bool G1_is_on_curve_proj(const G1_elem_proj *P)
{
    /* The identity is on the curve, return early. */
    if (G1_is_identity_proj(P)) {
        return true;
    }

    /* Check y^2 * z == x^3 + bz^3 */
    bool ret;
    mpz_t lhs, rhs, tmp;

    mpz_inits(lhs, rhs, tmp, NULL);

    mpz_pow_ui(lhs, P->y, 2);
    mpz_mul(lhs, lhs, P->z);

    mpz_pow_ui(rhs, P->x, 3);
    mpz_pow_ui(tmp, P->z, 3);
    mpz_mul_ui(tmp, tmp, 4);
    mpz_add(rhs, rhs, tmp);

    ret = fp_equiv(lhs, rhs);

    mpz_clears(lhs, rhs, tmp, NULL);

    return ret;
}

bool G1_is_identity_affine(const G1_elem_affine *P)
{
    return P->infinity;
}

bool G1_is_identity_proj(const G1_elem_proj *P)
{
    /**
     * In projective coordinates, a point is at infinity if P ~ (0 : 1 : 0).
     * On this curve, setting z = 0 yields x^3 = 0. Check if z == 0 to determine
     * if this point the identity (point at infinity).
     */
    return (mpz_cmp_si(P->z, 0) == 0);
}

bool G1_equiv_affine(const G1_elem_affine *P, const G1_elem_affine *Q)
{
    if (P->infinity && Q->infinity) {
        return true;
    }

    if (P->infinity || Q->infinity) {
        /* We know from above this must be an exclusive or case. */
        return false;
    }

    return fp_equiv(P->x, Q->x) && fp_equiv(P->y, Q->y);
}

bool G1_equiv_proj(const G1_elem_proj *P, const G1_elem_proj *Q)
{
    /**
     * Check if (x1 : y1 : z1) ~ (x2 : y2 : z2). Do this by checking
     * z1*x2 == z2*x1 && z1*y2 == z2*y1. Or, if they are both the
     * identity.
     */
    if (G1_is_identity_proj(P) && G1_is_identity_proj(Q)) {
        return true;
    }

    if (G1_is_identity_proj(P) || G1_is_identity_proj(Q)) {
        /* We know from above this must be an exclusive or case. */
        return false;
    }

    bool ret;
    mpz_t x_lhs, x_rhs, y_lhs, y_rhs;

    mpz_inits(x_lhs, x_rhs, y_lhs, y_rhs, NULL);

    fp_mul(x_lhs, P->z, Q->x);
    fp_mul(x_rhs, Q->z, P->x);
    fp_mul(y_lhs, P->z, Q->y);
    fp_mul(y_rhs, Q->z, P->y);

    ret = (fp_equiv(x_lhs, x_rhs) && fp_equiv(y_lhs, y_rhs));

    mpz_clears(x_lhs, x_rhs, y_lhs, y_rhs, NULL);

    return ret;
}

void G1_negate_affine(G1_elem_affine *r, const G1_elem_affine *P)
{
    mpz_set(r->x, P->x);

    if (P->infinity) {
        mpz_set(r->y, P->y);
        r->infinity = true;
    } else {
        mpz_neg(r->y, P->y);
        r->infinity = false;
    }
}

void G1_negate_proj(G1_elem_proj *r, const G1_elem_proj *P)
{
    mpz_set(r->x, P->x);
    mpz_neg(r->y, P->y);
    mpz_set(r->z, P->z);
}

void G1_add_proj(G1_elem_proj *r, const G1_elem_proj *P, const G1_elem_proj *Q)
{
    /* Algorithm 7: https://eprint.iacr.org/2015/1060.pdf */
    mpz_t t0, t1, t2, t3, t4,
          rx, ry, rz,
          b3;

    /* If P or Q are the identity, we can return early. */
    if (G1_is_identity_proj(P)) {
        mpz_set(r->x, Q->x);
        mpz_set(r->y, Q->y);
        mpz_set(r->z, Q->z);

        return;
    }

    if (G1_is_identity_proj(Q)) {
        mpz_set(r->x, P->x);
        mpz_set(r->y, P->y);
        mpz_set(r->z, P->z);

        return;
    }

    mpz_inits(t0, t1, t2, t3, t4, rx, ry, rz, NULL);

    /* b3 = 3*b, where b is the constant from the curve equation. It is 4 in this case. */
    mpz_init_set_si(b3, 12);

    fp_mul(t0, P->x, Q->x); // 1. t0 ← X1·X2
    fp_mul(t1, P->y, Q->y); // 2. t1 ← Y1·Y2
    fp_mul(t2, P->z, Q->z); // 3. t2 ← Z1·Z2
    fp_add(t3, P->x, P->y); // 4. t3 ← X1+Y1
    fp_add(t4, Q->x, Q->y); // 5. t4 ← X2+Y2
    fp_mul(t3, t3, t4);     // 6. t3 ← t3·t4
    fp_add(t4, t0, t1);     // 7. t4 ← t0+t1
    fp_sub(t3, t3, t4);     // 8. t3 ← t3−t4
    fp_add(t4, P->y, P->z); // 9. t4 ← Y1+Z1
    fp_add(rx, Q->y, Q->z); // 10. X3 ← Y2+Z2
    fp_mul(t4, t4, rx);     // 11. t4 ← t4·X3
    fp_add(rx, t1, t2);     // 12. X3 ← t1+t2
    fp_sub(t4, t4, rx);     // 13. t4 ← t4−X3
    fp_add(rx, P->x, P->z); // 14. X3 ← X1+Z1
    fp_add(ry, Q->x, Q->z); // 15. Y3 ← X2+Z2
    fp_mul(rx, rx, ry);     // 16. X3 ← X3·Y3
    fp_add(ry, t0, t2);     // 17. Y3 ← t0+t2
    fp_sub(ry, rx, ry);     // 18. Y3 ← X3−Y3
    fp_add(rx, t0, t0);     // 19. X3 ← t0+t0
    fp_add(t0, rx, t0);     // 20. t0 ← X3+t0
    fp_mul(t2, b3, t2);     // 21. t2 ← b3·t2
    fp_add(rz, t1, t2);     // 22. Z3 ← t1+t2
    fp_sub(t1, t1, t2);     // 23. t1 ← t1−t2
    fp_mul(ry, b3, ry);     // 24. Y3 ← b3·Y3
    fp_mul(rx, t4, ry);     // 25. X3 ← t4·Y3
    fp_mul(t2, t3, t1);     // 26. t2 ← t3·t1
    fp_sub(rx, t2, rx);     // 27. X3 ← t2−X3
    fp_mul(ry, ry, t0);     // 28. Y3 ← Y3·t0
    fp_mul(t1, t1, rz);     // 29. t1 ← t1·Z3
    fp_add(ry, t1, ry);     // 30. Y3 ← t1+Y3
    fp_mul(t0, t0, t3);     // 31. t0 ← t0·t3
    fp_mul(rz, rz, t4);     // 32. Z3 ← Z3·t4
    fp_add(rz, rz, t0);     // 33. Z3 ← Z3+t0

    mpz_set(r->x, rx);
    mpz_set(r->y, ry);
    mpz_set(r->z, rz);

    mpz_clears(t0, t1, t2, t3, t4, rx, ry, rz, b3, NULL);
}

void G1_double_proj(G1_elem_proj *r, const G1_elem_proj *P)
{
    /* Algorithm 9: https://eprint.iacr.org/2015/1060.pdf */
    mpz_t t0, t1, t2,
          rx, ry, rz,
          b3;

    /* If P is the identity, we can return early. */
    if (G1_is_identity_proj(P)) {
        mpz_set(r->x, P->x);
        mpz_set(r->y, P->y);
        mpz_set(r->z, P->z);

        return;
    }

    mpz_inits(t0, t1, t2, rx, ry, rz, NULL);

    /* b3 = 3*b, where b is the constant from the curve equation. It is 4 in this case. */
    mpz_init_set_si(b3, 12);

    fp_mul(t0, P->y, P->y); // 1. t0 ← Y·Y
    fp_add(rz, t0, t0);     // 2. Z3 ← t0+t0
    fp_add(rz, rz, rz);     // 3. Z3 ← Z3+Z3
    fp_add(rz, rz, rz);     // 4. Z3 ← Z3+Z3
    fp_mul(t1, P->y, P->z); // 5. t1 ← Y·Z
    fp_mul(t2, P->z, P->z); // 6. t2 ← Z·Z
    fp_mul(t2, b3, t2);     // 7. t2 ← b3·t2
    fp_mul(rx, t2, rz);     // 8. X3 ← t2·Z3
    fp_add(ry, t0, t2);     // 9. Y3 ← t0+t2
    fp_mul(rz, t1, rz);     // 10. Z3 ← t1·Z3
    fp_add(t1, t2, t2);     // 11. t1 ← t2+t2
    fp_add(t2, t1, t2);     // 12. t2 ← t1+t2
    fp_sub(t0, t0, t2);     // 13. t0 ← t0−t2
    fp_mul(ry, t0, ry);     // 14. Y3 ← t0·Y3
    fp_add(ry, rx, ry);     // 15. Y3 ← X3+Y3
    fp_mul(t1, P->x, P->y); // 16. t1 ← X·Y
    fp_mul(rx, t0, t1);     // 17. X3 ← t0·t1
    fp_add(rx, rx, rx);     // 18. X3 ← X3+X3

    mpz_set(r->x, rx);
    mpz_set(r->y, ry);
    mpz_set(r->z, rz);

    mpz_clears(t0, t1, t2, rx, ry, rz, b3, NULL);
}

void G1_add_mixed(G1_elem_proj *r, const G1_elem_proj *P, const G1_elem_affine *Q)
{
    /* Algorithm 8: https://eprint.iacr.org/2015/1060.pdf */
    mpz_t t0, t1, t2, t3, t4,
          rx, ry, rz,
          b3;

    /* If P or Q are the identity, we can return early. */
    if (G1_is_identity_proj(P)) {
        G1_affine2proj(r, Q);

        return;
    }

    if (G1_is_identity_affine(Q)) {
        mpz_set(r->x, P->x);
        mpz_set(r->y, P->y);
        mpz_set(r->z, P->z);

        return;
    }

    mpz_inits(t0, t1, t2, t3, t4, rx, ry, rz, NULL);

    /* b3 = 3*b, where b is the constant from the curve equation. It is 4 in this case. */
    mpz_init_set_si(b3, 12);

    fp_mul(t0, P->x, Q->x); // 1. t0 ← X1·X2
    fp_mul(t1, P->y, Q->y); // 2. t1 ← Y1·Y2
    fp_add(t3, Q->x, Q->y); // 3. t3 ← X2+Y2
    fp_add(t4, P->x, P->y); // 4. t4 ← X1+Y1
    fp_mul(t3, t3, t4);     // 5. t3 ← t3·t4
    fp_add(t4, t0, t1);     // 6. t4 ← t0+t1
    fp_sub(t3, t3, t4);     // 7. t3 ← t3−t4
    fp_mul(t4, Q->y, P->z); // 8. t4 ← Y2·Z1
    fp_add(t4, t4, P->y);   // 9. t4 ← t4+Y1
    fp_mul(ry, Q->x, P->z); // 10. Y3 ← X2·Z1
    fp_add(ry, ry, P->x);   // 11. Y3 ← Y3+X1
    fp_add(rx, t0, t0);     // 12. X3 ← t0+t0
    fp_add(t0, rx, t0);     // 13. t0 ← X3+t0
    fp_mul(t2, b3, P->z);   // 14. t2 ← b3·Z1
    fp_add(rz, t1, t2);     // 15. Z3 ← t1+t2
    fp_sub(t1, t1, t2);     // 16. t1 ← t1−t2
    fp_mul(ry, b3, ry);     // 17. Y3 ← b3·Y3
    fp_mul(rx, t4, ry);     // 18. X3 ← t4·Y3
    fp_mul(t2, t3, t1);     // 19. t2 ← t3·t1
    fp_sub(rx, t2, rx);     // 20. X3 ← t2−X3
    fp_mul(ry, ry, t0);     // 21. Y3 ← Y3·t0
    fp_mul(t1, t1, rz);     // 22. t1 ← t1·Z3
    fp_add(ry, t1, ry);     // 23. Y3 ← t1+Y3
    fp_mul(t0, t0, t3);     // 24. t0 ← t0·t3
    fp_mul(rz, rz, t4);     // 25. Z3 ← Z3·t4
    fp_add(rz, rz, t0);     // 26. Z3 ← Z3+t0

    mpz_set(r->x, rx);
    mpz_set(r->y, ry);
    mpz_set(r->z, rz);

    mpz_clears(t0, t1, t2, t3, t4, rx, ry, rz, b3, NULL);
}

void G1_mul_scalar(G1_elem_affine *r, const G1_elem_affine *P, const mpz_t m)
{
    G1_elem_proj r_proj, P_proj;
    mp_bitcnt_t c;

    G1_identity_init_proj(&P_proj);
    G1_identity_init_proj(&r_proj);

    G1_affine2proj(&P_proj, P);
    G1_affine2proj(&r_proj, P);

    c = (mp_bitcnt_t)mpz_sizeinbase(m, 2) - 2;
    for (;;) {
        G1_double_proj(&r_proj, &r_proj);

        if (mpz_tstbit(m, c)) {
            G1_add_proj(&r_proj, &r_proj, &P_proj);
        }

        if (!c) break;
        c--;
    }

    G1_proj2affine(r, &r_proj);

    G1_elem_free_proj(&P_proj);
    G1_elem_free_proj(&r_proj);
}
