#include <assert.h>
#include <stdio.h>

#include "finite_field.h"
#include "G1.h"

#define TEST_UNIMPL \
    printf("Test %s is not implemented!\n", __func__); \
    return false;

#define TEST_MAIN_INIT                  \
    int result, pass_count, fail_count; \
                                        \
    pass_count = fail_count = 0;

#define TEST_MAIN_RETURN                                                       \
    printf("G1 test results: %d passed, %d failed\n", pass_count, fail_count); \
                                                                               \
    return !!(fail_count);

#define TEST_RUN(test)                   \
    result = test;                       \
    if (!result) {                       \
        printf("FAIL\n\n");              \
        fail_count++;                    \
    } else {                             \
        printf("PASS\n\n");              \
        pass_count++;                    \
    }

bool test_is_identity()
{
    printf("Running %s...\n", __func__);
    bool b1, b2, b3, b4;
    G1_elem_proj e_proj, g_proj;
    G1_elem_affine e_affn, g_affn;

    G1_identity_init_proj(&e_proj);
    G1_identity_init_affine(&e_affn);
    G1_generator_init_proj(&g_proj);
    G1_generator_init_affine(&g_affn);

    b1 = G1_is_identity_proj(&e_proj);
    b2 = G1_is_identity_proj(&g_proj);
    b3 = G1_is_identity_affine(&e_affn);
    b4 = G1_is_identity_affine(&g_affn);

    printf("e_proj is identity: %d\n", b1);
    printf("g_proj is identity: %d\n", b2);
    printf("e_affn is identity: %d\n", b3);
    printf("g_affn is identity: %d\n", b4);

    G1_elem_free_proj(&e_proj);
    G1_elem_free_proj(&g_proj);
    G1_elem_free_affine(&e_affn);
    G1_elem_free_affine(&g_affn);

    return b1 && !b2 && b3 && !b4;
}

bool test_is_on_curve()
{
    printf("Running %s...\n", __func__);
    G1_elem_proj e_proj, g_proj, not_on_curve;
    G1_elem_affine e_affn, g_affn;

    G1_identity_init_proj(&e_proj);
    G1_identity_init_affine(&e_affn);
    G1_generator_init_proj(&g_proj);
    G1_generator_init_affine(&g_affn);
    G1_elem_proj_from_str(&not_on_curve,
        "0x0719ed2d3f59c903d00180853a54ca52f536e48a548eaa9ff02a64643a66be9a1476725ac1fa89831072ba6cedbe2101",
        "0x055dd7c08f4fc4a706f9d4f777af34c8fce3fcd442dd46c717a2db4c9d5162c5ae9d746af5d1d5ae03014edd1b8eaffe",
        "0x0719ed2d3f59c903d00180853a54ca52f536e48a548eaa9ff02a64643a66be9a1476725ac1fa89831072ba6cedbe2101"
    );

    assert(G1_is_on_curve_proj(&e_proj));
    assert(G1_is_on_curve_proj(&g_proj));
    assert(G1_is_on_curve_affine(&e_affn));
    assert(G1_is_on_curve_affine(&g_affn));
    assert(!G1_is_on_curve_proj(&not_on_curve));

    G1_elem_free_proj(&e_proj);
    G1_elem_free_proj(&g_proj);
    G1_elem_free_proj(&not_on_curve);
    G1_elem_free_affine(&e_affn);
    G1_elem_free_affine(&g_affn);

    return true;
}

bool test_affine2proj()
{
    printf("Running %s...\n", __func__);
    G1_elem_affine e_affn, g_affn;
    G1_elem_proj e_proj, g_proj, tmp;

    G1_identity_init_affine(&e_affn);
    G1_generator_init_affine(&g_affn);

    G1_elem_proj_from_str(&e_proj, "0x0", "0x0", "0x0");
    G1_elem_proj_from_str(&g_proj, "0x0", "0x0", "0x0");

    G1_affine2proj(&e_proj, &e_affn);
    G1_affine2proj(&g_proj, &g_affn);

    assert(G1_is_on_curve_proj(&e_proj));
    assert(G1_is_on_curve_proj(&g_proj));

    assert(G1_is_identity_proj(&e_proj));

    G1_generator_init_proj(&tmp);
    assert(G1_equiv_proj(&tmp, &g_proj));

    G1_elem_free_proj(&e_proj);
    G1_elem_free_proj(&g_proj);
    G1_elem_free_proj(&tmp);
    G1_elem_free_affine(&e_affn);
    G1_elem_free_affine(&g_affn);

    return true;
}

bool test_proj2affine()
{
    printf("Running %s...\n", __func__);
    G1_elem_affine e_affn, g_affn, tmp;
    G1_elem_proj e_proj, g_proj;

    G1_identity_init_proj(&e_proj);
    G1_generator_init_proj(&g_proj);

    G1_elem_affine_from_str(&e_affn, "0x0", "0x0");
    G1_elem_affine_from_str(&g_affn, "0x0", "0x0");

    G1_proj2affine(&e_affn, &e_proj);
    G1_proj2affine(&g_affn, &g_proj);

    assert(G1_is_on_curve_affine(&e_affn));
    assert(G1_is_identity_affine(&e_affn));

    G1_generator_init_affine(&tmp);

    assert(G1_equiv_affine(&tmp, &g_affn));
    assert(G1_is_on_curve_affine(&g_affn));

    G1_elem_free_proj(&e_proj);
    G1_elem_free_proj(&g_proj);
    G1_elem_free_affine(&tmp);
    G1_elem_free_affine(&e_affn);
    G1_elem_free_affine(&g_affn);

    return true;
}

bool test_add_proj()
{
    printf("Running %s...\n", __func__);
    bool ret;
    G1_elem_proj P, actual, expect;

    G1_elem_proj_from_str(&P,
        "0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb",
        "0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1",
        "0x1"
    );

    G1_elem_proj_from_str(&expect,
        "0x05dff4ac6726c6cb9b6d4dac3f33e92c062e48a6104cc52f6e7f23d4350c60bd7803e16723f9f1478a13c2b29f4325ad",
        "0x14e4b429606d02bc3c604c0410e5fc01d6093a00bb3e2bc9395952af0b6a0dbd599a8782a1bea48a2aa4d8e1b1df7ca5",
        "0x0430df56ea4aba6928180e61b1f2cb8f962f5650798fdf279a55bee62edcdb27c04c720ae01952ac770553ef06aadf22"
    );

    G1_identity_init_proj(&actual);

    G1_add_proj(&actual, &actual, &actual);
    assert(G1_is_on_curve_proj(&actual));
    assert(G1_is_identity_proj(&actual));

    G1_add_proj(&actual, &P, &actual);
    assert(G1_is_on_curve_proj(&actual));
    assert(!G1_is_identity_proj(&actual));

    ret = G1_equiv_proj(&actual, &P);
    assert(ret);

    G1_add_proj(&actual, &P, &P);

    ret = G1_equiv_proj(&actual, &expect);
    printf("Actual:\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, actual.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, actual.y));
    printf("\tz: %s\n", mpz_get_str(NULL, 16, actual.z));
    printf("Expected:\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, expect.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, expect.y));
    printf("\tz: %s\n", mpz_get_str(NULL, 16, expect.z));

    G1_elem_free_proj(&P);
    G1_elem_free_proj(&actual);
    G1_elem_free_proj(&expect);

    return ret;
}

bool test_add_mixed()
{
    printf("Running %s...\n", __func__);
    bool ret;
    G1_elem_affine Q;
    G1_elem_proj P, actual, expect;

    G1_elem_proj_from_str(&P,
        "0x191b633a5aaef4b0de0c66c95e9f5a7218bf0dde55cf16e4c600f5d7308e1c7d974471889fad50f2f244e8edfb988ba4",
        "0x057810823f03699b4e1497b3691fa7b1c2deb8cd5ba748927e921b3a2d19f2c4eae75664cc7c5f76f65b24b470dbdc51",
        "0x0d397f2ba295771076a85b84bafd4a721237631be7721cd89727e128d5b6705cac7db029cb9604ff6884c19b42bb871a"
    );

    G1_elem_affine_from_str(&Q,
        "0x0572cbea904d67468808c8eb50a9450c9721db309128012543902d0ac358a62ae28f75bb8f1c7c42c39a8c5529bf0f4e",
        "0x166a9d8cabc673a322fda673779d8e3822ba3ecb8670e461f73bb9021d5fd76a4c56d9d4cd16bd1bba86881979749d28"
    );

    G1_elem_proj_from_str(&expect,
        "0x11fe1d770cdf4bd88a845b2688c95b2e47bc9a78518165f044d6eb1faa217b67daee27d5e26fc4bb21423d0b80bde83f",
        "0x181e9afd2aa89c9f99f1ca7dba34645eb3996ea770650feb7796748c99e28a4a99550a809fed9fe800dbb32cce6f98e3",
        "0x077c8a7ac3188fa2ab3558545d0753c78bfc7491d5d16b7c979edf8e30f6ffe361ec68a7d0861165b1cf3cac19261931"
    );

    G1_identity_init_proj(&actual);

    assert(G1_is_on_curve_proj(&P));
    assert(!G1_is_identity_proj(&P));
    assert(G1_is_on_curve_affine(&Q));
    assert(!G1_is_identity_affine(&Q));

    G1_add_mixed(&actual, &P, &Q);

    assert(G1_is_on_curve_proj(&actual));
    assert(!G1_is_identity_proj(&actual));

    ret = G1_equiv_proj(&actual, &expect);
    printf("Actual (NOT NORMALIZED):\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, actual.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, actual.y));
    printf("\tz: %s\n", mpz_get_str(NULL, 16, actual.z));
    printf("Expected (NOT NORMALIZED):\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, expect.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, expect.y));
    printf("\tz: %s\n", mpz_get_str(NULL, 16, expect.z));

    G1_elem_free_proj(&P);
    G1_elem_free_affine(&Q);
    G1_elem_free_proj(&actual);
    G1_elem_free_proj(&expect);

    return ret;
}

bool test_double_proj()
{
    printf("Running %s...\n", __func__);
    bool ret;
    G1_elem_proj P, actual, expect;

    G1_elem_proj_from_str(&P,
        "0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb",
        "0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1",
        "0x1"
    );

    G1_elem_proj_from_str(&expect,
        "0x05dff4ac6726c6cb9b6d4dac3f33e92c062e48a6104cc52f6e7f23d4350c60bd7803e16723f9f1478a13c2b29f4325ad",
        "0x14e4b429606d02bc3c604c0410e5fc01d6093a00bb3e2bc9395952af0b6a0dbd599a8782a1bea48a2aa4d8e1b1df7ca5",
        "0x0430df56ea4aba6928180e61b1f2cb8f962f5650798fdf279a55bee62edcdb27c04c720ae01952ac770553ef06aadf22"
    );

    G1_identity_init_proj(&actual);

    G1_double_proj(&actual, &actual);
    assert(G1_is_on_curve_proj(&actual));
    assert(G1_is_identity_proj(&actual));

    G1_double_proj(&actual, &P);

    assert(G1_is_on_curve_proj(&actual));
    assert(!G1_is_identity_proj(&actual));

    ret = G1_equiv_proj(&actual, &expect);
    printf("Actual:\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, actual.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, actual.y));
    printf("\tz: %s\n", mpz_get_str(NULL, 16, actual.z));
    printf("Expected:\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, expect.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, expect.y));
    printf("\tz: %s\n", mpz_get_str(NULL, 16, expect.z));

    G1_elem_free_proj(&P);
    G1_elem_free_proj(&actual);
    G1_elem_free_proj(&expect);

    return ret;
}

bool test_mul_scalar()
{
    printf("Running %s...\n", __func__);
    bool ret;
    G1_elem_affine P, actual, expect;
    mpz_t m;

    mpz_init_set_str(m, "0x0c9af8fcfd18dc6a102260d25e1de7918505fac80ef519d22b568298a56da71b", 0);
    G1_elem_affine_from_str(&expect,
        "0x055877e8edd635641c0b869e39903b0324ad8fc53736fc7d7f9795a539f556d0f19281364f22a97073172e9fe9cc9c7a",
        "0x08539abef12abbecb497cf5eec349803e8b7afe5a7acdad2e4a8fe8519fdf4274f4b19438081db80a0275bc70f822c25"
    );

    G1_generator_init_affine(&P);
    G1_identity_init_affine(&actual);

    G1_mul_scalar(&actual, &P, m);

    ret = G1_equiv_affine(&actual, &expect);
    printf("Actual:\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, actual.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, actual.y));
    printf("Expected:\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, expect.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, expect.y));

    G1_elem_free_affine(&P);
    G1_elem_free_affine(&actual);
    G1_elem_free_affine(&expect);
    mpz_clear(m);

    return ret;
}

int main() {
    TEST_MAIN_INIT
    fp_params_init();

    TEST_RUN(test_is_identity());
    TEST_RUN(test_is_on_curve());
    TEST_RUN(test_double_proj());
    TEST_RUN(test_add_proj());
    TEST_RUN(test_add_mixed());
    TEST_RUN(test_affine2proj());
    TEST_RUN(test_proj2affine());
    TEST_RUN(test_mul_scalar());

    fp_params_free();
    TEST_MAIN_RETURN
}
