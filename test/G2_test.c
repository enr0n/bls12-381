#include <assert.h>
#include <stdio.h>

#include "finite_field.h"
#include "BLS12_381.h"

#define TEST_UNIMPL \
    printf("Test %s is not implemented!\n", __func__); \
    return false;

#define TEST_MAIN_INIT                  \
    int result, pass_count, fail_count; \
                                        \
    pass_count = fail_count = 0;

#define TEST_MAIN_RETURN                                                       \
    printf("G2 test results: %d passed, %d failed\n", pass_count, fail_count); \
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
    G2_elem_proj e_proj, g_proj;
    G2_elem_affine e_affn, g_affn;

    G2_identity_init_proj(&e_proj);
    G2_identity_init_affine(&e_affn);
    G2_generator_init_proj(&g_proj);
    G2_generator_init_affine(&g_affn);

    b1 = G2_is_identity_proj(&e_proj);
    b2 = G2_is_identity_proj(&g_proj);
    b3 = G2_is_identity_affine(&e_affn);
    b4 = G2_is_identity_affine(&g_affn);

    printf("e_proj is identity: %d\n", b1);
    printf("g_proj is identity: %d\n", b2);
    printf("e_affn is identity: %d\n", b3);
    printf("g_affn is identity: %d\n", b4);

    G2_elem_free_proj(&e_proj);
    G2_elem_free_proj(&g_proj);
    G2_elem_free_affine(&e_affn);
    G2_elem_free_affine(&g_affn);

    return b1 && !b2 && b3 && !b4;
}

bool test_is_on_curve()
{
    printf("Running %s...\n", __func__);
    G2_elem_proj e_proj, g_proj, not_on_curve;
    G2_elem_affine e_affn, g_affn;

    G2_identity_init_proj(&e_proj);
    G2_identity_init_affine(&e_affn);
    G2_generator_init_proj(&g_proj);
    G2_generator_init_affine(&g_affn);
    G2_elem_proj_from_str(&not_on_curve,
        "0x0719ed2d3f59c903d00180853a54ca52f536e48a548eaa9ff02a64643a66be9a1476725ac1fa89831072ba6cedbe2101",
        "0x088b6c373a5e0c3a87455eadae0a619f705f2543ffd2110642af7c12363e28005649e6cacdcd3ad30d14583fec46f004",
        "0x0e9893cec412d3fb347b43df6397e4031fabd1e4ff929b0d7d98e2773a2eae42dd75b309621da48342a409805637f5a9",
        "0x0ec97568ce78c1c0eb58a8b19884c4fc99fef923ae6e6b87db03eb63e8b914ea1d04fc2e91d126a1d248c586f6748fec",
        "0x0719ed2d3f59c903d00180853a54ca52f536e48a548eaa9ff02a64643a66be9a1476725ac1fa89831072ba6cedbe2101",
        "0x088b6c373a5e0c3a87455eadae0a619f705f2543ffd2110642af7c12363e28005649e6cacdcd3ad30d14583fec46f004"
    );

    assert(G2_is_on_curve_proj(&e_proj));
    assert(G2_is_on_curve_proj(&g_proj));
    assert(G2_is_on_curve_affine(&e_affn));
    assert(G2_is_on_curve_affine(&g_affn));
    assert(!G2_is_on_curve_proj(&not_on_curve));

    G2_elem_free_proj(&e_proj);
    G2_elem_free_proj(&g_proj);
    G2_elem_free_proj(&not_on_curve);
    G2_elem_free_affine(&e_affn);
    G2_elem_free_affine(&g_affn);

    return true;
}

bool test_affine2proj()
{
    printf("Running %s...\n", __func__);
    G2_elem_affine e_affn, g_affn;
    G2_elem_proj e_proj, g_proj, tmp;

    G2_identity_init_affine(&e_affn);
    G2_generator_init_affine(&g_affn);

    G2_elem_proj_from_str(&e_proj, "0x0", "0x0", "0x0", "0x0", "0x0", "0x0");
    G2_elem_proj_from_str(&g_proj, "0x0", "0x0", "0x0", "0x0", "0x0", "0x0");

    G2_affine2proj(&e_proj, &e_affn);
    G2_affine2proj(&g_proj, &g_affn);

    assert(G2_is_on_curve_proj(&e_proj));
    assert(G2_is_on_curve_proj(&g_proj));

    assert(G2_is_identity_proj(&e_proj));

    G2_generator_init_proj(&tmp);
    assert(G2_equiv_proj(&tmp, &g_proj));

    G2_elem_free_proj(&e_proj);
    G2_elem_free_proj(&g_proj);
    G2_elem_free_proj(&tmp);
    G2_elem_free_affine(&e_affn);
    G2_elem_free_affine(&g_affn);

    return true;
}

bool test_proj2affine()
{
    printf("Running %s...\n", __func__);
    G2_elem_affine e_affn, g_affn, tmp;
    G2_elem_proj e_proj, g_proj;

    G2_identity_init_proj(&e_proj);
    G2_generator_init_proj(&g_proj);

    G2_elem_affine_from_str(&e_affn, "0x0", "0x0", "0x0", "0x0");
    G2_elem_affine_from_str(&g_affn, "0x0", "0x0", "0x0", "0x0");

    G2_proj2affine(&e_affn, &e_proj);
    G2_proj2affine(&g_affn, &g_proj);

    assert(G2_is_on_curve_affine(&e_affn));
    assert(G2_is_identity_affine(&e_affn));

    G2_generator_init_affine(&tmp);

    assert(G2_equiv_affine(&tmp, &g_affn));
    assert(G2_is_on_curve_affine(&g_affn));

    G2_elem_free_proj(&e_proj);
    G2_elem_free_proj(&g_proj);
    G2_elem_free_affine(&tmp);
    G2_elem_free_affine(&e_affn);
    G2_elem_free_affine(&g_affn);

    return true;
}

bool test_add_proj()
{
    printf("Running %s...\n", __func__);
    bool ret;
    G2_elem_proj P, Q, actual, expect;

    G2_elem_proj_from_str(&P,
        "0x000ba06412053c785514ff947c1f6c97723af6b2bd0731a69319e242fe005a5d15a4790a868fe48d7d291fd859fae5c1",
        "0x06701931c92f4fd5386dfbfcdc7f85a308fd6eb79cd18611f1ce14325192241f205649bd44653853b658f425e151e969",
        "0x10c4a709de555d14d708b116ab5d6a049ce271d00f3ac5ec923b5e94caab5f36164d8af5a4926c69c899dd17fc45fa50",
        "0x11dcd8003a25b76785bf91fbfdbd57dd2e14ec61f656ba268b40df9a4a0b0ee588ae0fcf40bd0dd4fbb7346f71a07551",
        "0x0f5ea3e124e08402b14c61f0924ae07931b2d9a2ed59a6121914de24d755f322bb1d02fa161c50b1dfcfbeaa75c9685d",
        "0x1674c0c46a310a000ebcf4c540e970def59770ceb51e9ff18dd00ac5916aa6525707783eb202c1e1ab554d14c814a2cb"
    );

    G2_elem_proj_from_str(&Q,
        "0x0d0621aa460598a2ae95dc8773963bbef6937d711ee9ad39b894d74cf4eb2839affef27184a609b212871cb44b2b5768",
        "0x19b3389e245a6a0debaf8fcab87c9b6a752a861651901954d0c438c682f6a3133d6b7dc7643e0328af8649e21cc5abfd",
        "0x065ae9215806e8a55fd2d9ec4af9d2d448599cdb85d9080b2c9b4766434c33d103730c92c30a69d0602a8804c2a7c65f",
        "0x0e9c342d8a6d4b3a1cbd02c7bdc0e0aa304de41a04569ae33184419e66bbc0271c361c973962955ba6405f0e51beb98b",
        "0x14701735f9267b3559664f0bb6170c92f088db17a9af5a3d7fff70bf69e7fcc224e705e2f271e23dde0bc5c5af268cc8",
        "0x064a76d550996df129652183ad6335e04c940dfc143713673a7586fba1facf0c24536b224dd6b2963335595bb55c2865"
    );

    G2_elem_proj_from_str(&expect,
        "0x0dc5f6ed4c30e8521fa82c3e484dc3a574025a3a88fce59c436171b7474f22975cd692a5f96310ed495055e5a2554bd3",
        "0x13ef18e991288589d5c5525fd5a7d6402dd47aa35494fcb88ad8ee8f635d0e8e92900063e44d6ca9a91ca2f2e268c657",
        "0x1069704ae34cdf2ec371e2afd74ac9ea81f31e30ce30116e29e7624a358fbbc2ed7c8f8a17ae017e87c50a80998c78f5",
        "0x15b6984d5b51ff19de8b0efe5d2e4128eaf6b390eabbb4f34ad99a918892592cccf8769fbff829447eacbc6ed79f5d13",
        "0x109c002863f129aeb23503e9a14c38ed7373865a5083deb8fc0b0ea5a98496c5a05f0cac929220cdff5407d238e87ed9",
        "0x0aaa022d1361871ac6a9f396ee0b1f82c18fedbde63adfede38c1c0e44b3381621c8ee2e8370adb97ee67a17253d4e64"
    );

    G2_identity_init_proj(&actual);

    G2_add_proj(&actual, &actual, &actual);
    assert(G2_is_on_curve_proj(&actual));
    assert(G2_is_identity_proj(&actual));

    G2_add_proj(&actual, &P, &actual);
    assert(G2_is_on_curve_proj(&actual));
    assert(!G2_is_identity_proj(&actual));

    ret = G2_equiv_proj(&actual, &P);
    assert(ret);

    G2_add_proj(&actual, &P, &Q);

    ret = G2_equiv_proj(&actual, &expect);
    printf("Actual:\n");
    printf("\tx0: %s\n", mpz_get_str(NULL, 16, actual.x->a));
    printf("\tx1: %s\n", mpz_get_str(NULL, 16, actual.x->b));
    printf("\ty0: %s\n", mpz_get_str(NULL, 16, actual.y->a));
    printf("\ty1: %s\n", mpz_get_str(NULL, 16, actual.y->b));
    printf("\tz0: %s\n", mpz_get_str(NULL, 16, actual.z->a));
    printf("\tz1: %s\n", mpz_get_str(NULL, 16, actual.z->b));
    printf("Expected:\n");
    printf("\tx0: %s\n", mpz_get_str(NULL, 16, expect.x->a));
    printf("\tx1: %s\n", mpz_get_str(NULL, 16, expect.x->b));
    printf("\ty0: %s\n", mpz_get_str(NULL, 16, expect.y->a));
    printf("\ty1: %s\n", mpz_get_str(NULL, 16, expect.y->b));
    printf("\tz0: %s\n", mpz_get_str(NULL, 16, expect.z->a));
    printf("\tz1: %s\n", mpz_get_str(NULL, 16, expect.z->b));

    G2_elem_free_proj(&P);
    G2_elem_free_proj(&actual);
    G2_elem_free_proj(&expect);

    return ret;
}

bool test_add_mixed()
{
    printf("Running %s...\n", __func__);
    bool ret;
    G2_elem_affine Q;
    G2_elem_proj P, actual, expect, expect2;

    G2_elem_proj_from_str(&P,
        "0x000ba06412053c785514ff947c1f6c97723af6b2bd0731a69319e242fe005a5d15a4790a868fe48d7d291fd859fae5c1",
        "0x06701931c92f4fd5386dfbfcdc7f85a308fd6eb79cd18611f1ce14325192241f205649bd44653853b658f425e151e969",
        "0x10c4a709de555d14d708b116ab5d6a049ce271d00f3ac5ec923b5e94caab5f36164d8af5a4926c69c899dd17fc45fa50",
        "0x11dcd8003a25b76785bf91fbfdbd57dd2e14ec61f656ba268b40df9a4a0b0ee588ae0fcf40bd0dd4fbb7346f71a07551",
        "0x0f5ea3e124e08402b14c61f0924ae07931b2d9a2ed59a6121914de24d755f322bb1d02fa161c50b1dfcfbeaa75c9685d",
        "0x1674c0c46a310a000ebcf4c540e970def59770ceb51e9ff18dd00ac5916aa6525707783eb202c1e1ab554d14c814a2cb"
    );

    G2_elem_affine_from_str(&Q,
        "0x1638533957d540a9d2370f17cc7ed5863bc0b995b8825e0ee1ea1e1e4d00dbae81f14b0bf3611b78c952aacab827a053",
        "0x0a4edef9c1ed7f729f520e47730a124fd70662a904ba1074728114d1031e1572c6c886f6b57ec72a6178288c47c33577",
        "0x0468fb440d82b0630aeb8dca2b5256789a66da69bf91009cbfe6bd221e47aa8ae88dece9764bf3bd999d95d71e4c9899",
        "0x0f6d4552fa65dd2638b361543f887136a43253d9c66c411697003f7a13c308f5422e1aa0a59c8967acdefd8b6e36ccf3"
    );

    G2_elem_proj_from_str(&expect,
        "0x16f99a87ec7a6c50ad0e97f46b02f336f85233faf9d402cd9c2c9e455039d613ff6c46a9a8d98251c72720b19e747bf7",
        "0x0159d2a02719b207bf6fa788f66555f5619adf1f472c9363c3cbd8f3f67dc90a2e45ec921fdabec259abc61677902d0a",
        "0x13c94fffea9d27e1628a3d509d949b94b154d65c32b438d5ebb2a68b5603308027e1fe2df0997417938275536b8a2db2",
        "0x171fa4f9061efd0cefff327c7cc758344d5351b68a104c50eab057c544e0a8397d547ec9a5c261dece19b60a46893206",
        "0x050426f6df8716c8c9527beb74ec6d0106e5af6cc9ec11a1daede3ab32cebe379bb41333c4491d95cbc966afb1354667",
        "0x106b46bfe550c8f2d5ea956500ee47063e4b698614debda2d66d287365c6a136ad23b44127c7f23e575aeb384c2833e6"
    );

    // Make sure that the same point, in a different representation also counts as equivalent.
    G2_elem_proj_from_str(&expect2,
        "0x0dc5f6ed4c30e8521fa82c3e484dc3a574025a3a88fce59c436171b7474f22975cd692a5f96310ed495055e5a2554bd3",
        "0x13ef18e991288589d5c5525fd5a7d6402dd47aa35494fcb88ad8ee8f635d0e8e92900063e44d6ca9a91ca2f2e268c657",
        "0x1069704ae34cdf2ec371e2afd74ac9ea81f31e30ce30116e29e7624a358fbbc2ed7c8f8a17ae017e87c50a80998c78f5",
        "0x15b6984d5b51ff19de8b0efe5d2e4128eaf6b390eabbb4f34ad99a918892592cccf8769fbff829447eacbc6ed79f5d13",
        "0x109c002863f129aeb23503e9a14c38ed7373865a5083deb8fc0b0ea5a98496c5a05f0cac929220cdff5407d238e87ed9",
        "0x0aaa022d1361871ac6a9f396ee0b1f82c18fedbde63adfede38c1c0e44b3381621c8ee2e8370adb97ee67a17253d4e64"
    );

    G2_identity_init_proj(&actual);

    assert(G2_is_on_curve_proj(&P));
    assert(!G2_is_identity_proj(&P));
    assert(G2_is_on_curve_affine(&Q));
    assert(!G2_is_identity_affine(&Q));

    G2_add_mixed(&actual, &P, &Q);

    assert(G2_is_on_curve_proj(&actual));
    assert(!G2_is_identity_proj(&actual));

    ret = G2_equiv_proj(&actual, &expect);
    printf("Actual:\n");
    printf("\tx0: %s\n", mpz_get_str(NULL, 16, actual.x->a));
    printf("\tx1: %s\n", mpz_get_str(NULL, 16, actual.x->b));
    printf("\ty0: %s\n", mpz_get_str(NULL, 16, actual.y->a));
    printf("\ty1: %s\n", mpz_get_str(NULL, 16, actual.y->b));
    printf("\tz0: %s\n", mpz_get_str(NULL, 16, actual.z->a));
    printf("\tz1: %s\n", mpz_get_str(NULL, 16, actual.z->b));
    printf("Expected:\n");
    printf("\tx0: %s\n", mpz_get_str(NULL, 16, expect.x->a));
    printf("\tx1: %s\n", mpz_get_str(NULL, 16, expect.x->b));
    printf("\ty0: %s\n", mpz_get_str(NULL, 16, expect.y->a));
    printf("\ty1: %s\n", mpz_get_str(NULL, 16, expect.y->b));
    printf("\tz0: %s\n", mpz_get_str(NULL, 16, expect.z->a));
    printf("\tz1: %s\n", mpz_get_str(NULL, 16, expect.z->b));

    ret = ret && G2_equiv_proj(&actual, &expect2);
    printf("Actual:\n");
    printf("\tx0: %s\n", mpz_get_str(NULL, 16, actual.x->a));
    printf("\tx1: %s\n", mpz_get_str(NULL, 16, actual.x->b));
    printf("\ty0: %s\n", mpz_get_str(NULL, 16, actual.y->a));
    printf("\ty1: %s\n", mpz_get_str(NULL, 16, actual.y->b));
    printf("\tz0: %s\n", mpz_get_str(NULL, 16, actual.z->a));
    printf("\tz1: %s\n", mpz_get_str(NULL, 16, actual.z->b));
    printf("Expected (not normalized):\n");
    printf("\tx0: %s\n", mpz_get_str(NULL, 16, expect2.x->a));
    printf("\tx1: %s\n", mpz_get_str(NULL, 16, expect2.x->b));
    printf("\ty0: %s\n", mpz_get_str(NULL, 16, expect2.y->a));
    printf("\ty1: %s\n", mpz_get_str(NULL, 16, expect2.y->b));
    printf("\tz0: %s\n", mpz_get_str(NULL, 16, expect2.z->a));
    printf("\tz1: %s\n", mpz_get_str(NULL, 16, expect2.z->b));

    G2_elem_free_proj(&P);
    G2_elem_free_affine(&Q);
    G2_elem_free_proj(&actual);
    G2_elem_free_proj(&expect);
    G2_elem_free_proj(&expect2);

    return ret;
}

bool test_double_proj()
{
    printf("Running %s...\n", __func__);
    bool ret;
    G2_elem_proj P, actual, expect;

    G2_elem_proj_from_str(&P,
        "0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8",
        "0x13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e",
        "0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801",
        "0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be",
        "0x1",
        "0x0"
    );

    G2_elem_proj_from_str(&expect,
        "0x0d0621aa460598a2ae95dc8773963bbef6937d711ee9ad39b894d74cf4eb2839affef27184a609b212871cb44b2b5768",
        "0x19b3389e245a6a0debaf8fcab87c9b6a752a861651901954d0c438c682f6a3133d6b7dc7643e0328af8649e21cc5abfd",
        "0x065ae9215806e8a55fd2d9ec4af9d2d448599cdb85d9080b2c9b4766434c33d103730c92c30a69d0602a8804c2a7c65f",
        "0x0e9c342d8a6d4b3a1cbd02c7bdc0e0aa304de41a04569ae33184419e66bbc0271c361c973962955ba6405f0e51beb98b",
        "0x14701735f9267b3559664f0bb6170c92f088db17a9af5a3d7fff70bf69e7fcc224e705e2f271e23dde0bc5c5af268cc8",
        "0x064a76d550996df129652183ad6335e04c940dfc143713673a7586fba1facf0c24536b224dd6b2963335595bb55c2865"
    );

    G2_identity_init_proj(&actual);

    G2_double_proj(&actual, &actual);
    assert(G2_is_on_curve_proj(&actual));
    assert(G2_is_identity_proj(&actual));

    G2_double_proj(&actual, &P);

    assert(G2_is_on_curve_proj(&actual));
    assert(!G2_is_identity_proj(&actual));

    ret = G2_equiv_proj(&actual, &expect);
    printf("Actual:\n");
    printf("\tx0: %s\n", mpz_get_str(NULL, 16, actual.x->a));
    printf("\tx1: %s\n", mpz_get_str(NULL, 16, actual.x->b));
    printf("\ty0: %s\n", mpz_get_str(NULL, 16, actual.y->a));
    printf("\ty1: %s\n", mpz_get_str(NULL, 16, actual.y->b));
    printf("\tz0: %s\n", mpz_get_str(NULL, 16, actual.z->a));
    printf("\tz1: %s\n", mpz_get_str(NULL, 16, actual.z->b));
    printf("Expected:\n");
    printf("\tx0: %s\n", mpz_get_str(NULL, 16, expect.x->a));
    printf("\tx1: %s\n", mpz_get_str(NULL, 16, expect.x->b));
    printf("\ty0: %s\n", mpz_get_str(NULL, 16, expect.y->a));
    printf("\ty1: %s\n", mpz_get_str(NULL, 16, expect.y->b));
    printf("\tz0: %s\n", mpz_get_str(NULL, 16, expect.z->a));
    printf("\tz1: %s\n", mpz_get_str(NULL, 16, expect.z->b));

    G2_elem_free_proj(&P);
    G2_elem_free_proj(&actual);
    G2_elem_free_proj(&expect);

    return ret;
}

bool test_mul_scalar()
{
    printf("Running %s...\n", __func__);
    bool ret;
    G2_elem_affine P, actual, expect;
    mpz_t m;

    mpz_init_set_str(m, "0x0c9af8fcfd18dc6a102260d25e1de7918505fac80ef519d22b568298a56da71b", 0);
    G2_elem_affine_from_str(&expect,
        "0x08459e7c2ee2c69f369f2851e477e0119e34b01b9eadd8c58ca8410d7396e093a33a369c36b125ffa380be87c1619030",
        "0x08640e1bb61022ade73b70d8ab178bc37ae7ee600b5e1668238dbccc5e980f7f4c8f3af04a152f506d2110ca43ca6997",
        "0x0403053e2bdf3968be0da9204e920ee5ca93e943ee17eb6c6bd1da894594b051e706fed86602b55264f0971ef7f883e8",
        "0x0450c250f7ed0dc53a726238437126336d9b3b2a137f170da17c4a25535039993d1d610614ba4b0accc56e6bc18f6da7"
    );

    G2_generator_init_affine(&P);
    G2_identity_init_affine(&actual);

    G2_mul_scalar(&actual, &P, m);

    ret = G2_equiv_affine(&actual, &expect);
    printf("Actual:\n");
    printf("\tx0: %s\n", mpz_get_str(NULL, 16, actual.x->a));
    printf("\tx1: %s\n", mpz_get_str(NULL, 16, actual.x->b));
    printf("\ty0: %s\n", mpz_get_str(NULL, 16, actual.y->a));
    printf("\ty1: %s\n", mpz_get_str(NULL, 16, actual.y->b));
    printf("Expected:\n");
    printf("\tx0: %s\n", mpz_get_str(NULL, 16, expect.x->a));
    printf("\tx1: %s\n", mpz_get_str(NULL, 16, expect.x->b));
    printf("\ty0: %s\n", mpz_get_str(NULL, 16, expect.y->a));
    printf("\ty1: %s\n", mpz_get_str(NULL, 16, expect.y->b));

    G2_elem_free_affine(&P);
    G2_elem_free_affine(&actual);
    G2_elem_free_affine(&expect);
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
