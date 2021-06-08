#include <stdio.h>
#include <stdlib.h>

#include "params.h"
#include "finite_field.h"

#define FP6_DEFINE_TEST(op,                                   \
                        expect_a0_str, expect_a1_str,         \
                        expect_b0_str, expect_b1_str,         \
                        expect_c0_str, expect_c1_str)         \
    int res;                                                  \
    fp6_elem expect, actual;                                  \
                                                              \
    fp6_elem_init(&actual);                                   \
                                                              \
    fp6_elem_from_str(                                        \
        &expect,                                              \
        expect_a0_str, expect_a1_str,                         \
        expect_b0_str, expect_b1_str,                         \
        expect_c0_str, expect_c1_str                          \
    );                                                        \
                                                              \
    op;                                                       \
                                                              \
    res = fp6_equal(&actual, &expect);                        \
    printf("Actual:    \n");                                  \
    printf("\ta0: %s\n", mpz_get_str(NULL, 16, actual.a->a)); \
    printf("\ta1: %s\n", mpz_get_str(NULL, 16, actual.a->b)); \
    printf("\tb0: %s\n", mpz_get_str(NULL, 16, actual.b->a)); \
    printf("\tb1: %s\n", mpz_get_str(NULL, 16, actual.b->b)); \
    printf("\tc0: %s\n", mpz_get_str(NULL, 16, actual.c->a)); \
    printf("\tc1: %s\n", mpz_get_str(NULL, 16, actual.c->b)); \
    printf("Expected:\n");                                    \
    printf("\ta0: %s\n", mpz_get_str(NULL, 16, expect.a->a)); \
    printf("\ta1: %s\n", mpz_get_str(NULL, 16, expect.a->b)); \
    printf("\tb0: %s\n", mpz_get_str(NULL, 16, expect.b->a)); \
    printf("\tb1: %s\n", mpz_get_str(NULL, 16, expect.b->b)); \
    printf("\tc0: %s\n", mpz_get_str(NULL, 16, expect.c->a)); \
    printf("\tc1: %s\n", mpz_get_str(NULL, 16, expect.c->b)); \
    fp6_elem_clear(&actual);                                  \
    fp6_elem_clear(&expect);                                  \
                                                              \
    return res;

void test_elems_init(fp6_elem *e1, fp6_elem *e2)
{
    fp6_elem_from_str(
        e1,
        "0x8c0ed57c",
        "0x8c0ed563",
        "0x8c0ed562",
        "0x8c0ed561",
        "0x8c0ed560",
        "0x8c0ed567"
    );

    fp6_elem_from_str(
        e2,
        "0x8c0ed566",
        "0x8c0ed565",
        "0x8c0ed564",
        "0x8c0ed56b",
        "0x8c0ed56a",
        "0x8c0ed569"
    );
}

int test_fp6_add(const fp6_elem *e1, const fp6_elem *e2, const mpz_t p)
{
    FP6_DEFINE_TEST(
        fp6_add(&actual, e1, e2, p),
        "0x0000000000000000000000000000000000000000000000000000000000000000000000000000000000000001181daae2",
        "0x0000000000000000000000000000000000000000000000000000000000000000000000000000000000000001181daac8",
        "0x0000000000000000000000000000000000000000000000000000000000000000000000000000000000000001181daac6",
        "0x0000000000000000000000000000000000000000000000000000000000000000000000000000000000000001181daacc",
        "0x0000000000000000000000000000000000000000000000000000000000000000000000000000000000000001181daaca",
        "0x0000000000000000000000000000000000000000000000000000000000000000000000000000000000000001181daad0"
    );
}

int test_fp6_mul(const fp6_elem *e1, const fp6_elem *e2, const mpz_t p)
{
    FP6_DEFINE_TEST(
        fp6_mul(&actual, e1, e2, p),
        "0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153fffe877e16fb74df198a",
        "0x00000000000000000000000000000000000000000000000000000000000000000000000000000001cbc15d96ae5f0654",
        "0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffff20be8b7f5e9be1fc",
        "0x00000000000000000000000000000000000000000000000000000000000000000000000000000001cbc15d99f6b80757",
        "0x0000000000000000000000000000000000000000000000000000000000000000000000000000000000000007a8cfac17",
        "0x00000000000000000000000000000000000000000000000000000000000000000000000000000001cbc15d9d3f11079e"
    );
}

int test_fp6_square(const fp6_elem *e1,  const mpz_t p)
{
    FP6_DEFINE_TEST(
        fp6_square(&actual, e1, p),
        "0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153fffe877e1715b7a71e48",
        "0x00000000000000000000000000000000000000000000000000000000000000000000000000000001cbc15d947e23b0f6",
        "0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffff20be8b9310b1e3e4",
        "0x00000000000000000000000000000000000000000000000000000000000000000000000000000001cbc15d9366060593",
        "0x0000000000000000000000000000000000000000000000000000000000000000000000000000000000000014ca33ac19",
        "0x00000000000000000000000000000000000000000000000000000000000000000000000000000001cbc15d9b0ed5b24c"
    );
}

int test_fp6_inv(const fp6_elem *e1,  const mpz_t p)
{
    FP6_DEFINE_TEST(
        fp6_inv(&actual, e1, p),
        "0x0f5d3888ce8fb7e81051c8ab459ec57c5cac1cc5bad8d497f5d941f454303adb51e1f13856f01eb3b6e96a349bfd506c",
        "0x00b5066943316367b6f802b3e2feb164112cf652c56a8fe4395f22843335c3f6db7d125937cfa97b25f59ab0a0ccae5d",
        "0x1308087806edee690d656f71d22ae367dfc4963d04e46a2853d2e3d80c1b7494544e1cba6032f0db129086d6da4a1d61",
        "0x18b68ea30cef69b18b39c24d25fff7c55b0ac9f2255e96f66d8aac1249925ad6c5f810e86dacdf69633e1b646ce82c64",
        "0x07b1fa73b4b5c51f32ea6b76a0fa01e75c758e06983c5f6813890ef66fa9285f69c75be4c28fb6fb742f8f4cb2e91545",
        "0x10fc74a9af1ceb669cf38cd6277e1cd4f7f4c10b46219aed14de3948b4b64a318f174c3b5e403449dcbc2aab5f712c57"
    );
}

int main()
{
    int result, pass_count, fail_count;

    pass_count = fail_count = 0;

    mpz_t p;
    mpz_init_set_str(p, BLS12_381_MODULUS, 0);

    fp6_elem e1, e2;
    test_elems_init(&e1, &e2);

    printf("Running test_fp6_add...\n");
    result = test_fp6_add(&e1, &e2, p);
    if (!result) {
        printf("FAIL\n\n");
        fail_count++;
    } else {
        printf("PASS\n\n");
        pass_count++;
    }

    printf("Running test_fp6_mul...\n");
    result = test_fp6_mul(&e1, &e2, p);
    if (!result) {
        printf("FAIL\n\n");
        fail_count++;
    } else {
        printf("PASS\n\n");
        pass_count++;
    }

    printf("Running test_fp6_square...\n");
    result = test_fp6_square(&e1, p);
    if (!result) {
        printf("FAIL\n\n");
        fail_count++;
    } else {
        printf("PASS\n\n");
        pass_count++;
    }

    printf("Running test_fp6_inv...\n");
    result = test_fp6_inv(&e1, p);
    if (!result) {
        printf("FAIL\n\n");
        fail_count++;
    } else {
        printf("PASS\n\n");
        pass_count++;
    }

    fp6_elem_clear(&e1);
    fp6_elem_clear(&e2);
    mpz_clear(p);

    printf("FP6 test results: %d passed, %d failed\n", pass_count, fail_count);

    return !!(fail_count);
}
