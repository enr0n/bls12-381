#include <stdio.h>

#include "finite_field.h"

#define FP2_DEFINE_TEST(op, expect_a_str, expect_b_str)   \
    int res;                                              \
    fp2_elem expect, actual;                              \
                                                          \
    fp2_elem_init(&actual);                               \
                                                          \
    fp2_elem_from_str(                                    \
        &expect,                                          \
        expect_a_str,                                     \
        expect_b_str                                      \
    );                                                    \
                                                          \
    op;                                                   \
                                                          \
    res = fp2_equal(&actual, &expect);                    \
    printf("Actual:\n");                                  \
    printf("\ta: %s\n", mpz_get_str(NULL, 16, actual.a)); \
    printf("\tb: %s\n", mpz_get_str(NULL, 16, actual.b)); \
    printf("Expected:\n");                                \
    printf("\ta: %s\n", mpz_get_str(NULL, 16, expect.a)); \
    printf("\tb: %s\n", mpz_get_str(NULL, 16, expect.b)); \
    fp2_elem_free(&actual);                              \
    fp2_elem_free(&expect);                              \
                                                          \
    return res;

void test_elems_init(fp2_elem *e1, fp2_elem *e2)
{
    fp2_elem_from_str(
        e1,
        "0x05676605b164404c1f9a3973b47e86eed77e3d22fde82372eb797db6a0d3e198f79fd3bfb4e508d39b9dd53d757d9e0b",
        "0x076c4f72966ea34a1339330a64b25e90133b12340aacd558050f5b7657507c98d72d6c6a7ffed2a6d2c3b37af2255f1f"
    );

    fp2_elem_from_str(
        e2,
        "0x02ab3cf3b87533a113526d3c74cc69fc5995c4925eb653770cc2c6599c849ebe54a84d3e24b2cd2ac11f49a1d114e02e",
        "0x070b77d9a4c794c74a10412939f981cd9e7c1f77056120b99e2856332bff06a96e549b6a2f18e9361a6060800e1d9dad"
    );
}

int test_fp2_add(const fp2_elem *e1, const fp2_elem *e2)
{
    FP2_DEFINE_TEST(
        fp2_add(&actual, e1, e2),
        "0x0812a2f969d973ed32eca6b0294af0eb311401b55c9e76e9f83c44103d5880574c4820fdd997d5fe5cbd1edf46927e39",
        "0x0e77c74c3b3638115d4974339eabe05db1b731ab100df611a337b1a9834f8342458207d4af17bbdced2413fb0042fccc"
    );
}

int test_fp2_mul(const fp2_elem *e1, const fp2_elem *e2)
{
    FP2_DEFINE_TEST(
        fp2_mul(&actual, e1, e2),
        "0x0cd7fcd385025ea22f9d4c8aa420419eac3092d7dfe30b296ced11e0b92c38c16434d0b18073a1081927106895de01c6",
        "0x0fedd0b9a7c48b1eb1e6a82914d957b6b2535dec296908035a1f11d316ba4b04deead89d6477f30caade4aa90fc6097d"
    );
}

int test_fp2_square(const fp2_elem *e1)
{
    FP2_DEFINE_TEST(
        fp2_square(&actual, e1),
        "0x0dded60ef884ac4e16a180aeccabc1c2ff8488210444c1fb3499574e110ac627403d567932b6c72197530cb1956cc190",
        "0x0d307452085efda84f6e3b72aac08b412d82b42bbe3cd5a0ccfa5d6cc3d731a7c727737b5634eb282ec5554c08e3effa"
    );
}

int test_fp2_inv(const fp2_elem *e1)
{
    FP2_DEFINE_TEST(
        fp2_inv(&actual, e1),
        "0x06d454959dcd2001b9b0eaabbdf8129470f641dfe0a8827c96e558fba25038c6cc18e806b9a2e15a9e3c9ff079dd778f",
        "0x099a20ff4b2a19767f903dfb8718f991c4ff5091c6551a278f4fe8df60ab0118de863bc5aa844c9170c1e96d216deff4"
    );
}

int main()
{
    int result, pass_count, fail_count;

    pass_count = fail_count = 0;

    fp_params_init();

    fp2_elem e1, e2;
    test_elems_init(&e1, &e2);

    printf("Running test_fp2_add...\n");
    result = test_fp2_add(&e1, &e2);
    if (!result) {
        printf("FAIL\n\n");
        fail_count++;
    } else {
        printf("PASS\n\n");
        pass_count++;
    }

    printf("Running test_fp2_mul...\n");
    result = test_fp2_mul(&e1, &e2);
    if (!result) {
        printf("FAIL\n\n");
        fail_count++;
    } else {
        printf("PASS\n\n");
        pass_count++;
    }

    printf("Running test_fp2_square...\n");
    result = test_fp2_square(&e1);
    if (!result) {
        printf("FAIL\n\n");
        fail_count++;
    } else {
        printf("PASS\n\n");
        pass_count++;
    }

    printf("Running test_fp2_inv...\n");
    result = test_fp2_inv(&e1);
    if (!result) {
        printf("FAIL\n\n");
        fail_count++;
    } else {
        printf("PASS\n\n");
        pass_count++;
    }

    fp2_elem_free(&e1);
    fp2_elem_free(&e2);

    fp_params_free();

    printf("FP2 test results: %d passed, %d failed\n", pass_count, fail_count);

    return !!(fail_count);
}
