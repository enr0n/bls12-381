#include <stdio.h>
#include <stdlib.h>

#include "finite_field.h"
#include "G1.h"
#include "G2.h"
#include "pairing.h"

int test_final_exponentiation()
{
    int res;
    fp12_elem expect, actual, f;

    const char *fa[6] = {
        "0x05194f5785436c8debf0eb2bab4c6ef3de7dc0633c85769173777b782bf897fa45025fd03e7be941123c4ee19910e62e",
        "0x0e760e96f911ae38a6042da82d7b0e30787864e725e9d5462d224c91c4497104d838d566d894564bc19e09d8af706c3f",
        "0x0751a051e0beb4a0e2351a7527d813b371e189056307d718a446e4016a3df787568a842f3401768dc03b966bd1db90ac",
        "0x0cbc592a19a3f60c9938676b257b9c01ed9d708f9428b29e272a811d13d734485970d9d3f1c097b12bfa3d1678096b1d",
        "0x159ef660e2d84185f55c0ccae1dd7f8f71b12c0beb7a431fede9e62794d9154e9a0ce4715f64b032492459076224c99b",
        "0x087d1320fe5bad5c2d8e12c49e6aff41a0b80e1497bbe85682e22ed853f256041bdf97ef02bdb5d80a5f9bc31d85f25e"
    };
    const char *fb[6] = {
        "0x05d928cb508feeb3329e51aa0bec4f33ba865a22da5a4e97eb31b78c0150c0c6134f0f94bd0154b28430ee4c6052e82b",
        "0x159bfbbdc31bb5cb0082c59e5f744773335ef1fdddb8ed86a1c23f61f18800b647ff7dae335fb9ab5fcf2188cb64d72d",
        "0x1431225e128c5e2bfafb9eba23746150907688583f52e07fcde4cc93452b0c2bcd0f0893b48a696c403c6980d0940741",
        "0x07508024863ec263bded120e45deb29c1f1303a056b279e116cb5fdb03013db19f81e78fa2b2b409cb2ce8e3ba96f4e6",
        "0x1868172fbbeb861d69c6c10f315c273d08312812c643dbf60588d0de3d2c4b3e9b21acd402f7ddee53f1c4797646ba96",
        "0x1562633d4f2387ff79a0f625a6989072296a946ca6bbfa3fef879defde15ed96d205b2eebb454f48fb76fa8a845bcba7"
    };
    fp12_elem_from_str(&f, fa, fb);

    const char *ra[6] = {
        "0x1250ebd871fc0a92a7b2d83168d0d727272d441befa15c503dd8e90ce98db3e7b6d194f60839c508a84305aaca1789b6",
        "0x089a1c5b46e5110b86750ec6a532348868a84045483c92b7af5af689452eafabf1a8943e50439f1d59882a98eaa0170f",
        "0x1368bb445c7c2d209703f239689ce34c0378a68e72a6b3b216da0e22a5031b54ddff57309396b38c881c4c849ec23e87",
        "0x193502b86edb8857c273fa075a50512937e0794e1e65a7617c90d8bd66065b1fffe51d7a579973b1315021ec3c19934f",
        "0x01b2f522473d171391125ba84dc4007cfbf2f8da752f7c74185203fcca589ac719c34dffbbaad8431dad1c1fb597aaa5",
        "0x018107154f25a764bd3c79937a45b84546da634b8f6be14a8061e55cceba478b23f7dacaa35c8ca78beae9624045b4b6"
    };
    const char *rb[6] = {
        "0x19f26337d205fb469cd6bd15c3d5a04dc88784fbb3d0b2dbdea54d43b2b73f2cbb12d58386a8703e0f948226e47ee89d",
        "0x06fba23eb7c5af0d9f80940ca771b6ffd5857baaf222eb95a7d2809d61bfe02e1bfd1b68ff02f0b8102ae1c2d5d5ab1a",
        "0x11b8b424cd48bf38fcef68083b0b0ec5c81a93b330ee1a677d0d15ff7b984e8978ef48881e32fac91b93b47333e2ba57",
        "0x03350f55a7aefcd3c31b4fcb6ce5771cc6a0e9786ab5973320c806ad360829107ba810c5a09ffdd9be2291a0c25a99a2",
        "0x04c581234d086a9902249b64728ffd21a189e87935a954051c7cdba7b3872629a4fafc05066245cb9108f0242d0fe3ef",
        "0x0f41e58663bf08cf068672cbd01a7ec73baca4d72ca93544deff686bfd6df543d48eaa24afe47e1efde449383b676631"
    };
    fp12_elem_from_str(&expect, ra, rb);
    fp12_elem_init(&actual);

    final_exponentiation(&actual, &f);

    res = fp12_equal(&actual, &expect);
    printf("Actual:    \n");
    printf("\ta0: %s\n", mpz_get_str(NULL, 16, actual.a->a->a));
    printf("\ta1: %s\n", mpz_get_str(NULL, 16, actual.a->a->b));
    printf("\ta2: %s\n", mpz_get_str(NULL, 16, actual.a->b->a));
    printf("\ta3: %s\n", mpz_get_str(NULL, 16, actual.a->b->b));
    printf("\ta4: %s\n", mpz_get_str(NULL, 16, actual.a->c->a));
    printf("\ta5: %s\n", mpz_get_str(NULL, 16, actual.a->c->b));
    printf("\tb0: %s\n", mpz_get_str(NULL, 16, actual.b->a->a));
    printf("\tb1: %s\n", mpz_get_str(NULL, 16, actual.b->a->b));
    printf("\tb2: %s\n", mpz_get_str(NULL, 16, actual.b->b->a));
    printf("\tb3: %s\n", mpz_get_str(NULL, 16, actual.b->b->b));
    printf("\tb4: %s\n", mpz_get_str(NULL, 16, actual.b->c->a));
    printf("\tb5: %s\n", mpz_get_str(NULL, 16, actual.b->c->b));
    printf("Expected:\n");
    printf("\ta0: %s\n", mpz_get_str(NULL, 16, expect.a->a->a));
    printf("\ta1: %s\n", mpz_get_str(NULL, 16, expect.a->a->b));
    printf("\ta2: %s\n", mpz_get_str(NULL, 16, expect.a->b->a));
    printf("\ta3: %s\n", mpz_get_str(NULL, 16, expect.a->b->b));
    printf("\ta4: %s\n", mpz_get_str(NULL, 16, expect.a->c->a));
    printf("\ta5: %s\n", mpz_get_str(NULL, 16, expect.a->c->b));
    printf("\tb0: %s\n", mpz_get_str(NULL, 16, expect.b->a->a));
    printf("\tb1: %s\n", mpz_get_str(NULL, 16, expect.b->a->b));
    printf("\tb2: %s\n", mpz_get_str(NULL, 16, expect.b->b->a));
    printf("\tb3: %s\n", mpz_get_str(NULL, 16, expect.b->b->b));
    printf("\tb4: %s\n", mpz_get_str(NULL, 16, expect.b->c->a));
    printf("\tb5: %s\n", mpz_get_str(NULL, 16, expect.b->c->b));
    fp12_elem_clear(&actual);
    fp12_elem_clear(&expect);

    return res;
}

int test_miller_loop_point_acc()
{
    int res;
    G1_elem_affine P;
    G2_elem_affine Q, actual, expect;
    fp12_elem tmp;

    fp12_elem_init(&tmp);

    G1_generator_init_affine(&P);
    G2_generator_init_affine(&Q);
    G2_identity_init_affine(&actual);

    G2_elem_affine_from_str(&expect,
        "0x149ee6d25a1c8648c86f8673946935a6c41c801f2eb23fd9171fd41f31b7ec5c4a85f5b285f02508d928c716f9c1bc67",
        "0x06e63a1561fa8259ba6f6e15257a5cd4ab79697e45abf5042aa7cf28b80aa054b18d535b71cf70dafd2b694040d19479",
        "0x0fcf61319e51e41eed89a6b0ff12b4da6f664b92422160a52083df7086215506947457947cd8917791ac308465a849c6",
        "0x14c164377443f5accb1eb0dfbb4910e98a49e53d8ef016c76df99c2649d3f49b1a3d8b6a55e5881f8ea0ecfdc3de9b2f"
    );

    miller_loop(&tmp, &actual, &P, &Q);

    res = G2_equiv_affine(&actual, &expect);
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

    fp12_elem_clear(&tmp);

    G1_elem_free_affine(&P);
    G2_elem_free_affine(&Q);
    G2_elem_free_affine(&actual);
    G2_elem_free_affine(&expect);

    return res;
}

int test_ate()
{
    int res;
    fp12_elem actual, expect;
    G1_elem_affine P;
    G2_elem_affine Q;

    const char *ra[6] = {
        "0x1250ebd871fc0a92a7b2d83168d0d727272d441befa15c503dd8e90ce98db3e7b6d194f60839c508a84305aaca1789b6",
        "0x089a1c5b46e5110b86750ec6a532348868a84045483c92b7af5af689452eafabf1a8943e50439f1d59882a98eaa0170f",
        "0x1368bb445c7c2d209703f239689ce34c0378a68e72a6b3b216da0e22a5031b54ddff57309396b38c881c4c849ec23e87",
        "0x193502b86edb8857c273fa075a50512937e0794e1e65a7617c90d8bd66065b1fffe51d7a579973b1315021ec3c19934f",
        "0x01b2f522473d171391125ba84dc4007cfbf2f8da752f7c74185203fcca589ac719c34dffbbaad8431dad1c1fb597aaa5",
        "0x018107154f25a764bd3c79937a45b84546da634b8f6be14a8061e55cceba478b23f7dacaa35c8ca78beae9624045b4b6"
    };
    const char *rb[6] = {
        "0x19f26337d205fb469cd6bd15c3d5a04dc88784fbb3d0b2dbdea54d43b2b73f2cbb12d58386a8703e0f948226e47ee89d",
        "0x06fba23eb7c5af0d9f80940ca771b6ffd5857baaf222eb95a7d2809d61bfe02e1bfd1b68ff02f0b8102ae1c2d5d5ab1a",
        "0x11b8b424cd48bf38fcef68083b0b0ec5c81a93b330ee1a677d0d15ff7b984e8978ef48881e32fac91b93b47333e2ba57",
        "0x03350f55a7aefcd3c31b4fcb6ce5771cc6a0e9786ab5973320c806ad360829107ba810c5a09ffdd9be2291a0c25a99a2",
        "0x04c581234d086a9902249b64728ffd21a189e87935a954051c7cdba7b3872629a4fafc05066245cb9108f0242d0fe3ef",
        "0x0f41e58663bf08cf068672cbd01a7ec73baca4d72ca93544deff686bfd6df543d48eaa24afe47e1efde449383b676631"
    };
    fp12_elem_from_str(&expect, ra, rb);
    fp12_elem_init(&actual);

    G1_generator_init_affine(&P);
    G2_generator_init_affine(&Q);

    ate(&actual, &Q, &P);

    res = fp12_equal(&actual, &expect);
    printf("Actual:    \n");
    printf("\ta0: %s\n", mpz_get_str(NULL, 16, actual.a->a->a));
    printf("\ta1: %s\n", mpz_get_str(NULL, 16, actual.a->a->b));
    printf("\ta2: %s\n", mpz_get_str(NULL, 16, actual.a->b->a));
    printf("\ta3: %s\n", mpz_get_str(NULL, 16, actual.a->b->b));
    printf("\ta4: %s\n", mpz_get_str(NULL, 16, actual.a->c->a));
    printf("\ta5: %s\n", mpz_get_str(NULL, 16, actual.a->c->b));
    printf("\tb0: %s\n", mpz_get_str(NULL, 16, actual.b->a->a));
    printf("\tb1: %s\n", mpz_get_str(NULL, 16, actual.b->a->b));
    printf("\tb2: %s\n", mpz_get_str(NULL, 16, actual.b->b->a));
    printf("\tb3: %s\n", mpz_get_str(NULL, 16, actual.b->b->b));
    printf("\tb4: %s\n", mpz_get_str(NULL, 16, actual.b->c->a));
    printf("\tb5: %s\n", mpz_get_str(NULL, 16, actual.b->c->b));
    printf("Expected:\n");
    printf("\ta0: %s\n", mpz_get_str(NULL, 16, expect.a->a->a));
    printf("\ta1: %s\n", mpz_get_str(NULL, 16, expect.a->a->b));
    printf("\ta2: %s\n", mpz_get_str(NULL, 16, expect.a->b->a));
    printf("\ta3: %s\n", mpz_get_str(NULL, 16, expect.a->b->b));
    printf("\ta4: %s\n", mpz_get_str(NULL, 16, expect.a->c->a));
    printf("\ta5: %s\n", mpz_get_str(NULL, 16, expect.a->c->b));
    printf("\tb0: %s\n", mpz_get_str(NULL, 16, expect.b->a->a));
    printf("\tb1: %s\n", mpz_get_str(NULL, 16, expect.b->a->b));
    printf("\tb2: %s\n", mpz_get_str(NULL, 16, expect.b->b->a));
    printf("\tb3: %s\n", mpz_get_str(NULL, 16, expect.b->b->b));
    printf("\tb4: %s\n", mpz_get_str(NULL, 16, expect.b->c->a));
    printf("\tb5: %s\n", mpz_get_str(NULL, 16, expect.b->c->b));
    fp12_elem_clear(&actual);
    fp12_elem_clear(&expect);

    G1_elem_free_affine(&P);
    G2_elem_free_affine(&Q);

    return res;
}

int main()
{
    int result, pass_count, fail_count;

    pass_count = fail_count = 0;

    BLS12_381_init();

    printf("Running test_final_exponentiation...\n");
    result = test_final_exponentiation();
    if (!result) {
        printf("FAIL\n\n");
        fail_count++;
    } else {
        printf("PASS\n\n");
        pass_count++;
    }

    printf("Running test_miller_loop_point_acc...\n");
    result = test_miller_loop_point_acc();
    if (!result) {
        printf("FAIL\n\n");
        fail_count++;
    } else {
        printf("PASS\n\n");
        pass_count++;
    }

    printf("Running test_ate...\n");
    result = test_ate();
    if (!result) {
        printf("FAIL\n\n");
        fail_count++;
    } else {
        printf("PASS\n\n");
        pass_count++;
    }

    BLS12_381_free();

    printf("Pairing test results: %d passed, %d failed\n", pass_count, fail_count);

    return !!(fail_count);
}
