#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "finite_field.h"

#define FP12_DEFINE_TEST(op,                                     \
                        expect_a0_str, expect_a1_str,            \
                        expect_a2_str, expect_a3_str,            \
                        expect_a4_str, expect_a5_str,            \
                        expect_b0_str, expect_b1_str,            \
                        expect_b2_str, expect_b3_str,            \
                        expect_b4_str, expect_b5_str)            \
    int res;                                                     \
    fp12_elem expect, actual;                                    \
                                                                 \
    fp12_elem_init(&actual);                                     \
                                                                 \
    const char *expect_a[6] = {                                  \
          expect_a0_str, expect_a1_str,                          \
          expect_a2_str, expect_a3_str,                          \
          expect_a4_str, expect_a5_str                           \
    };                                                           \
                                                                 \
    const char *expect_b[6] = {                                  \
          expect_b0_str, expect_b1_str,                          \
          expect_b2_str, expect_b3_str,                          \
          expect_b4_str, expect_b5_str                           \
    };                                                           \
                                                                 \
    fp12_elem_from_str(&expect, expect_a, expect_b);             \
                                                                 \
    op;                                                          \
                                                                 \
    res = fp12_equal(&actual, &expect);                          \
    printf("Actual:    \n");                                     \
    printf("\ta0: %s\n", mpz_get_str(NULL, 16, actual.a->a->a)); \
    printf("\ta1: %s\n", mpz_get_str(NULL, 16, actual.a->a->b)); \
    printf("\ta2: %s\n", mpz_get_str(NULL, 16, actual.a->b->a)); \
    printf("\ta3: %s\n", mpz_get_str(NULL, 16, actual.a->b->b)); \
    printf("\ta4: %s\n", mpz_get_str(NULL, 16, actual.a->c->a)); \
    printf("\ta5: %s\n", mpz_get_str(NULL, 16, actual.a->c->b)); \
    printf("\tb0: %s\n", mpz_get_str(NULL, 16, actual.b->a->a)); \
    printf("\tb1: %s\n", mpz_get_str(NULL, 16, actual.b->a->b)); \
    printf("\tb2: %s\n", mpz_get_str(NULL, 16, actual.b->b->a)); \
    printf("\tb3: %s\n", mpz_get_str(NULL, 16, actual.b->b->b)); \
    printf("\tb4: %s\n", mpz_get_str(NULL, 16, actual.b->c->a)); \
    printf("\tb5: %s\n", mpz_get_str(NULL, 16, actual.b->c->b)); \
    printf("Expected:\n");                                       \
    printf("\ta0: %s\n", mpz_get_str(NULL, 16, expect.a->a->a)); \
    printf("\ta1: %s\n", mpz_get_str(NULL, 16, expect.a->a->b)); \
    printf("\ta2: %s\n", mpz_get_str(NULL, 16, expect.a->b->a)); \
    printf("\ta3: %s\n", mpz_get_str(NULL, 16, expect.a->b->b)); \
    printf("\ta4: %s\n", mpz_get_str(NULL, 16, expect.a->c->a)); \
    printf("\ta5: %s\n", mpz_get_str(NULL, 16, expect.a->c->b)); \
    printf("\tb0: %s\n", mpz_get_str(NULL, 16, expect.b->a->a)); \
    printf("\tb1: %s\n", mpz_get_str(NULL, 16, expect.b->a->b)); \
    printf("\tb2: %s\n", mpz_get_str(NULL, 16, expect.b->b->a)); \
    printf("\tb3: %s\n", mpz_get_str(NULL, 16, expect.b->b->b)); \
    printf("\tb4: %s\n", mpz_get_str(NULL, 16, expect.b->c->a)); \
    printf("\tb5: %s\n", mpz_get_str(NULL, 16, expect.b->c->b)); \
    fp12_elem_free(&actual);                                    \
    fp12_elem_free(&expect);                                    \
                                                                 \
    return res;

void test_elems_init(fp12_elem *e1, fp12_elem *e2)
{
    const char *e1a[6] = {
        "0x15f30bb146420a3596196bc00af6c7091b1a912d3fe7024390015b157ce58b8cd3b853fc617f50270c4dafd0b8364404",
        "0x154e7e0da1249a829901d1abbda7c91f95378f62afdff255ada13ece005b9f200223cf7241ef7224541007e5e06f8f56",
        "0x1554a7b734fd6a38a44245f5a4802051ce3cd4e43ef0d2245c67706ce8bcadce52b7f555b942eda2c5f0ec0e52955a3f",
        "0x129a9681d8972bf82d2d412892e1965d884373226ac49a18b0e0429a6724493bbb7cb552488b0916b200d2cb39ed87b3",
        "0x00329b0efda5da0fe6f68ae49534a765298d16166d1d4a2a6431d5152ea73c01f4099634fc82edc01744cb997f49fa44",
        "0x06549512e88b6a95eb31e043e5edba7e9e64a488b140186ca782becc4cc4c16cf3190036543bb8a8403bbdfe4b0f6f5c"
    };

    const char *e1b[6] = {
        "0x13481c6fffc20eb036e011b81b1d3ce16fbc76865078fef1a6866d1f31f2a2a5a102d326b37dedfdc0c5196258fb405f",
        "0x1967190a72709dff96d913851bfc0cd2ad3ff4f84e2e65f835ff59469e8a8c234ea6287f925088b30c9f9c61d9dc7d05",
        "0x0192a396c9f3dd698610e6663ae63a148303a67bb2fe23b93a2541af555eef94f4968726c731c9dd1b4609c2e6c498e2",
        "0x0911cb58b939e97b259e92cbc3f762d0a4392863a9625facf0b70b09ff27c6e4e5a029aaf840f3d4a86e7ef3834db63b",
        "0x0e73f51654c73f1674f30fe5a450bf7881f4e6e272ae2c41e573fabf2f6717dea83f8a7112db1dbe03b8173e57c57383",
        "0x04cc1157baea01e9d4da36d60eb9fe96d8308a9f3f4a38256cf6b460fd3b3152dad05f740a8fb550d9ad14acf6b61d57"
    };
    fp12_elem_from_str(e1, e1a, e1b);

    const char *e2a[6] = {
        "0x08eb56c0a9dc196203f6493a107198ad1ea4b6559c1c760d4a31a6d0cafb0a31534c9aa6df8994e07f4877b4cc90fcf6",
        "0x03a4235753d501c0a8a80c20ef7d6742f8945f128f7682a66dc146e7a3452e057a2268566c4015f00cd841cd916aebad",
        "0x1536e376b042ac9c91c708ce541ec0febfc442931ea206ea47be83f27dcca7c6aa944b7e7391e3d509c3442686bd7f82",
        "0x149bb326456b03e2739079be25d52a021b25710fd1df6357c3e82b7545fbe220a26dc2d30c8a03af75a874158c7a358d",
        "0x14f3fdf39f3aab18d657d23acfa6ea19b8f25a0d932bdfd037f18c3263ece00934630916f93c7c220ce5489d32992449",
        "0x0ab75a4fb139b593d0f95d1141b8f91ee1da9c02de1ecb866b3f2e0c6fe2dad46648c117c7a071ab58b1680002f4181b"
    };

    const char *e2b[6] = {
        "0x10037f4b9d2667c806ebe566021be42ed3069e5461c99d5074fae7a6d32581da10865513e4dceab3f65bb6104b82723d",
        "0x00e5f5a8675f33fee5cc6665a0091c499b4840d88d71ccdbc7f0aced12965be1916ed600e431059f47fa8b524cc21128",
        "0x185b3bd20f55b08ade36986e680e429fab32c96066d0e157e56f585662db45243589bfc60f0cb828ed80a58a93806dbb",
        "0x18028b15d1e9786464831e24047a36fe0f80812c3030f2e18b018bf04d2749b8e77381789662c9505a0094e7424026f0",
        "0x1005b726c0f3eb3b2783bfbfb9f49edfad4be32094adf427da33f4dd377d24755063164f70d4b9ef75b9cdcfe485ea70",
        "0x0cd6ccb109f7cc80c01fcd00a4428a81ec32b65ac96cfb2896a8698de0588ec05bd264c05f46d27f8e45c56f6506beae"
    };
    fp12_elem_from_str(e2, e2a, e2b);
}

int test_fp12_add(const fp12_elem *e1, const fp12_elem *e2)
{
    FP12_DEFINE_TEST(
        fp12_add(&actual, e1, e2),
        "0x04dd5087b69e3cfd4ef40d43d81cb2ded547fbfde87e659173022f45512f9f9a0858eea48fb4e507d197278584c7964f",
        "0x18f2a164f4f99c4341a9ddccad2530628dcbee753f5674fc1b6285b5a3a0cd257c4637c8ae2f881460e849b371da7b03",
        "0x108a7943abc0303aeaeda70db55334792989cbf26a0dc64f3cf521be6fd85f70dea040d57b80d17815b53034d9532f16",
        "0x0d3537bde482494055a21330756b13883ef198ad491eeab10d979b6eb66f35383f3e7826a3c10cc66daa46e0c6681295",
        "0x152699029ce08528bd4e5d1f64db917ee27f7024004929fa9c23614792941c0b286c9f4bf5bf69e2242a1436b1e31e8d",
        "0x110bef6299c52029bc2b3d5527a6b39d803f408b8f5ee3f312c1ecd8bca79c415961c14e1bdc2a5398ed25fe4e038777",
        "0x094a89d163688fddf2b04f67d9ed7438de4bc955bebd8982b45082250e672e5b92dd283be706d8b1fd21cf72a47e07f1",
        "0x004bfcc8a04feb643189d23478b97c44e410ea4be81b201496bf3392ba6ff1e0c168fe81c52d8e529a9b27b4269ee382",
        "0x19eddf68d9498df464477ed4a2f47cb42e366fdc19cf05111f949a05b83a34b92a2046ecd63e820608c6af4d7a45069d",
        "0x0713448451a37b453f0609398525ecf74f425e0ae60e3fcf1487c459559e1a79ae67ab24dd4fbd25487013dac58e3280",
        "0x04789a52dc3b43b7515b27ef1af9b180cac97e7e13d70daa58771cfb7033462fd9f6a0c1d25bd7adbf72e50e3c4bb348",
        "0x11a2de08c4e1ce6a94fa03d6b2fc8918c46340fa08b7334e039f1deedd93c01336a2c43469d687d067f2da1c5bbcdc05"
    );
}

int test_fp12_mul(const fp12_elem *e1, const fp12_elem *e2)
{
    FP12_DEFINE_TEST(
        fp12_mul(&actual, e1, e2),
        "0x0129afc7c1ca7d1f962c83320b1f917416f67d39c5555d6d0fe36fee9ed7268b644b855036b05d016ef6d3a246026046",
        "0x0d6e3fba17f36cfe1f88ee187892b650827fb693538213ed785551294534f7771d7438b13a289031ea7b2f75b64c2a27",
        "0x006364f3ed0b46136ff623364595b668196a032814bf422cc3ebec42d8ed984ca1a9a1a0ad14a537241d0dfcb0c1dd95",
        "0x02680bd3eb89aa87261c6d6b16a9a0b157d65994699f9109f035ed5c33c7e3cf74339bb002b09c81638ede38bda4d1b1",
        "0x07879ada3a5032e92bd6fd23c9aca9af8d1cace444b8c0cead456442f2486f952dd69959a25bfe4a9d7ed1617fc85e35",
        "0x14d9b49d750e3234fb4076bc9cd84563dd67c7e02935248385ad7c5a8179042f0a1bfee3538ce6f66529c6a06d6bdfb5",
        "0x07d324bb6fb13e2058f0f3dea6fd9b627e4d29ad84d0bc2813c45049273f48b60a7ae85ae05ae77d38e7a62e3d152d17",
        "0x00343f072a0d670d454ff64f58d53cdc7a65b2a6f4b0d32652f256e7857a41f26dc76fc2b868ae193254d51119b12553",
        "0x092a1658eb9b3d3c0bc9b2479f0a9f22259b1fa599f0eae47094f04788e3b6b7adcb86b1f127c9e555d5a14348b1e843",
        "0x19305585246b802f64de7a7ad24dc21c93caf38f95595e3247a61da36c017d1fed2840254bcf80b3e9df9dadbdc2be01",
        "0x0ce8082c150f325da559b64d4bf6cdecf57b4cd08c9ba8904100fb926ef0f1197a0b5d7278752486a6317f63318bd6c4",
        "0x11261ebccab20e9c9228b218c35a4bb1e5abbfc01535987518f0e76af842e6b99461b6f8cbf95e8698c063809d81b15d"
    );
}

int test_fp12_square(const fp12_elem *e1)
{
    FP12_DEFINE_TEST(
        fp12_square(&actual, e1),
        "0x190868e8fac47fda5c2a8131bffeac4a1762f21089ce303f0f48601472c36dcab0484a3ba2c744321f95cb9f6a233e0b",
        "0x17654a404c677dea53ed0be5ff5bbdd785d950de47759c532b00951eccdc0e59fff2114e85a3ffbb8bc6c7a4b08efd1b",
        "0x03e0993de3f59a473df3ffbcf296799bb6432f7f07b54857deb383b8214710c8f22c05a503079318004c74d9e44cdcc6",
        "0x08893437da4595deb892c894262c505344391dbfe4ffaab56bb7ea831e2172bf7be94f7ecb8062c56ba4b3c0b23fefac",
        "0x0309a24597b14dd0de296a188e7aa14c772a860f3bb4a871313b74999af95d3dd52cb141fc2859af6fe68a902986dfaf",
        "0x1118aea76bfe11f72620aab9c43c3d9c4dbdef9b82f5afaa63d0d39fdcba35ea1d3127eb014f978b3c187ce7226785bb",
        "0x00d3d7547d23ed8c972e8e833fbd7c83ec8abb1c11a662e3415be8fc54d142822ded146d6ee0c0fe38cc089ed51ad410",
        "0x0b8c4e0618b6ba68dc721dda927b69f6935df0efd59c4da17964129728bb4e1a1f1ada7e0cc161e3895540117f6a20ca",
        "0x176fdcaa4311d056920a3fc7bbe2a1dab3766f1c5e8c0dda28990e63443f2e646815c997c3c9ab13ec22c1d5e91e931c",
        "0x09fb834104bc84a0f1361af5cee3f88b91370e2548df45fd668529cd30ec04d0c237a0e495eed432fceb6447bd2344e0",
        "0x08bccbdc6b2621760a67a65adf69b386cb07fadf02fd262f61f58583a8ac9c8a08000d056c1d82cf420205f5b4184e23",
        "0x11ece3887f4e66254df1e1e48bf1b8bb491cbddb32d4a846c006e334c83abc7fd8345c6412f07c1dd7f0e1c62365e483"
    );
}

int test_fp12_inv(const fp12_elem *e1)
{
    FP12_DEFINE_TEST(
        fp12_inv(&actual, e1),
        "0x0cede514d1f8a3751e7ce444093ffa1d35dd2e30310806b407ba1948247ab9e52cfaf49915575d1e75e944a00c0ca2bd",
        "0x0efa0737debdcad772f304546151a73e923c2613be65112290d6fef6af7089b4dc2073224faac826eebc9384ef55abe8",
        "0x0931f46c846f1599e1f77e11d4af739c96ad5d57b0a276071dc2e998a7bae3413cb7dbc67c25e9152f095786a857e8c6",
        "0x052ee2047680480938e2b69a5aa24d5ed7ec76e8d7c3a70242e3bc1ee0bfffdebc339b22ae1f3f79c9e602fd97994234",
        "0x00dadaaf8101f5cb9876260230437046ee8a070ca7aead923260504b7ab40a75f37abde788ae19e92712ba357cca90e0",
        "0x19104f41e277c40b1776202d24d01c1625593aba6a07e567bcb0050eb7a069844acb64c80b72cd15b749266ca862c77a",
        "0x16008a2f3d3c5125ac2768b6f8f67393f271f04af22fe58d3b8a011eba3c141bb4810dc3dc5c485a9770974c6e079ccb",
        "0x12307506f2af0a464313d0b1b58a89fad5ed0afbcd7b5309644cfb01969185658bf9a5cbe7bd5708fb673fadbbb3a40d",
        "0x00622d3ac69696f367a7e5b5e4baa290891aa5f44a1a763a15fbea3e4f74e052586d8bf81d106bb3e1ec9db6172cdef6",
        "0x0d7a59b37e2e353ef6bb04dbd56841f04b50189996550602d93df7e5c6b511c6735031ba43337e069ff7df8bf03ffa67",
        "0x13e736ea0183ac3af18929935f3e65e41de0f6056d77f1c2f193a1883a3bb405d602acc3083638ee3d3e2afdf74404dc",
        "0x0cc556f48855d00fb42260b84756fcedd9eeabbbdd600550b016c0c20052afbbba1f4678e61a608a42b74f7b93e5978f"
    );
}

int test_fp12_frobenius(const fp12_elem *e1)
{
    FP12_DEFINE_TEST(
        fp12_frobenius(&actual, e1),
        "0x15f30bb146420a3596196bc00af6c7091b1a912d3fe7024390015b157ce58b8cd3b853fc617f50270c4dafd0b8364404",
        "0x04b293dc985b4c17b219d60a85a3e3b7cf3fbc2243a52069b98f93d2f65557041c88308c6f648ddb65eef81a1f901b55",
        "0x02481c0453dde3b97f0da08d9b51b759443d35c9f851ba0c0a3790b5cb83cd9338da21ccad32815731e8bfcce5c67e93",
        "0x09750a21a5bf579a8547a4962d93d17d5b4524a6875b7a0520d7075dd6e3351ecb10f65d1ab9769515a17ccaa80a2709",
        "0x12be8c09e110109b513c0f93e6063cffc191bcb0f78cf5e7f935c49a8b64be6574ca83735303807f65a52e7b45d66c6e",
        "0x10b50724c658a333c47eedfb0248c118da9a7b250c571972e92eee26a4828c20d3fe94f583d13bc3b3fb53a2e240f8e1",
        "0x0d5e4cde27d75477c3c3b859eaf60bb2523808bf53ee4d149721aedd46e061d533856fdcdb8028ae3aae7fc1c148ec7d",
        "0x0d4607d984afaebd55039de595953cc03eaf171f453dc04f54f134e849358c5e6f91970fe1ac1ffdefc32d06d92c89d6",
        "0x1565f1dc8dcf6165429bf9f7267864fb4dd752d030f253c03404308d1112c5f3992d9fde39f7b7b25342694ae50388a8",
        "0x160405304c584fdfa0963f6b74e102f3afc7a14889af7679dbac658fb8f6fb854076ba6ac787182ede7c55d20994fa80",
        "0x0d3bd9e1ec2dc95390cc5a49fb7664c7c08ca96a34acf8cdb07c9bb0b1d5e375df619a5c1637ba70291a29a57f3dba60",
        "0x1320fc638c14abd4f7256d8d690d3ce4e6bf14221e8719ff7c1dce5ae7f98d5671a69d7003ce6757ab6e7cdafa4ddd30"
    );
}

int test_fp12_pow(const fp12_elem *e1)
{
    int res;
    fp12_elem expect, actual;
    mpz_t exp;

    fp12_elem_init(&actual);
    fp12_elem_init(&expect);
    mpz_init(exp);

    // expect = e1^6
    fp12_square(&expect, e1);
    fp12_mul(&expect, &expect, e1);
    fp12_square(&expect, &expect);

    mpz_set_ui(exp, 6);
    fp12_pow(&actual, e1, exp);

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

    res = fp12_equal(&actual, &expect);

    // expect = e1^-6
    fp12_inv(&expect, &expect);

    mpz_set_si(exp, -6);
    fp12_pow(&actual, e1, exp);

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

    res = res && fp12_equal(&actual, &expect);

    // expect = e1^5
    fp12_square(&expect, e1);
    fp12_square(&expect, &expect);
    fp12_mul(&expect, &expect, e1);

    mpz_set_ui(exp, 5);
    fp12_pow(&actual, e1, exp);

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

    res = fp12_equal(&actual, &expect);

    // expect = e1^-5
    fp12_inv(&expect, &expect);

    mpz_set_si(exp, -5);
    fp12_pow(&actual, e1, exp);

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

    res = res && fp12_equal(&actual, &expect);

    fp12_elem_free(&actual);
    fp12_elem_free(&expect);
    mpz_clear(exp);

    return res;
}

int main()
{
    int result, pass_count, fail_count;

    pass_count = fail_count = 0;

    fp_params_init();

    fp12_elem e1, e2;
    test_elems_init(&e1, &e2);

    printf("Running test_fp12_add...\n");
    result = test_fp12_add(&e1, &e2);
    if (!result) {
        printf("FAIL\n\n");
        fail_count++;
    } else {
        printf("PASS\n\n");
        pass_count++;
    }

    printf("Running test_fp12_mul...\n");
    result = test_fp12_mul(&e1, &e2);
    if (!result) {
        printf("FAIL\n\n");
        fail_count++;
    } else {
        printf("PASS\n\n");
        pass_count++;
    }

    printf("Running test_fp12_square...\n");
    result = test_fp12_square(&e1);
    if (!result) {
        printf("FAIL\n\n");
        fail_count++;
    } else {
        printf("PASS\n\n");
        pass_count++;
    }

    printf("Running test_fp12_inv...\n");
    result = test_fp12_inv(&e1);
    if (!result) {
        printf("FAIL\n\n");
        fail_count++;
    } else {
        printf("PASS\n\n");
        pass_count++;
    }

    printf("Running test_fp12_frobenius...\n");
    result = test_fp12_frobenius(&e1);
    if (!result) {
        printf("FAIL\n\n");
        fail_count++;
    } else {
        printf("PASS\n\n");
        pass_count++;
    }

    printf("Running test_fp12_pow...\n");
    result = test_fp12_pow(&e1);
    if (!result) {
        printf("FAIL\n\n");
        fail_count++;
    } else {
        printf("PASS\n\n");
        pass_count++;
    }

    fp12_elem_free(&e1);
    fp12_elem_free(&e2);

    fp_params_free();

    printf("FP12 test results: %d passed, %d failed\n", pass_count, fail_count);

    return !!(fail_count);
}
