#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "finite_field.h"

#include "hash.h"

#include "BLS12_381.h"

#define TEST_MAIN_INIT                  \
    int result, pass_count, fail_count; \
                                        \
    pass_count = fail_count = 0;

#define TEST_MAIN_RETURN                                                         \
    printf("Hash test results: %d passed, %d failed\n", pass_count, fail_count); \
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

/* Test vectors from https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-11#appendix-K.1 */
#define EXPAND_MESSAGE_XMD_TEST_DST "QUUX-V01-CS02-with-expander"
#define G1_TEST_DST "QUUX-V01-CS02-with-BLS12381G1_XMD:SHA-256_SSWU_RO_"
#define G2_TEST_DST "QUUX-V01-CS02-with-BLS12381G2_XMD:SHA-256_SSWU_RO_"

const char *test_vec_msg[10] = {
    "",
    "abc",
    "abcdef0123456789",
    "q128_qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq",
    "a512_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
    "",
    "abc",
    "abcdef0123456789",
    "q128_qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq",
    "a512_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
};

const uint32_t test_vec_len_in_bytes[10] = {
    0x20,
    0x20,
    0x20,
    0x20,
    0x20,
    0x80,
    0x80,
    0x80,
    0x80,
    0x80
};

const char *test_vec_uniform_bytes[10] = {
    "f659819a6473c1835b25ea59e3d38914c98b374f0970b7e4c92181df928fca88",
    "1c38f7c211ef233367b2420d04798fa4698080a8901021a795a1151775fe4da7",
    "8f7e7b66791f0da0dbb5ec7c22ec637f79758c0a48170bfb7c4611bd304ece89",
    "72d5aa5ec810370d1f0013c0df2f1d65699494ee2a39f72e1716b1b964e1c642",
    "3b8e704fc48336aca4c2a12195b720882f2162a4b7b13a9c350db46f429b771b",
    "8bcffd1a3cae24cf9cd7ab85628fd111bb17e3739d3b53f89580d217aa79526f1708354a76a402d3569d6a9d19ef3de4d0b991e4f54b9f20dcde9b95a66824cbdf6c1a963a1913d43fd7ac443a02fc5d9d8d77e2071b86ab114a9f34150954a7531da568a1ea8c760861c0cde2005afc2c114042ee7b5848f5303f0611cf297f",
    "fe994ec51bdaa821598047b3121c149b364b178606d5e72bfbb713933acc29c186f316baecf7ea22212f2496ef3f785a27e84a40d8b299cec56032763eceeff4c61bd1fe65ed81decafff4a31d0198619c0aa0c6c51fca15520789925e813dcfd318b542f8799441271f4db9ee3b8092a7a2e8d5b75b73e28fb1ab6b4573c192",
    "c9ec7941811b1e19ce98e21db28d22259354d4d0643e301175e2f474e030d32694e9dd5520dde93f3600d8edad94e5c364903088a7228cc9eff685d7eaac50d5a5a8229d083b51de4ccc3733917f4b9535a819b445814890b7029b5de805bf62b33a4dc7e24acdf2c924e9fe50d55a6b832c8c84c7f82474b34e48c6d43867be",
    "48e256ddba722053ba462b2b93351fc966026e6d6db493189798181c5f3feea377b5a6f1d8368d7453faef715f9aecb078cd402cbd548c0e179c4ed1e4c7e5b048e0a39d31817b5b24f50db58bb3720fe96ba53db947842120a068816ac05c159bb5266c63658b4f000cbf87b1209a225def8ef1dca917bcda79a1e42acd8069",
    "396962db47f749ec3b5042ce2452b619607f27fd3939ece2746a7614fb83a1d097f554df3927b084e55de92c7871430d6b95c2a13896d8a33bc48587b1f66d21b128a1a8240d5b0c26dfe795a1a842a0807bb148b77c2ef82ed4b6c9f7fcb732e7f94466c8b51e52bf378fba044a31f5cb44583a892f5969dcd73b3fa128816e"
};

const char *test_vec_G1_field_elems_0[5] = {
    "0x0ba14bd907ad64a016293ee7c2d276b8eae71f25a4b941eece7b0d89f17f75cb3ae5438a614fb61d6835ad59f29c564f",
    "0x0d921c33f2bad966478a03ca35d05719bdf92d347557ea166e5bba579eea9b83e9afa5c088573c2281410369fbd32951",
    "0x062d1865eb80ebfa73dcfc45db1ad4266b9f3a93219976a3790ab8d52d3e5f1e62f3b01795e36834b17b70e7b76246d4",
    "0x010476f6a060453c0b1ad0b628f3e57c23039ee16eea5e71bb87c3b5419b1255dc0e5883322e563b84a29543823c0e86",
    "0x0a8ffa7447f6be1c5a2ea4b959c9454b431e29ccc0802bc052413a9c5b4f9aac67a93431bd480d15be1e057c8a08e8c6"
};

const char *test_vec_G1_field_elems_1[5] = {
    "0x019b9bd7979f12657976de2884c7cce192b82c177c80e0ec604436a7f538d231552f0d96d9f7babe5fa3b19b3ff25ac9",
    "0x003574a00b109ada2f26a37a91f9d1e740dffd8d69ec0c35e1e9f4652c7dba61123e9dd2e76c655d956e2b3462611139",
    "0x0cdc3e2f271f29c4ff75020857ce6c5d36008c9b48385ea2f2bf6f96f428a3deb798aa033cd482d1cdc8b30178b08e3a",
    "0x0b1a912064fb0554b180e07af7e787f1f883a0470759c03c1b6509eb8ce980d1670305ae7b928226bb58fdc0a419f46e",
    "0x05d487032f602c90fa7625dbafe0f4a49ef4a6b0b33d7bb349ff4cf5410d297fd6241876e3e77b651cfc8191e40a68b7"
};

const char *test_vec_G2_field_elems_0_0[5] = {
    "0x03dbc2cce174e91ba93cbb08f26b917f98194a2ea08d1cce75b2b9cc9f21689d80bd79b594a613d0a68eb807dfdc1cf8",
    "0x15f7c0aa8f6b296ab5ff9c2c7581ade64f4ee6f1bf18f55179ff44a2cf355fa53dd2a2158c5ecb17d7c52f63e7195771",
    "0x0313d9325081b415bfd4e5364efaef392ecf69b087496973b229303e1816d2080971470f7da112c4eb43053130b785e1",
    "0x025820cefc7d06fd38de7d8e370e0da8a52498be9b53cba9927b2ef5c6de1e12e12f188bbc7bc923864883c57e49e253",
    "0x190b513da3e66fc9a3587b78c76d1d132b1152174d0b83e3c1114066392579a45824c5fa17649ab89299ddd4bda54935"
};

const char *test_vec_G2_field_elems_0_1[5] = {
    "0x05a2acec64114845711a54199ea339abd125ba38253b70a92c876df10598bd1986b739cad67961eb94f7076511b3b39a",
    "0x01c8067bf4c0ba709aa8b9abc3d1cef589a4758e09ef53732d670fd8739a7274e111ba2fcaa71b3d33df2a3a0c8529dd",
    "0x062f84cb21ed89406890c051a0e8b9cf6c575cf6e8e18ecf63ba86826b0ae02548d83b483b79e48512b82a6c0686df8f",
    "0x034147b77ce337a52e5948f66db0bab47a8d038e712123bb381899b6ab5ad20f02805601e6104c29df18c254b8618c7b",
    "0x12ab625b0fe0ebd1367fe9fac57bb1168891846039b4216b9d94007b674de2d79126870e88aeef54b2ec717a887dcf39"
};

const char *test_vec_G2_field_elems_1_0[5] = {
    "0x02f99798e8a5acdeed60d7e18e9120521ba1f47ec090984662846bc825de191b5b7641148c0dbc237726a334473eee94",
    "0x187111d5e088b6b9acfdfad078c4dacf72dcd17ca17c82be35e79f8c372a693f60a033b461d81b025864a0ad051a06e4",
    "0x1739123845406baa7be5c5dc74492051b6d42504de008c635f3535bb831d478a341420e67dcc7b46b2e8cba5379cca97",
    "0x0930315cae1f9a6017c3f0c8f2314baa130e1cf13f6532bff0a8a1790cd70af918088c3db94bda214e896e1543629795",
    "0x0e6a42010cf435fb5bacc156a585e1ea3294cc81d0ceb81924d95040298380b164f702275892cedd81b62de3aba3f6b5"
};

const char *test_vec_G2_field_elems_1_1[5] = {
    "0x145a81e418d4010cc027a68f14391b30074e89e60ee7a22f87217b2f6eb0c4b94c9115b436e6fa4607e95a98de30a435",
    "0x08b852331c96ed983e497ebc6dee9b75e373d923b729194af8e72a051ea586f3538a6ebb1e80881a082fa2b24df9f566",
    "0x01897665d9cb5db16a27657760bbea7951f67ad68f8d55f7113f24ba6ddd82caef240a9bfa627972279974894701d975",
    "0x10c4df2cacf67ea3cb3108b00d4cbd0b3968031ebc8eac4b1ebcefe84d6b715fde66bef0219951ece29d1facc8a520ef",
    "0x117d9a0defc57a33ed208428cb84e54c85a6840e7648480ae428838989d25d97a0af8e3255be62b25c2a85630d2dddd8"
};

const char *test_vec_G1_mapped_points_0_x[5] = {
    "0x11a3cce7e1d90975990066b2f2643b9540fa40d6137780df4e753a8054d07580db3b7f1f03396333d4a359d1fe3766fe",
    "0x125435adce8e1cbd1c803e7123f45392dc6e326d292499c2c45c5865985fd74fe8f042ecdeeec5ecac80680d04317d80",
    "0x08834484878c217682f6d09a4b51444802fdba3d7f2df9903a0ddadb92130ebbfa807fffa0eabf257d7b48272410afff",
    "0x0cbd7f84ad2c99643fea7a7ac8f52d63d66cefa06d9a56148e58b984b3dd25e1f41ff47154543343949c64f88d48a710",
    "0x0cf97e6dbd0947857f3e578231d07b309c622ade08f2c08b32ff372bd90db19467b2563cc997d4407968d4ac80e154f8"
};

const char *test_vec_G1_mapped_points_0_y[5] = {
    "0x0eeaf6d794e479e270da10fdaf768db4c96b650a74518fc67b04b03927754bac66f3ac720404f339ecdcc028afa091b7",
    "0x0e8828948c989126595ee30e4f7c931cbd6f4570735624fd25aef2fa41d3f79cfb4b4ee7b7e55a8ce013af2a5ba20bf2",
    "0x0b318f7ecf77f45a0f038e62d7098221d2dbbca2a394164e2e3fe953dc714ac2cde412d8f2d7f0c03b259e6795a2508e",
    "0x052c00e4ed52d000d94881a5638ae9274d3efc8bc77bc0e5c650de04a000b2c334a9e80b85282a00f3148dfdface0865",
    "0x127f0cddf2613058101a5701f4cb9d0861fd6c2a1b8e0afe194fccf586a3201a53874a2761a9ab6d7220c68661a35ab3"
};

const char *test_vec_G1_mapped_points_1_x[5] = {
    "0x160003aaf1632b13396dbad518effa00fff532f604de1a7fc2082ff4cb0afa2d63b2c32da1bef2bf6c5ca62dc6b72f9c",
    "0x11def93719829ecda3b46aa8c31fc3ac9c34b428982b898369608e4f042babee6c77ab9218aad5c87ba785481eff8ae4",
    "0x158418ed6b27e2549f05531a8281b5822b31c3bf3144277fbb977f8d6e2694fedceb7011b3c2b192f23e2a44b2bd106e",
    "0x06493fb68f0d513af08be0372f849436a787e7b701ae31cb964d968021d6ba6bd7d26a38aaa5a68e8c21a6b17dc8b579",
    "0x092f1acfa62b05f95884c6791fba989bbe58044ee6355d100973bf9553ade52b47929264e6ae770fb264582d8dce512a"
};

const char *test_vec_G1_mapped_points_1_y[5] = {
    "0x0d8bb2d14e20cf9f6036152ed386d79189415b6d015a20133acb4e019139b94e9c146aaad5817f866c95d609a361735e",
    "0x0007c9cef122ccf2efd233d6eb9bfc680aa276652b0661f4f820a653cec1db7ff69899f8e52b8e92b025a12c822a6ce6",
    "0x1879074f344471fac5f839e2b4920789643c075792bec5af4282c73f7941cda5aa77b00085eb10e206171b9787c4169f",
    "0x02e98f2ccf5802b05ffaac7c20018bc0c0b2fd580216c4aa2275d2909dc0c92d0d0bdc979226adeb57a29933536b6bb4",
    "0x028e6d0169a72cfedb737be45db6c401d3adfb12c58c619c82b93a5dfcccef12290de530b0480575ddc8397cda0bbebf"
};

const char *test_vec_G1_hash_out_x[5] = {
    "0x052926add2207b76ca4fa57a8734416c8dc95e24501772c814278700eed6d1e4e8cf62d9c09db0fac349612b759e79a1",
    "0x03567bc5ef9c690c2ab2ecdf6a96ef1c139cc0b2f284dca0a9a7943388a49a3aee664ba5379a7655d3c68900be2f6903",
    "0x11e0b079dea29a68f0383ee94fed1b940995272407e3bb916bbf268c263ddd57a6a27200a784cbc248e84f357ce82d98",
    "0x15f68eaa693b95ccb85215dc65fa81038d69629f70aeee0d0f677cf22285e7bf58d7cb86eefe8f2e9bc3f8cb84fac488",
    "0x082aabae8b7dedb0e78aeb619ad3bfd9277a2f77ba7fad20ef6aabdc6c31d19ba5a6d12283553294c1825c4b3ca2dcfe"
};

const char *test_vec_G1_hash_out_y[5] = {
    "0x08ba738453bfed09cb546dbb0783dbb3a5f1f566ed67bb6be0e8c67e2e81a4cc68ee29813bb7994998f3eae0c9c6a265",
    "0x0b9c15f3fe6e5cf4211f346271d7b01c8f3b28be689c8429c85b67af215533311f0b8dfaaa154fa6b88176c229f2885d",
    "0x03a87ae2caf14e8ee52e51fa2ed8eefe80f02457004ba4d486d6aa1f517c0889501dc7413753f9599b099ebcbbd2d709",
    "0x1807a1d50c29f430b8cafc4f8638dfeeadf51211e1602a5f184443076715f91bb90a48ba1e370edce6ae1062f5e6dd38",
    "0x05b84ae5a942248eea39e1d91030458c40153f3b654ab7872d779ad1e942856a20c438e8d99bc8abfbf74729ce1f7ac8"
};

bool test_message_expand_xmd(int test_num)
{
    printf("Running %s #%d\n", __func__, test_num);

    bool ret;
    uint32_t len_in_bytes;
    octet_string *bytes;
    char *out;

    len_in_bytes = test_vec_len_in_bytes[test_num];

    octet_string_alloc(&bytes, len_in_bytes);

    expand_message_xmd(bytes, test_vec_msg[test_num],
                       EXPAND_MESSAGE_XMD_TEST_DST, len_in_bytes);

    out = octet_string_to_str(bytes);
    ret = strcmp(out, test_vec_uniform_bytes[test_num]) == 0;

    printf("Actual:\n\t%s\n", out);
    printf("Expected:\n\t%s\n", test_vec_uniform_bytes[test_num]);

    octet_string_free(bytes);
    free(out);

    return ret;
}

bool test_hash_to_field_fp(int test_num)
{
    printf("Running %s #%d\n", __func__, test_num);

    bool ret;
    mpz_t *elems;
    mpz_t expect;

    mpz_init(expect);

    elems = hash_to_field_fp(test_vec_msg[test_num], G1_TEST_DST, 2);

    mpz_set_str(expect, test_vec_G1_field_elems_0[test_num], 0);
    ret = (mpz_cmp(expect, elems[0]) == 0);
    printf("u[0]\nActual:\n\t%s\n", mpz_get_str(NULL, 16, elems[0]));
    printf("Expected:\n\t%s\n", mpz_get_str(NULL, 16, expect));

    mpz_set_str(expect, test_vec_G1_field_elems_1[test_num], 0);
    ret = ret && (mpz_cmp(expect, elems[1]) == 0);
    printf("u[1]\nActual:\n\t%s\n", mpz_get_str(NULL, 16, elems[1]));
    printf("Expected:\n\t%s\n", mpz_get_str(NULL, 16, expect));

    mpz_clear(expect);

    return ret;
}

bool test_hash_to_field_fp2(int test_num)
{
    printf("Running %s #%d\n", __func__, test_num);

    bool ret;
    fp2_elem **elems;
    fp2_elem expect;

    elems = hash_to_field_fp2(test_vec_msg[test_num], G2_TEST_DST, 2);

    fp2_elem_from_str(&expect,
            test_vec_G2_field_elems_0_0[test_num],
            test_vec_G2_field_elems_0_1[test_num]);

    ret = fp2_equal(&expect, elems[0]);
    printf("u[0]\nActual:\n");
    printf("\ta: %s\n", mpz_get_str(NULL, 16, elems[0]->a));
    printf("\tb: %s\n", mpz_get_str(NULL, 16, elems[0]->b));
    printf("Expected:\n");
    printf("\ta: %s\n", mpz_get_str(NULL, 16, expect.a));
    printf("\tb: %s\n", mpz_get_str(NULL, 16, expect.b));
    fp2_elem_free(&expect);

    fp2_elem_from_str(&expect,
            test_vec_G2_field_elems_1_0[test_num],
            test_vec_G2_field_elems_1_1[test_num]);

    ret = ret && fp2_equal(&expect, elems[1]);
    printf("u[1]\nActual:\n");
    printf("\ta: %s\n", mpz_get_str(NULL, 16, elems[1]->a));
    printf("\tb: %s\n", mpz_get_str(NULL, 16, elems[1]->b));
    printf("Expected:\n");
    printf("\ta: %s\n", mpz_get_str(NULL, 16, expect.a));
    printf("\tb: %s\n", mpz_get_str(NULL, 16, expect.b));
    fp2_elem_free(&expect);

    fp2_elem_free(elems[0]);
    fp2_elem_free(elems[1]);
    free(elems);

    return ret;
}

bool test_map_to_curve_G1(int test_num)
{
    printf("Running %s #%d\n", __func__, test_num);

    bool ret;
    mpz_t u;
    G1_elem_affine actual, expect;

    mpz_init(u);
    G1_identity_init_affine(&actual);

    mpz_set_str(u, test_vec_G1_field_elems_0[test_num], 0);
    G1_elem_affine_from_str(&expect,
                           test_vec_G1_mapped_points_0_x[test_num],
                           test_vec_G1_mapped_points_0_y[test_num]);

    map_to_curve_G1(&actual, u);

    ret = G1_equiv_affine(&actual, &expect);
    printf("Actual:\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, actual.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, actual.y));
    printf("Expected:\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, expect.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, expect.y));
    G1_elem_free_affine(&expect);

    mpz_set_str(u, test_vec_G1_field_elems_1[test_num], 0);
    G1_elem_affine_from_str(&expect,
                           test_vec_G1_mapped_points_1_x[test_num],
                           test_vec_G1_mapped_points_1_y[test_num]);

    map_to_curve_G1(&actual, u);

    ret = ret && G1_equiv_affine(&actual, &expect);
    printf("Actual:\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, actual.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, actual.y));
    printf("Expected:\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, expect.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, expect.y));
    G1_elem_free_affine(&expect);

    G1_elem_free_affine(&actual);
    mpz_clear(u);

    return ret;
}

bool test_hash_to_G1(int test_num)
{
    printf("Running %s #%d\n", __func__, test_num);

    bool ret;
    G1_elem_affine actual, expect;

    G1_identity_init_affine(&actual);

    G1_elem_affine_from_str(&expect,
                           test_vec_G1_hash_out_x[test_num],
                           test_vec_G1_hash_out_y[test_num]);

    BLS12_381_hash_to_G1(&actual, test_vec_msg[test_num], G1_TEST_DST);

    ret = G1_equiv_affine(&actual, &expect);
    printf("Actual:\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, actual.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, actual.y));
    printf("Expected:\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, expect.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, expect.y));

    G1_elem_free_affine(&expect);
    G1_elem_free_affine(&actual);

    return ret;
}

int main()
{
    TEST_MAIN_INIT;

    fp_params_init();

    for (int i = 0; i < 10; i++) {
        TEST_RUN(test_message_expand_xmd(i));
    }

    for (int i = 0; i < 5; i++) {
        TEST_RUN(test_hash_to_field_fp(i));
    }

    for (int i = 0; i < 5; i++) {
        TEST_RUN(test_hash_to_field_fp2(i));
    }

    for (int i = 0; i < 5; i++) {
        TEST_RUN(test_map_to_curve_G1(i));
    }

    for (int i = 0; i < 5; i++) {
        TEST_RUN(test_hash_to_G1(i));
    }

    fp_params_free();

    TEST_MAIN_RETURN;
}
