#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "hash.h"

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

int main()
{
    TEST_MAIN_INIT;

    for (int i = 0; i < 10; i++) {
        TEST_RUN(test_message_expand_xmd(i));
    }

    for (int i = 0; i < 5; i++) {
        TEST_RUN(test_hash_to_field_fp(i));
    }

    TEST_MAIN_RETURN;
}
