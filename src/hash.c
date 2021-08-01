#include <string.h>
#include <stdint.h>

#include <openssl/sha.h>

#include "octet_string.h"

#include "BLS12_381.h"

#define DST_MAX_LEN 255

#define B_IN_BYTES (256 / 8)
#define R_IN_BYTES 64

#define CEIL(x,y) 1 + ((x - 1) / y)

#define HASH(ctx,dst,src)                     \
    SHA256_Init(&ctx);                        \
    SHA256_Update(&ctx, src->data, src->len); \
    SHA256_Final(dst, &ctx);

octet_string *I2OSP(octet_string *o, int x, int x_len)
{
    octet_string *tmp;

    octet_string_alloc(&tmp, x_len);

    for (int i = 1; i <= x_len && x != 0; i++) {
        tmp->data[x_len - 1] = x % 256;
        x /= 256;
    }
    tmp->len = x_len;

    octet_strcpy(o, tmp);
    octet_string_free(tmp);

    return o;
}

int expand_message_xmd(octet_string *bytes, const char *msg, const char *DST, uint32_t len_in_bytes)
{
    uint32_t ell, DST_len;
    octet_string *DST_prime, *msg_prime, *uniform_bytes,
                 *tmp, *arg;

    SHA256_CTX sha_ctx;
    uint8_t b0[SHA256_DIGEST_LENGTH], bi[SHA256_DIGEST_LENGTH];

    ell = CEIL(len_in_bytes, B_IN_BYTES);
    if (ell > 255) {
        return 1;
    }

    DST_len = strnlen(DST, DST_MAX_LEN);

    octet_string_alloc(&uniform_bytes, SHA256_DIGEST_LENGTH * ell);
    octet_string_alloc(&DST_prime, DST_len + 2);
    octet_string_alloc(&msg_prime, R_IN_BYTES + DST_prime->len + strlen(msg) + 2);
    octet_string_alloc(&tmp, SHA256_DIGEST_LENGTH + DST_len + 1);
    octet_string_alloc(&arg, 0);

    octet_string_appendn(DST_prime, (uint8_t*)DST, DST_len);
    octet_strcat(DST_prime, I2OSP(arg, DST_len, 1));

    octet_strcat(msg_prime, I2OSP(arg, 0, R_IN_BYTES));
    octet_string_appendn(msg_prime, (uint8_t*)msg, strlen(msg));

    octet_strcat(msg_prime, I2OSP(arg, len_in_bytes, 2));
    octet_string_append(msg_prime, 0x0);
    octet_strcat(msg_prime, DST_prime);

    HASH(sha_ctx, b0, msg_prime);

    octet_string_appendn(tmp, b0, SHA256_DIGEST_LENGTH);
    octet_string_append(tmp, 0x1);
    octet_strcat(tmp, DST_prime);

    HASH(sha_ctx, bi, tmp);

    octet_string_appendn(uniform_bytes, bi, SHA256_DIGEST_LENGTH);

    if (ell < 2) goto out;

    for (int i = 2; i <= ell; i++) {
        octet_string_reset(tmp);

        for (int j = 0; j < SHA256_DIGEST_LENGTH; j++) {
            octet_string_append(tmp, b0[j] ^ bi[j]);
        }

        octet_strcat(tmp, I2OSP(arg, i, 1));
        octet_strcat(tmp, DST_prime);

        HASH(sha_ctx, bi, tmp);

        octet_string_appendn(uniform_bytes, bi, SHA256_DIGEST_LENGTH);
    }

out:
    octet_strncpy(bytes, uniform_bytes, len_in_bytes);

    octet_string_free(DST_prime);
    octet_string_free(msg_prime);
    octet_string_free(uniform_bytes);
    octet_string_free(tmp);
    octet_string_free(arg);

    return 0;
}
