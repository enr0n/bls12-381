#ifndef _HASH_H_
#define _HASH_H_

#include <stdint.h>

#include <gmp.h>

#include "octet_string.h"

int expand_message_xmd(octet_string *bytes, const char *msg, const char *DST, uint32_t len_in_bytes);

mpz_t *hash_to_field_fp(const char *msg, const char *DST, uint32_t count);
fp2_elem **hash_to_field_fp2(const char *msg, const char *DST, uint32_t count);

void map_to_curve_simple_swu_3mod4(mpz_t xn, mpz_t xd, mpz_t y, const mpz_t u);
void map_to_curve_G1(G1_elem_affine *P, const mpz_t u);

#endif /* _HASH_H_ */
