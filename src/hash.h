#ifndef _HASH_H_
#define _HASH_H_

#include <stdint.h>

#include <octet_string.h>

int expand_message_xmd(octet_string *bytes, const char *msg, const char *DST, uint32_t len_in_bytes);

#endif /* _HASH_H_ */
