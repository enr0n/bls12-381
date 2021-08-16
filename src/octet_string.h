#ifndef _OCTET_STRING_H_
#define _OCTET_STRING_H_

#include <stdint.h>

typedef struct {
    uint8_t *data;

    uint32_t len;
    uint32_t cap;
} octet_string;

void octet_string_alloc(octet_string **o, uint32_t cap);
void octet_string_free(octet_string *o);
void octet_string_reset(octet_string *o);

char *octet_string_to_str(const octet_string *o);

octet_string *octet_string_append(octet_string *dest, uint8_t byte);
octet_string *octet_string_appendn(octet_string *dest, const uint8_t *bytes, uint32_t n);

octet_string *octet_substr(octet_string *dest, const octet_string *src,
                           uint32_t start, uint32_t len);

octet_string *octet_strcpy(octet_string *dest, const octet_string *src);
octet_string *octet_strcat(octet_string *dest, const octet_string *src);
octet_string *octet_strncpy(octet_string *dest, const octet_string *src, uint32_t n);

#endif /* _OCTET_STRING_H_ */
