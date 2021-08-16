#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "octet_string.h"

void octet_string_alloc(octet_string **po, uint32_t cap)
{
    octet_string *o;
    o = calloc(1, sizeof(octet_string*));

    o->data = calloc(cap, sizeof(uint8_t));
    o->cap = cap;
    o->len = 0;

    *po = o;
}

void octet_string_free(octet_string *o)
{
    if (!o) {
        return;
    }

    if (o->data) {
        free(o->data);
    }

    free(o);
}

void octet_string_reset(octet_string *o)
{
    for (int i = 0; i < (int)o->cap; i++) {
        o->data[i] = 0x0;
    }
    o->len = 0;
}

char *octet_string_to_str(const octet_string *o)
{
    char *str;

    str = calloc(2 * o->len + 1, sizeof(char));
    for (int i = 0; i < (int)o->len; i++) {
        sprintf(&str[2*i], "%02x", o->data[i]);
    }
    str[2 * o->len] = '\0';

    return str;
}

static inline void octet_string_realloc(octet_string *o, uint32_t cap)
{
    uint8_t *new;

    new = realloc(o->data, cap);
    assert(new != NULL);

    o->data = new;
    o->cap = cap;
}

octet_string *octet_string_append(octet_string *dest, uint8_t byte)
{
    assert(dest->cap >= dest->len);

    if (dest->cap == dest->len) {
        octet_string_realloc(dest, dest->cap + 1);
    }

    dest->data[dest->len] = byte;
    dest->len++;

    return dest;
}

octet_string *octet_string_appendn(octet_string *dest, const uint8_t *bytes, uint32_t n)
{
    if (dest->cap < dest->len + n) {
        octet_string_realloc(dest, dest->len + n);
    }

    for (int i = 0; i < (int)n; i++) {
        dest->data[dest->len + i] = bytes[i];
    }
    dest->len += n;

    return dest;
}

octet_string *octet_substr(octet_string *dest, const octet_string *src,
                           uint32_t start, uint32_t len)
{
    if (dest->cap < len) {
        octet_string_realloc(dest, len);
    }

    for (int i = 0; i < (int)len; i++) {
        dest->data[i] = src->data[start + i];
    }
    dest->len = len;

    return dest;
}

octet_string *octet_strcpy(octet_string *dest, const octet_string *src)
{
    if (src->len > dest->cap) {
        octet_string_realloc(dest, src->cap);
    }

    memcpy(dest->data, src->data, src->len);
    dest->len = src->len;

    assert(dest->cap >= dest->len);

    return dest;
}

octet_string *octet_strcat(octet_string *dest, const octet_string *src)
{
    if (src->len + dest->len > dest->cap) {
        octet_string_realloc(dest, src->len + dest->len);
    }

    for (int i = 0; i < (int)src->len; i++) {
        dest->data[dest->len + i] = src->data[i];
    }
    dest->len += src->len;

    assert(dest->cap >= dest->len);

    return dest;
}

octet_string *octet_strncpy(octet_string *dest, const octet_string *src, uint32_t n)
{
    if (n > dest->cap) {
        octet_string_realloc(dest, n);
    }

    memcpy(dest->data, src->data, n);
    dest->len = n;

    return dest;
}
