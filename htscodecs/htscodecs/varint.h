#ifndef VARINT_H
#define VARINT_H

#include <stdint.h>

static inline int u32tou7(uint8_t *cp, uint32_t i) {
    uint8_t *op = cp;

    do {
	*cp++ = (i&0x7f) + ((i>=0x80)<<7);
	i >>= 7;
    } while (i);

    return cp-op;
}

static inline int u7tou32(uint8_t *cp, uint8_t *cp_end, uint32_t *i) {
    uint8_t *op = cp, c;
    uint32_t j = 0, s = 0;

    if (cp >= cp_end) {
	*i = 0;
	return 0;
    }

    do {
	c = *cp++;
	j |= (c & 0x7f) << s;
	s += 7;
    } while ((c & 0x80) && cp < cp_end);

    *i = j;
    return cp-op;
}

#endif /* VARINT_H */
