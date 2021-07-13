#ifndef _PAIRING_H_
#define _PAIRING_H_

#include "finite_field.h"

void BLS12_381_init();
void BLS12_381_free();

void final_exponentiation(fp12_elem *r, const fp12_elem *f);

#endif /* _PAIRING_H_ */
