#ifndef _PAIRING_H_
#define _PAIRING_H_

#include "finite_field.h"
#include "G1.h"
#include "G2.h"

void BLS12_381_init();
void BLS12_381_free();

void final_exponentiation(fp12_elem *r, const fp12_elem *f);

void miller_loop(fp12_elem *r, G2_elem_affine *R, const G1_elem_affine *P, const G2_elem_affine *Q);

void ate(fp12_elem *r, const G2_elem_affine *Q, const G1_elem_affine *P);
#endif /* _PAIRING_H_ */
