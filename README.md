# libBLS12381

## Overview

### WARNING: This library has not been reviewed or audited. It was written for educational purposes.

This library implements the curve BLS12-381. BLS12-381 is a pairing-friendly
curve from the BLS (Barreto, Lynn, Scott) family of curves with embedding
degree k=12. BLS12 curves are parameterized by the polynomials

    p(t) = (t - 1)^2 * (t^4 - t^2 + 1) / 3 + t
    r(t) = t^4 - t^2 + 1.

The curve equation is

    E: y^2 = x^3 + 4,

which admits a a sextic twist

    E': y^2 = x^3 + 4(u + 1) where u

is neither a quadratic- nor cubic- resiude mod p, and p is the prime
characteristic of finite field F\_p.

Our pairing is defined e: G\_2 x G\_1 -> G\_T, where G\_1, G\_2 and G\_T each have
prime order r.  G\_1 is defined over the finite field F\_p, G\_2 is defined over
the quadratic extension F\_p^2, and G\_T is a subgroup of the full extension
field F\_p^12. The fields are constructed using the following tower:

     GF(p^2)  = GF(p)[u] / (u^2 + 1)
     GF(p^6)  = GF(p^2)[v] / (v^3 - u - 1)
     GF(p^12) = GF(p^6)[w] / (w^2 - v)

Thus, we have G\_1 is a subgroup of E(F\_p), G\_2 is a subgroup of E'(F\_p^2), and
G\_T is the r roots of unity from the multiplicative group of F\_p^12.

BLS12-381's sextic twist is M-type. This means that the twisting isomorphism
from E -> E' is more efficient than the untwisting isomorphism from E' -> E.
For this reason, rather than untwisting points on G\_2 during the pairing
computation, we twist points in G\_1, thus computing the whole pairing on the
curve. This is called the twisted ate pairing.

See [`include/BLS12_381.h`]() for curve parameters and more.

## Compiling

libBLS12381 depends on `libgmp` for multi-precision integer arithmetic, and `libcrypto` for SHA-256.

To build, run `make` from the top-level directory. The unit tests can be run with `make -C test`.

## Usage

There is only one header, `BLS12_381.h`. Before calling other library
functions, the `BLS12_381_init()` function must be called. When the curve is no
longer needed, `BLS12_381_free()` must be called.

### Working with Points in `G_1` and `G_2`

For simplicity, this section will only refer to elements in `G_1` because the APIs
are identical. For example, the `BLS12_381_hash_to_G1` has a corresponding
`BLS12_381_hash_to_G2`.

Points, or elements, in `G_1` can be represented in either affine or projective
coordinates, which are represented by the `G1_elem_affine` and `G1_elem_proj`
types. Recall that the projetive coordinates `(X, Y, Z)` correspond to the
affine coordinates `(X/Z, Y/Z)`.

There are several ways in the API to initialize an element in `G_1`. The simplest approach
is to initialize a point as the point at infinity (i.e. the identity element). This can be
done with the `G1_identity_init_affine` and `G1_identity_init_proj` functions. Another simple
approach is to initialize an element as a constant generator (see `BLS12_381.h` for the parameters).
This is done with `G1_generator_init_affine` and `G1_generator_init_proj`.

If a set of coordinates are known to give a point in `G_1`, a point can be initialized with
`G1_elem_affine_from_str` and `G1_elem_proj_from_str`.

Finally, `BLS12_381_hash_to_G1` implements the [IETF hash-to-curve draft-standard](https://datatracker.ietf.org/doc/draft-irtf-cfrg-hash-to-curve/). This provides a way to hash strings to a point on `G_1`.

Addition and scalar multiplication is done with the `G1_mul_scalar`,
`G1_add_proj`, `G1_add_mixed`, and `G1_double_proj` functions. See
[`include/BLS12_381.h`]() for the complete API.

### Pairing Computation

Pairing computation is done with `BLS12_381_pairing`. This gives an element in the field `F_p^12`. Field arithmetic
operations are provided with `fp12_add`, `fp12_mul`, etc. See [`include/BLS12_381.h`]() for the complete API.

## Note on Constant Time

As this library was initially written for educational purposes, the
implementation favored straight-forward, correct code over optimized or
constant-time code. As such, several operations (including scalar
multiplication) are not constant-time.  Future work will be focused on constant-time
re-factoring, and optimizations.
