//! Utils for matrix and vector operations

use ark_ec::{pairing::Pairing, AffineRepr, CurveGroup, ScalarMul, VariableBaseMSM};
use rayon::prelude::*;

/// MSM between a scalar vector and a G1 vector
pub fn inner_product<PE: Pairing>(a: &[PE::ScalarField], b: &[PE::G1Affine]) -> PE::G1Affine {
    PE::G1::msm_unchecked(b, &a).into_affine()
}

/// Scale given vector `v` by scalar `a`
pub fn scale_vector<PE: Pairing>(
    a: PE::ScalarField,
    v: &[PE::ScalarField],
) -> Vec<PE::ScalarField> {
    cfg_iter!(v).map(|x| a * x).collect()
}

/// Given a group element `g` and vector `multiples` of scalars, returns a
/// vector with elements `v_i * g`
pub fn multiples_of_g<G: AffineRepr>(g: &G, multiples: &[G::ScalarField]) -> Vec<G> {
    g.into_group().batch_mul(multiples)
}
