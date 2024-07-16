//! zkSNARK for Linear Subspaces as defined in appendix D of the paper.
//! Use to prove knowledge of openings of multiple Pedersen commitments. Can
//! also prove knowledge and equality of committed values in multiple
//! commitments. Note that this SNARK requires a trusted setup as the key
//! generation creates a trapdoor.

use crate::link::utils::*;
use ark_ec::{pairing::Pairing, AffineRepr, CurveGroup};
use ark_ff::{UniformRand, Zero};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    cfg_iter,
    marker::PhantomData,
    ops::{Mul, Neg},
    rand::Rng,
    vec::Vec,
};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Public params
#[derive(Clone, Default, PartialEq, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct PP<
    G1: Clone + Default + CanonicalSerialize + CanonicalDeserialize,
    G2: Clone + Default + CanonicalSerialize + CanonicalDeserialize,
> {
    pub l: usize, // # of rows
    pub t: usize, // # of cols
    pub g1: G1,
    pub g2: G2,
}

impl<
        G1: Clone + Default + CanonicalSerialize + CanonicalDeserialize,
        G2: Clone + Default + CanonicalSerialize + CanonicalDeserialize,
    > PP<G1, G2>
{
    pub fn new(l: usize, t: usize, g1: G1, g2: G2) -> PP<G1, G2> {
        PP { l, t, g1, g2 }
    }
}

/// Evaluation key
#[derive(Clone, Default, PartialEq, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct EK<G1: Clone + Default + CanonicalSerialize + CanonicalDeserialize> {
    pub p: Vec<G1>,
}

/// Verification key
#[derive(Clone, Default, PartialEq, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct VK<G2: Clone + Default + CanonicalSerialize + CanonicalDeserialize> {
    pub c: Vec<G2>,
    pub a: G2,
}

pub struct PESubspaceSnark<PE: Pairing> {
    pairing_engine_type: PhantomData<PE>,
}

// NB: Now the system is for y = Mx
impl<PE: Pairing> PESubspaceSnark<PE> {
    /// Matrix should be such that a column will have more than 1 non-zero item
    /// only if those values are equal. Eg for matrix below, h2 and h3
    /// commit to same value h1, 0, 0, 0
    /// 0, h2, 0, 0
    /// 0, h3, h4, 0
    pub fn keygen<R: Rng>(
        rng: &mut R,
        pp: &PP<PE::G1Affine, PE::G2Affine>,
        m: &[(usize, &[PE::G1Affine], usize, PE::G1Affine)],
    ) -> (EK<PE::G1Affine>, VK<PE::G2Affine>) {
        // `t` is the trapdoor
        let mut t: Vec<PE::ScalarField> = Vec::with_capacity(pp.l);
        for _ in 0..pp.l {
            t.push(PE::ScalarField::rand(rng));
        }

        let a = PE::ScalarField::rand(rng);

        let p = PE::G1::normalize_batch(
            &cfg_into_iter!(0..pp.t)
                .map(|i| {
                    m.iter().enumerate().map(|(u, (j, v, k, r))| {
                        if i >= *j && i < *j + v.len() {
                            v[i - j].mul(&t[u])
                        } else if i == *k {
                            r.mul(&t[u])
                        } else {
                            PE::G1::zero()
                        }
                    }).sum::<PE::G1>()
                })
                .collect::<Vec<_>>(),
        );

        let c = scale_vector::<PE>(a, &t);
        let ek = EK::<PE::G1Affine> { p };
        let vk = VK::<PE::G2Affine> {
            c: multiples_of_g::<PE::G2Affine>(&pp.g2, &c),
            a: pp.g2.mul(a).into_affine(),
        };
        (ek, vk)
    }

    pub fn prove(
        pp: &PP<PE::G1Affine, PE::G2Affine>,
        ek: &EK<PE::G1Affine>,
        w: &[PE::ScalarField],
    ) -> PE::G1Affine {
        assert!(pp.t >= w.len());
        inner_product::<PE>(w, &ek.p)
    }

    pub fn verify(
        pp: &PP<PE::G1Affine, PE::G2Affine>,
        vk: &VK<PE::G2Affine>,
        x: &[PE::G1Affine],
        pi: &PE::G1Affine,
    ) -> bool {
        assert_eq!(pp.l, x.len());
        assert!(vk.c.len() >= x.len());

        let mut a = x.to_vec();
        let mut b = cfg_iter!(vk.c[0..x.len()])
            .map(|b| PE::G2Prepared::from(*b))
            .collect::<Vec<_>>();
        a.push(*pi);
        b.push(PE::G2Prepared::from(vk.a.into_group().neg()));
        PE::multi_pairing(a, b).is_zero()
    }
}
