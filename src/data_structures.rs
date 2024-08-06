use crate::link::{EK, PP, VK};
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_serialize::*;
use ark_std::vec::Vec;

/// A proof in the Groth16 SNARK
#[derive(Clone, Debug, PartialEq, CanonicalSerialize, CanonicalDeserialize, Default)]
pub struct Proof<E: Pairing> {
    /// The `A` element in `G1`.
    pub a: E::G1Affine,
    /// The `B` element in `G2`.
    pub b: E::G2Affine,
    /// The `C` element in `G1`.
    pub c: E::G1Affine,

    /// The `D` element in `G1`. Commits to a subset of private inputs of the
    /// circuit
    pub d: E::G1Affine,
    /// cp_{link}
    // pub link_d: Vec<E::G1Affine>,
    /// proof of commitment opening equality between `cp_{link}` and `d`
    pub link_pi: E::G1Affine,
}

////////////////////////////////////////////////////////////////////////////////

/// A verification key in the Groth16 SNARK.
#[derive(Clone, Debug, PartialEq, CanonicalSerialize, CanonicalDeserialize, Default)]
pub struct VerifyingKey<E: Pairing> {
    /// The `alpha * G`, where `G` is the generator of `E::G1`.
    pub alpha_g1: E::G1Affine,
    /// The `alpha * H`, where `H` is the generator of `E::G2`.
    pub beta_g2: E::G2Affine,
    /// The `gamma * H`, where `H` is the generator of `E::G2`.
    pub gamma_g2: E::G2Affine,
    /// The `delta * H`, where `H` is the generator of `E::G2`.
    pub delta_g2: E::G2Affine,
    /// The `gamma^{-1} * (beta * a_i + alpha * b_i + c_i) * H`, where `H` is
    /// the generator of `E::G1`.
    pub gamma_abc_g1: (Vec<E::G1Affine>, Vec<E::G1Affine>),
    /// The element `eta*gamma^-1 * G` in `E::G1`.
    pub eta_gamma_inv_g1: E::G1Affine,

    /// Public parameters of the Subspace Snark
    pub link_pp: PP<E::G1Affine, E::G2Affine>,
    /// Verification key of the Subspace Snark
    pub link_vk: VK<E::G2Affine>,
}

/// Preprocessed verification key parameters that enable faster verification
/// at the expense of larger size in memory.
#[derive(Clone, Debug, PartialEq, CanonicalSerialize, CanonicalDeserialize, Default)]
pub struct PreparedVerifyingKey<E: Pairing> {
    /// The unprepared verification key.
    pub vk: VerifyingKey<E>,
    /// The element `e(alpha * G, beta * H)` in `E::GT`.
    pub alpha_g1_beta_g2: PairingOutput<E>,
    /// The element `- gamma * H` in `E::G2`, prepared for use in pairings.
    pub gamma_g2_neg_pc: E::G2Prepared,
    /// The element `- delta * H` in `E::G2`, prepared for use in pairings.
    pub delta_g2_neg_pc: E::G2Prepared,
}

impl<E: Pairing> From<PreparedVerifyingKey<E>> for VerifyingKey<E> {
    fn from(other: PreparedVerifyingKey<E>) -> Self {
        other.vk
    }
}

////////////////////////////////////////////////////////////////////////////////

/// The common elements for Proving Key for with and without CP_link
#[derive(Clone, Debug, PartialEq, CanonicalSerialize, CanonicalDeserialize, Default)]
pub struct ProvingKeyCommon<E: Pairing> {
    /// The element `beta * G` in `E::G1`.
    pub beta_g1: E::G1Affine,
    /// The element `delta * G` in `E::G1`.
    pub delta_g1: E::G1Affine,
    /// The element `eta*delta^-1 * G` in `E::G1`.
    pub eta_delta_inv_g1: E::G1Affine,
    /// The elements `a_i * G` in `E::G1`.
    pub a_query: Vec<E::G1Affine>,
    /// The elements `b_i * G` in `E::G1`.
    pub b_g1_query: Vec<E::G1Affine>,
    /// The elements `b_i * H` in `E::G2`.
    pub b_g2_query: Vec<E::G2Affine>,
    /// The elements `h_i * G` in `E::G1`.
    pub h_query: Vec<E::G1Affine>,
    /// The elements `l_i * G` in `E::G1`.
    pub l_query: Vec<E::G1Affine>,

    /// Evaluation key of cp_{link}
    pub link_ek: EK<E::G1Affine>,
}

/// The prover key for for the Groth16 zkSNARK.
#[derive(Clone, Debug, PartialEq, CanonicalSerialize, CanonicalDeserialize, Default)]
pub struct ProvingKey<E: Pairing> {
    /// The underlying verification key.
    pub vk: VerifyingKey<E>,
    pub common: ProvingKeyCommon<E>,
}
