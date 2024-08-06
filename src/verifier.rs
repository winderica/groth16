use ark_ec::{pairing::Pairing, AffineRepr, CurveGroup, VariableBaseMSM};

use crate::{link::PESubspaceSnark, r1cs_to_qap::R1CSToQAP, Groth16};

use super::{PreparedVerifyingKey, Proof, VerifyingKey};

use ark_relations::r1cs::{Result as R1CSResult, SynthesisError};
use ark_std::One;

use core::ops::Neg;

/// Prepare the verifying key `vk` for use in proof verification.
pub fn prepare_verifying_key<E: Pairing>(vk: &VerifyingKey<E>) -> PreparedVerifyingKey<E> {
    PreparedVerifyingKey {
        vk: vk.clone(),
        alpha_g1_beta_g2: E::pairing(vk.alpha_g1, vk.beta_g2),
        gamma_g2_neg_pc: vk.gamma_g2.into_group().neg().into_affine().into(),
        delta_g2_neg_pc: vk.delta_g2.into_group().neg().into_affine().into(),
    }
}

impl<E: Pairing, QAP: R1CSToQAP> Groth16<E, QAP> {
    /// Prepare proof inputs for use with [`verify_proof_with_prepared_inputs`],
    /// wrt the prepared verification key `pvk` and instance public inputs.
    pub fn prepare_inputs(
        pvk: &PreparedVerifyingKey<E>,
        public_inputs: &[E::ScalarField],
    ) -> R1CSResult<E::G1> {
        let mut public_inputs = public_inputs.to_vec();
        public_inputs.insert(0, E::ScalarField::one());
        if public_inputs.len() != pvk.vk.gamma_abc_g1.0.len() {
            return Err(SynthesisError::MalformedVerifyingKey);
        }

        Ok(E::G1::msm_unchecked(&pvk.vk.gamma_abc_g1.0, &public_inputs))
    }

    /// Verify a Groth16 proof `proof` against the prepared verification key
    /// `pvk` and prepared public inputs. This should be preferred over
    /// [`verify_proof`] if the instance's public inputs are
    /// known in advance.
    pub fn verify_proof_with_prepared_inputs(
        pvk: &PreparedVerifyingKey<E>,
        (proof, link_d): &(Proof<E>, Vec<E::G1Affine>),
        prepared_inputs: &E::G1,
    ) -> R1CSResult<bool> {
        let qap = E::multi_miller_loop(
            [
                <E::G1Affine as Into<E::G1Prepared>>::into(proof.a),
                (proof.d + prepared_inputs).into_affine().into(),
                proof.c.into(),
            ],
            [
                proof.b.into(),
                pvk.gamma_g2_neg_pc.clone(),
                pvk.delta_g2_neg_pc.clone(),
            ],
        );

        let test = E::final_exponentiation(qap).ok_or(SynthesisError::UnexpectedIdentity)?;

        Ok(test == pvk.alpha_g1_beta_g2
            && PESubspaceSnark::<E>::verify(
                &pvk.vk.link_pp,
                &pvk.vk.link_vk,
                &[link_d.clone(), vec![proof.d]].concat(),
                &proof.link_pi,
            ))
    }

    /// Verify a Groth16 proof `proof` against the prepared verification key
    /// `pvk`, with respect to the instance `public_inputs`.
    pub fn verify_proof(
        pvk: &PreparedVerifyingKey<E>,
        proof: &(Proof<E>, Vec<E::G1Affine>),
        public_inputs: &[E::ScalarField],
    ) -> R1CSResult<bool> {
        let prepared_inputs = Self::prepare_inputs(pvk, public_inputs)?;
        Self::verify_proof_with_prepared_inputs(pvk, proof, &prepared_inputs)
    }
}
