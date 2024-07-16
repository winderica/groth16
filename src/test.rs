use std::sync::Arc;

use crate::{prepare_verifying_key, Groth16};
use ark_crypto_primitives::snark::{CircuitSpecificSetupSNARK, SNARK};
use ark_ec::{pairing::Pairing, CurveGroup, VariableBaseMSM};
use ark_ff::Field;
use ark_relations::{
    lc,
    r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError},
};
use ark_std::{
    rand::{RngCore, SeedableRng},
    test_rng, UniformRand, Zero,
};

struct MySillyCircuit<F: Field> {
    a: Option<F>,
    b: Option<F>,
    c: Option<F>,
    d: Option<F>,
    e: Option<F>,
    f: Option<F>,
}

impl<ConstraintF: Field> ConstraintSynthesizer<ConstraintF> for MySillyCircuit<ConstraintF> {
    fn generate_constraints(
        self,
        cs: ConstraintSystemRef<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        let a = cs.new_witness_variable(|| self.a.ok_or(SynthesisError::AssignmentMissing))?;
        let b = cs.new_witness_variable(|| self.b.ok_or(SynthesisError::AssignmentMissing))?;
        let c = cs.new_witness_variable(|| self.c.ok_or(SynthesisError::AssignmentMissing))?;
        let d = cs.new_committed_variable(|| self.d.ok_or(SynthesisError::AssignmentMissing))?;
        let e = cs.new_committed_variable(|| self.e.ok_or(SynthesisError::AssignmentMissing))?;
        let f = cs.new_committed_variable(|| self.f.ok_or(SynthesisError::AssignmentMissing))?;
        let g = cs.new_input_variable(|| {
            let a = self.a.ok_or(SynthesisError::AssignmentMissing)?;
            let b = self.b.ok_or(SynthesisError::AssignmentMissing)?;
            let c = self.c.ok_or(SynthesisError::AssignmentMissing)?;
            let d = self.d.ok_or(SynthesisError::AssignmentMissing)?;
            let e = self.e.ok_or(SynthesisError::AssignmentMissing)?;
            let f = self.f.ok_or(SynthesisError::AssignmentMissing)?;

            Ok((a + b + c) * (d + e + f))
        })?;

        cs.enforce_constraint(lc!() + a + b + c, lc!() + d + e + f, lc!() + g)?;
        cs.enforce_constraint(lc!() + a + b + c, lc!() + d + e + f, lc!() + g)?;
        cs.enforce_constraint(lc!() + a + b + c, lc!() + d + e + f, lc!() + g)?;
        cs.enforce_constraint(lc!() + a + b + c, lc!() + d + e + f, lc!() + g)?;
        cs.enforce_constraint(lc!() + a + b + c, lc!() + d + e + f, lc!() + g)?;
        cs.enforce_constraint(lc!() + a + b + c, lc!() + d + e + f, lc!() + g)?;
        cs.enforce_constraint(lc!() + a + b + c, lc!() + d + e + f, lc!() + g)?;
        cs.enforce_constraint(lc!() + a + b + c, lc!() + d + e + f, lc!() + g)?;
        cs.enforce_constraint(lc!() + a + b + c, lc!() + d + e + f, lc!() + g)?;
        cs.enforce_constraint(lc!() + a + b + c, lc!() + d + e + f, lc!() + g)?;

        Ok(())
    }
}

fn test_prove_and_verify<E>(n_iters: usize)
where
    E: Pairing,
{
    let mut rng = ark_std::rand::rngs::StdRng::seed_from_u64(test_rng().next_u64());

    let generators = vec![E::G1Affine::rand(&mut rng), E::G1Affine::rand(&mut rng)];
    let g = &generators[..1];
    let r = E::G1Affine::rand(&mut rng);

    let pk = Groth16::<E>::generate_random_parameters_with_reduction(
        MySillyCircuit {
            a: None,
            b: None,
            c: None,
            d: None,
            e: None,
            f: None,
        },
        vec![
            (&g, r),
            (&g, r),
            (&g, r),
        ],
        &mut rng,
    )
    .unwrap();
    let vk = pk.vk.clone();
    let pvk = prepare_verifying_key::<E>(&vk);

    for _ in 0..n_iters {
        let a = E::ScalarField::rand(&mut rng);
        let b = E::ScalarField::rand(&mut rng);
        let c = E::ScalarField::rand(&mut rng);
        let d = E::ScalarField::rand(&mut rng);
        let e = E::ScalarField::rand(&mut rng);
        let f = E::ScalarField::rand(&mut rng);
        let g = (a + b + c) * (d + e + f);

        let (proof, _) = Groth16::<E>::prove(
            &pk,
            MySillyCircuit {
                a: Some(a),
                b: Some(b),
                c: Some(c),
                d: Some(d),
                e: Some(e),
                f: Some(f),
            },
            &mut rng,
        )
        .unwrap();

        let link_d = vec![
            E::G1::msm(&generators, &[d, E::ScalarField::zero()])
                .unwrap()
                .into_affine(),
            E::G1::msm(&generators, &[e, E::ScalarField::zero()])
                .unwrap()
                .into_affine(),
                E::G1::msm(&generators, &[f, E::ScalarField::zero()])
                .unwrap()
                .into_affine()
        ];
        let proof = (proof, link_d);

        assert!(Groth16::<E>::verify_with_processed_vk(&pvk, &[g], &proof).unwrap());
        assert!(!Groth16::<E>::verify_with_processed_vk(&pvk, &[a], &proof).unwrap());
    }
}

mod bls12_377 {
    use super::test_prove_and_verify;
    use ark_bls12_377::Bls12_377;

    #[test]
    fn prove_and_verify() {
        test_prove_and_verify::<Bls12_377>(5);
    }
}

mod bn254 {
    use super::test_prove_and_verify;

    use ark_bn254::Bn254;

    #[test]
    fn prove_and_verify() {
        test_prove_and_verify::<Bn254>(5);
    }
}

mod bls12_381 {
    use super::test_prove_and_verify;

    use ark_bls12_381::Bls12_381;

    #[test]
    fn prove_and_verify() {
        test_prove_and_verify::<Bls12_381>(5);
    }
}

mod bn_254 {
    use super::test_prove_and_verify;
    use ark_bn254::Bn254;

    #[test]
    fn prove_and_verify() {
        test_prove_and_verify::<Bn254>(100);
    }
}
