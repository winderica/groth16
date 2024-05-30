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
}

impl<ConstraintF: Field> ConstraintSynthesizer<ConstraintF> for MySillyCircuit<ConstraintF> {
    fn generate_constraints(
        self,
        cs: ConstraintSystemRef<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        let a = cs.new_witness_variable(|| self.a.ok_or(SynthesisError::AssignmentMissing))?;
        let b = cs.new_witness_variable(|| self.b.ok_or(SynthesisError::AssignmentMissing))?;
        let c = cs.new_committed_variable(|| self.c.ok_or(SynthesisError::AssignmentMissing))?;
        let d = cs.new_committed_variable(|| self.d.ok_or(SynthesisError::AssignmentMissing))?;
        let e = cs.new_input_variable(|| {
            let a = self.a.ok_or(SynthesisError::AssignmentMissing)?;
            let b = self.b.ok_or(SynthesisError::AssignmentMissing)?;
            let c = self.c.ok_or(SynthesisError::AssignmentMissing)?;
            let d = self.d.ok_or(SynthesisError::AssignmentMissing)?;

            Ok((a + b) * (c + d))
        })?;

        cs.enforce_constraint(lc!() + a + b, lc!() + c + d, lc!() + e)?;
        cs.enforce_constraint(lc!() + a + b, lc!() + c + d, lc!() + e)?;
        cs.enforce_constraint(lc!() + a + b, lc!() + c + d, lc!() + e)?;
        cs.enforce_constraint(lc!() + a + b, lc!() + c + d, lc!() + e)?;
        cs.enforce_constraint(lc!() + a + b, lc!() + c + d, lc!() + e)?;
        cs.enforce_constraint(lc!() + a + b, lc!() + c + d, lc!() + e)?;

        Ok(())
    }
}

fn test_prove_and_verify<E>(n_iters: usize)
where
    E: Pairing,
{
    let mut rng = ark_std::rand::rngs::StdRng::seed_from_u64(test_rng().next_u64());

    let (pk, vk) = Groth16::<E>::setup(
        MySillyCircuit {
            a: None,
            b: None,
            c: None,
            d: None,
        },
        &mut rng,
    )
    .unwrap();
    let pvk = prepare_verifying_key::<E>(&vk);

    for _ in 0..n_iters {
        let a = E::ScalarField::rand(&mut rng);
        let b = E::ScalarField::rand(&mut rng);
        let c = E::ScalarField::rand(&mut rng);
        let d = E::ScalarField::rand(&mut rng);
        let e = (a + b) * (c + d);

        let proof = Groth16::<E>::prove(
            &pk,
            MySillyCircuit {
                a: Some(a),
                b: Some(b),
                c: Some(c),
                d: Some(d),
            },
            &mut rng,
        )
        .unwrap();

        assert_eq!(
            proof.link_d,
            E::G1::msm(&pk.common.link_bases, &[c, d, E::ScalarField::zero()])
                .unwrap()
                .into_affine()
        );

        assert!(Groth16::<E>::verify_with_processed_vk(&pvk, &[e], &proof).unwrap());
        assert!(!Groth16::<E>::verify_with_processed_vk(&pvk, &[a], &proof).unwrap());
    }
}

mod bls12_377 {
    use super::test_prove_and_verify;
    use ark_bls12_377::Bls12_377;

    #[test]
    fn prove_and_verify() {
        test_prove_and_verify::<Bls12_377>(100);
    }
}

mod bw6_761 {
    use super::test_prove_and_verify;

    use ark_bw6_761::BW6_761;

    #[test]
    fn prove_and_verify() {
        test_prove_and_verify::<BW6_761>(1);
    }
}
