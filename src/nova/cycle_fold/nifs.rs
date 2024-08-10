#[cfg(test)]
mod tests {
    use ark_ff::PrimeField;
    use ark_pallas::{Fq, Fr, PallasConfig, Projective};
    use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef};
    use ark_std::UniformRand;
    use ark_vesta::VestaConfig;
    use rand::thread_rng;
    pub use crate::hash::pederson::PedersenCommitment;
    use crate::nova::commitment::CommitmentScheme;
    use crate::nova::cycle_fold::coprocessor::{Circuit, RelaxedR1CSInstance, setup_shape, synthesize};
    use crate::gadgets::r1cs::r1cs::{commit_T, commit_T_with_relaxed, RelaxedR1CSWitness};
    use crate::utils::cast_field_element;

    fn get_random_circuit() -> Circuit<PallasConfig> {
        let mut rng = ark_std::test_rng();
        let g1 = Projective::rand(&mut rng);
        let g2 = Projective::rand(&mut rng);

        let val = u64::rand(&mut rng);
        let r = <Fq as PrimeField>::BigInt::from(val).into();
        let r_scalar = unsafe { cast_field_element::<Fq, Fr>(&r) };

        let g_out = g1 * r_scalar + g2 ;

        Circuit {
            g1,
            g2,
            g_out,
            r,
        }
    }
    #[test]
    pub fn test_nifs_for_cycle_fold() {
        // define two circuits
        let c1 = get_random_circuit();
        let c2 = get_random_circuit();

        // define shape of the r1cs which is known, equal to 11
        let shape = setup_shape::<PallasConfig, VestaConfig>().unwrap();

        // pp which is Pederson's pp
        let pp = PedersenCommitment::<ark_vesta::Projective>::setup(shape.num_vars, b"test", &());

        // random challenge r
        let r: Fq = Fq::rand(&mut thread_rng());

        // counting the number of constraints on the second curve
        let cs = ConstraintSystem::<Fq>::new_ref();
        c1.clone().generate_constraints(cs.clone()).expect("error while generating constraints");
        println!("number of constraints on second curve: {}", cs.num_constraints());

        // define a running instance for second curve
        let (running_U, running_W) = synthesize::<
            PallasConfig,
            VestaConfig,
            PedersenCommitment<ark_vesta::Projective>,
        >(c1.clone(), &pp).unwrap();

        // turn it into relaxed instances
        let running_U = RelaxedR1CSInstance::from(&running_U);
        let running_W = RelaxedR1CSWitness::from_r1cs_witness(&shape, &running_W);

        // make sure it's satisfying
        shape.is_relaxed_satisfied(&running_U, &running_W, &pp).unwrap();

        // a new satisfying pair (u, w)
        let (u, w) = synthesize::<
            PallasConfig,
            VestaConfig,
            PedersenCommitment<ark_vesta::Projective>,
        >(c2, &pp).unwrap();

        // compute folding proof T, commitment_T
        let (T, commitment_T) = commit_T(&shape, &pp, &running_U, &running_W, &u, &w).unwrap();

        // fold instance/witnesses
        let folded_U = running_U.fold(&u, &commitment_T, &r).unwrap();
        let folded_W = running_W.fold(&w, &T, &r).unwrap();

        // check if they are satisfied
        shape.is_relaxed_satisfied(&folded_U, &folded_W, &pp).unwrap();
    }
}