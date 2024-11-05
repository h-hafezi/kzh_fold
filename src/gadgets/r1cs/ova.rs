use crate::commitment::{CommitmentScheme, Len};
use crate::gadgets::r1cs::r1cs::{Error, R1CSShape};
use ark_ec::CurveGroup;
use ark_ff::{AdditiveGroup};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelIterator,
};
use std::fmt;

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct OvaInstance<G: CurveGroup, C: CommitmentScheme<G>> {
    /// commitment to <witness||[0], [g]||[h]>
    pub commitment: C::Commitment,
    /// X is assumed to start with `u`.
    pub X: Vec<G::ScalarField>,
}

impl<G: CurveGroup, C: CommitmentScheme<G>> Clone for OvaInstance<G, C> {
    fn clone(&self) -> Self {
        Self {
            commitment: self.commitment,
            X: self.X.clone(),
        }
    }
}

impl<G: CurveGroup, C: CommitmentScheme<G>> fmt::Debug for OvaInstance<G, C>
where
    C::Commitment: fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("R1CSInstance")
            .field("commitment", &self.commitment)
            .field("X", &self.X)
            .finish()
    }
}

impl<G: CurveGroup, C: CommitmentScheme<G>> PartialEq for OvaInstance<G, C> {
    fn eq(&self, other: &Self) -> bool {
        self.commitment == other.commitment && self.X == other.X
    }
}

impl<G: CurveGroup, C: CommitmentScheme<G>> Eq for OvaInstance<G, C> where C::Commitment: Eq {}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct OvaWitness<G: CurveGroup> {
    pub W: Vec<G::ScalarField>,
}

/// in essence, it's the as OvaInstance, it's only for better prototypes
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct RelaxedOvaInstance<G: CurveGroup, C: CommitmentScheme<G>> {
    pub commitment: C::Commitment,
    /// X is assumed to start with `u`.
    pub X: Vec<G::ScalarField>,
}

impl<G: CurveGroup, C: CommitmentScheme<G>> Clone for RelaxedOvaInstance<G, C> {
    fn clone(&self) -> Self {
        Self {
            commitment: self.commitment,
            X: self.X.clone(),
        }
    }
}

impl<G: CurveGroup, C: CommitmentScheme<G>> fmt::Debug for RelaxedOvaInstance<G, C>
where
    C::Commitment: fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("RelaxedOvaInstance")
            .field("commitment", &self.commitment)
            .field("X", &self.X)
            .finish()
    }
}

impl<G: CurveGroup, C: CommitmentScheme<G>> PartialEq for RelaxedOvaInstance<G, C> {
    fn eq(&self, other: &Self) -> bool {
        self.commitment == other.commitment
            && self.X == other.X
    }
}

impl<G: CurveGroup, C: CommitmentScheme<G>> Eq for RelaxedOvaInstance<G, C> where C::Commitment: Eq {}


#[derive(Default, Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct RelaxedOvaWitness<G: CurveGroup> {
    pub W: Vec<G::ScalarField>,
    pub E: Vec<G::ScalarField>,
}

impl<G: CurveGroup> R1CSShape<G> {
    pub fn is_ova_satisfied<C: CommitmentScheme<G>>(
        &self,
        U: &OvaInstance<G, C>,
        W: &OvaWitness<G>,
        pp: &C::PP,
    ) -> Result<(), Error> {
        assert_eq!(W.W.len(), self.num_vars);
        assert_eq!(U.X.len(), self.num_io);
        assert_eq!(pp.len(), self.num_vars + self.num_constraints);

        let z = [U.X.as_slice(), W.W.as_slice()].concat();
        let Az = self.A.multiply_vec(&z);
        let Bz = self.B.multiply_vec(&z);
        let Cz = self.C.multiply_vec(&z);

        if ark_std::cfg_into_iter!(0..self.num_constraints).any(|idx| Az[idx] * Bz[idx] != Cz[idx])
        {
            return Err(Error::NotSatisfied);
        }

        // commit to W || E
        let concat: Vec<G::ScalarField> = {
            let mut res = Vec::new();
            res.extend(W.W.clone());
            res.extend(vec![G::ScalarField::ZERO; self.num_constraints]);
            res
        };

        if U.commitment != C::commit(pp, &concat) {
            return Err(Error::NotSatisfied);
        }

        Ok(())
    }

    pub fn is_relaxed_ova_satisfied<C: CommitmentScheme<G>>(
        &self,
        U: &RelaxedOvaInstance<G, C>,
        W: &RelaxedOvaWitness<G>,
        pp: &C::PP,
    ) -> Result<(), Error> {
        assert_eq!(W.W.len(), self.num_vars);
        assert_eq!(U.X.len(), self.num_io);
        assert_eq!(W.E.len(), self.num_constraints);
        assert_eq!(pp.len(), self.num_vars + self.num_constraints);

        let z = [U.X.as_slice(), W.W.as_slice()].concat();
        let Az = self.A.multiply_vec(&z);
        let Bz = self.B.multiply_vec(&z);
        let Cz = self.C.multiply_vec(&z);

        let u = U.X[0];

        if ark_std::cfg_into_iter!(0..self.num_constraints)
            .any(|idx| Az[idx] * Bz[idx] != u * Cz[idx] + W.E[idx])
        {
            return Err(Error::NotSatisfied);
        }

        let concat = {
            let mut res = Vec::new();
            res.extend(W.W.clone());
            res.extend(W.E.clone());
            res
        };

        let commitment = C::commit(pp, concat.as_slice());

        if U.commitment != commitment {
            return Err(Error::NotSatisfied);
        }

        Ok(())
    }
}

impl<G: CurveGroup, C: CommitmentScheme<G>> RelaxedOvaInstance<G, C> {
    pub fn new(shape: &R1CSShape<G>) -> Self {
        Self {
            commitment: C::Commitment::default(),
            X: vec![G::ScalarField::ZERO; shape.num_io],
        }
    }

    /// Folds an incoming **non-relaxed** [`OvaInstance`] into the current one.
    pub fn fold(
        &self,
        U2: &OvaInstance<G, C>,
        com_T: &C::Commitment,
        r: &G::ScalarField,
    ) -> Result<Self, Error> {
        let (X1, comm_W1) = (&self.X, self.commitment);
        let (X2, comm_W2) = (&U2.X, &U2.commitment);

        let X: Vec<G::ScalarField> = ark_std::cfg_iter!(X1)
            .zip(X2)
            .map(|(a, b)| *a + *r * *b)
            .collect();
        let commitment_W = comm_W1 + (*comm_W2 + *com_T) * *r;

        Ok(Self { commitment: commitment_W, X })
    }
}

pub fn commit_T<G: CurveGroup, C: CommitmentScheme<G>>(
    shape: &R1CSShape<G>,
    pp: &C::PP,
    U1: &RelaxedOvaInstance<G, C>,
    W1: &RelaxedOvaWitness<G>,
    U2: &OvaInstance<G, C>,
    W2: &OvaWitness<G>,
) -> Result<(Vec<G::ScalarField>, C::Commitment), Error> {
    assert_eq!(pp.len(), shape.num_constraints);

    let z1 = [&U1.X, &W1.W[..]].concat();
    let Az1 = shape.A.multiply_vec(&z1);
    let Bz1 = shape.B.multiply_vec(&z1);
    let Cz1 = shape.C.multiply_vec(&z1);

    let z2 = [&U2.X, &W2.W[..]].concat();
    let Az2 = shape.A.multiply_vec(&z2);
    let Bz2 = shape.B.multiply_vec(&z2);
    let Cz2 = shape.C.multiply_vec(&z2);

    // Circle-product.
    let Az1_Bz2: Vec<G::ScalarField> = ark_std::cfg_iter!(Az1)
        .zip(&Bz2)
        .map(|(&a, &b)| a * b)
        .collect();
    let Az2_Bz1: Vec<G::ScalarField> = ark_std::cfg_iter!(Az2)
        .zip(&Bz1)
        .map(|(&a, &b)| a * b)
        .collect();

    // Scalar product.
    // u2 = 1 since U2 is non-relaxed instance, thus no multiplication required for Cz1.
    let u1 = U1.X[0];
    let u1_Cz2: Vec<G::ScalarField> = ark_std::cfg_into_iter!(Cz2).map(|cz2| u1 * cz2).collect();

    // Compute cross-term.
    let T: Vec<G::ScalarField> = ark_std::cfg_into_iter!(0..Az1_Bz2.len())
        .map(|i| Az1_Bz2[i] + Az2_Bz1[i] - u1_Cz2[i] - Cz1[i])
        .collect();

    let comm_T = C::commit(pp, T.as_slice());

    Ok((T, comm_T))
}

impl<G: CurveGroup, C: CommitmentScheme<G>> From<&OvaInstance<G, C>> for RelaxedOvaInstance<G, C> {
    fn from(instance: &OvaInstance<G, C>) -> Self {
        Self {
            commitment: instance.commitment,
            X: instance.X.clone(),
        }
    }
}

impl<G: CurveGroup> RelaxedOvaWitness<G> {
    pub fn zero(shape: &R1CSShape<G>) -> Self {
        Self {
            W: vec![G::ScalarField::ZERO; shape.num_vars],
            E: vec![G::ScalarField::ZERO; shape.num_constraints],
        }
    }

    /// Initializes a new [`RelaxedOvaWitness`] from an [`OvaWitness`].
    pub fn from(shape: &R1CSShape<G>, witness: &OvaWitness<G>) -> Self {
        Self {
            W: witness.W.clone(),
            E: vec![G::ScalarField::ZERO; shape.num_constraints],
        }
    }

    /// Folds an incoming **non-relaxed** [`OvaWitness`] into the current one.
    pub fn fold(
        &self,
        W2: &OvaWitness<G>,
        T: &[G::ScalarField],
        r: &G::ScalarField,
    ) -> Result<Self, Error> {
        let (W1, E1) = (&self.W, &self.E);
        let W2 = &W2.W;

        if W1.len() != W2.len() {
            return Err(Error::InvalidWitnessLength);
        }

        let W: Vec<G::ScalarField> = ark_std::cfg_iter!(W1)
            .zip(W2)
            .map(|(a, b)| *a + *r * *b)
            .collect();

        // Note that W2 is not relaxed, thus E2 = 0.
        let E: Vec<G::ScalarField> = ark_std::cfg_iter!(E1)
            .zip(T)
            .map(|(a, b)| *a + *r * *b)
            .collect();
        Ok(Self { W, E })
    }
}

#[cfg(test)]
pub(crate) mod tests {
    #![allow(non_upper_case_globals)]
    #![allow(clippy::needless_range_loop)]

    use crate::commitment::CommitmentScheme;
    use crate::gadgets::r1cs::ova::{commit_T, OvaInstance, OvaWitness, RelaxedOvaInstance, RelaxedOvaWitness};
    use crate::gadgets::r1cs::r1cs::tests::{to_field_elements, to_field_sparse, A, B, C};
    use crate::gadgets::r1cs::r1cs::{R1CSInstance, R1CSShape, R1CSWitness};
    use crate::hash::pederson::PedersenCommitment;
    use ark_bls12_381::{Fr, G1Projective as G};
    use ark_ff::Zero;
    use ark_std::UniformRand;
    use rand::thread_rng;

    #[test]
    fn is_satisfied() {
        let (a, b, c) = {
            (
                to_field_sparse::<G>(A),
                to_field_sparse::<G>(B),
                to_field_sparse::<G>(C),
            )
        };

        const NUM_CONSTRAINTS: usize = 4;
        const NUM_WITNESS: usize = 4;
        const NUM_PUBLIC: usize = 2;

        let pp = PedersenCommitment::<G>::setup(
            NUM_WITNESS + NUM_CONSTRAINTS,
            b"test",
            &(),
        );

        let shape = R1CSShape::<G>::new(
            NUM_CONSTRAINTS,
            NUM_WITNESS,
            NUM_PUBLIC,
            &a,
            &b,
            &c,
        ).unwrap();

        let X = to_field_elements::<G>(&[1, 35]);
        let W = to_field_elements::<G>(&[3, 9, 27, 30]);

        // zero vector E
        let E = vec![Fr::zero(); NUM_CONSTRAINTS];

        // commit to W || E
        let concat = {
            let mut res = Vec::new();
            res.extend(W.clone());
            res.extend(E.clone());
            res
        };

        let commitment = PedersenCommitment::<G>::commit(&pp, &concat);

        // R1CS instance/witness
        let instance = R1CSInstance::<G, PedersenCommitment<G>>::new(
            &shape,
            &commitment,
            &X,
        ).unwrap();
        let witness = R1CSWitness::<G>::new(&shape, &W).unwrap();

        // convert ova instance/witness
        let ova_instance = OvaInstance { commitment: instance.commitment_W, X: instance.X };


        let ova_witness = OvaWitness {
            W: witness.W,
        };

        // check ova instance/witness is satisfied
        shape.is_ova_satisfied::<PedersenCommitment<G>>(
            &ova_instance,
            &ova_witness,
            &pp,
        ).expect("satisfy failed");

        // convert to relaxed ova instance/witness
        let mut relaxed_ova_instance = RelaxedOvaInstance::from(&ova_instance);
        let mut relaxed_ova_witness = RelaxedOvaWitness::from(&shape, &ova_witness);

        // check it's still satisfied
        shape.is_relaxed_ova_satisfied::<PedersenCommitment<G>>(
            &relaxed_ova_instance,
            &relaxed_ova_witness,
            &pp,
        ).expect("satisfy failed");

        // fold multiple times with the instance/witness
        for _ in 0..3 {
            let (T, com_T) = commit_T(
                &shape,
                &pp[shape.num_vars..].to_vec(),
                &relaxed_ova_instance,
                &relaxed_ova_witness,
                &ova_instance,
                &ova_witness,
            ).unwrap();

            let r = Fr::rand(&mut thread_rng());
            relaxed_ova_instance = relaxed_ova_instance.fold(&ova_instance, &com_T, &r).unwrap();
            relaxed_ova_witness = relaxed_ova_witness.fold(&ova_witness, &T, &r).unwrap();
        }

        // now check it is still satisfied
        shape.is_relaxed_ova_satisfied::<PedersenCommitment<G>>(
            &relaxed_ova_instance,
            &relaxed_ova_witness,
            &pp,
        ).expect("satisfy failed");
    }
}
