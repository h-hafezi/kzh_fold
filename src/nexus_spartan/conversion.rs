/*use crate::commitment::CommitmentScheme;
use crate::gadgets::r1cs::{R1CSInstance, R1CSShape, R1CSWitness};
use crate::gadgets::sparse::SparseMatrix;
use crate::nexus_spartan::crr1cs::{CRR1CSInstance, CRR1CSKey, CRR1CSShape, CRR1CSWitness};
use crate::nexus_spartan::errors::R1CSError;
use crate::nexus_spartan::polycommitments::PolyCommitmentScheme;
use crate::nexus_spartan::{Assignment, Instance};
use crate::polynomial::multilinear_poly::MultilinearPolynomial;
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ff::PrimeField;
use std::error::Error;
use std::fmt::Display;
use crate::constant_for_curves::G1Affine;

#[derive(Debug)]
pub enum ConversionError {
    ConversionError(R1CSError),
}

impl From<R1CSError> for ConversionError {
    fn from(error: R1CSError) -> Self {
        Self::ConversionError(error)
    }
}

impl Error for ConversionError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            ConversionError::ConversionError(e) => Some(e),
        }
    }
}

impl Display for ConversionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ConversionError(e) => write!(f, "Conversion error: {e}"),
        }
    }
}

impl<G> TryFrom<R1CSShape<G>> for CRR1CSShape<G::ScalarField>
where
    G: SWCurveConfig,
{
    type Error = ConversionError;
    fn try_from(shape: R1CSShape<G>) -> Result<Self, Self::Error> {
        let R1CSShape {
            num_constraints,
            num_vars,
            // This includes the leading `u` entry
            num_io,
            A,
            B,
            C,
        } = shape;
        // Spartan arranges the R1CS matrices using Z = [w, u, x], rather than [u, x, w]
        let rearrange =
            |matrix: SparseMatrix<G::ScalarField>| -> Vec<(usize, usize, G::ScalarField)> {
                matrix.iter().map(|(row, col, val)|
                    // this is a witness entry
                    if col >= num_io {
                        (row, col - num_io, val)
                    } else {
                        // this is an IO entry
                        (row, col + num_vars, val)
                    }).collect()
            };
        Ok(CRR1CSShape {
            inst: Instance::new(
                num_constraints,
                num_vars,
                // Spartan does not include the leading `u` entry in `num_inputs`.
                num_io - 1,
                rearrange(A).as_slice(),
                rearrange(B).as_slice(),
                rearrange(C).as_slice(),
            )?,
        })
    }
}

impl<PC, E> CRR1CSInstance<E, PC>
where
    E: Pairing,
    PC: PolyCommitmentScheme<E>,
{
    fn convert<G: SWCurveConfig, VC: CommitmentScheme<Projective<G>>>(
        instance: R1CSInstance<G, VC>,
        witness: R1CSWitness<G>,
        key: &PC::PolyCommitmentKey,
    ) -> Self where
        E: Pairing<G1Affine=Affine<G>, ScalarField=G::ScalarField>,
    {
        let R1CSInstance { X, .. } = instance;

        let poly_W = MultilinearPolynomial::new(witness.W);
        let comm_W = PC::commit(&poly_W, &key);

        CRR1CSInstance {
            input: Assignment::new(&X[1..]).unwrap(),
            comm_W,
        }
    }
}

impl<F: PrimeField> CRR1CSWitness<F> {
    fn convert<G: SWCurveConfig<ScalarField=F>>(witness: R1CSWitness<G>) -> Self {
        let R1CSWitness { W } = witness;
        CRR1CSWitness { W: Assignment::new(&W).unwrap() }
    }
}


#[cfg(test)]
mod tests {
    use ark_ec::CurveConfig;
    use crate::constant_for_curves::{ScalarField, E, G1};
    use crate::gadgets::r1cs::r1cs::{commit_T, R1CSInstance, R1CSWitness};
    use crate::hash::pederson::PedersenCommitment;
    use crate::nexus_spartan::crr1cs::{is_sat, CRR1CSInstance, CRR1CSKey, CRR1CSShape, CRR1CSWitness};
    use crate::nova::util_test::setup_test_r1cs;
    use ark_ec::short_weierstrass::Projective;
    // use ark_spartan::crr1csproof::{is_sat, CRR1CSKey};
    use ark_std::One;
    use rand::thread_rng;
    use crate::nexus_spartan::polycommitments::PolyCommitmentScheme;
    use crate::pcs::multilinear_pcs;
    use crate::pcs::multilinear_pcs::{PolyCommit, SRS};
    use crate::polynomial::multilinear_poly::MultilinearPolynomial;

    type Pedersen = PedersenCommitment<Projective<G1>>;

    #[test]
    fn test_conversion() {
        let (shape, U, W, pp) = setup_test_r1cs::<G1, Pedersen>(3, None, &());
        println!("{}", W.W.len());
        // check that the folded instance-witness pair is still satisfying
        shape.is_satisfied(&U, &W, &pp).expect("r1cs sat failed");

        // convert to the corresponding Spartan types
        let shape = CRR1CSShape::<<ark_bn254::g1::Config as CurveConfig>::ScalarField>::try_from(shape).unwrap();
        let min_num_vars = shape.get_num_vars();
        let SRS: SRS<E> = MultilinearPolynomial::setup(2, &mut thread_rng()).unwrap();
        let key: CRR1CSKey<E, MultilinearPolynomial<ScalarField>> = CRR1CSKey::new(&SRS, shape.get_num_cons(), shape.get_num_vars());
        let instance: CRR1CSInstance<E, MultilinearPolynomial<ScalarField>> = CRR1CSInstance::convert(U, W.clone(), &key.keys.ck);
        let witness = CRR1CSWitness::<ScalarField>::convert(W);


        // check that the Spartan instance-witness pair is still satisfying
        // assert!(is_sat(&shape, &instance, &witness, &key).unwrap());
    }
}
 */