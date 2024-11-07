use crate::commitment::{Commitment, CommitmentScheme};
use crate::gadgets::non_native::util::convert_field_one_to_field_two;
use crate::gadgets::r1cs::r1cs_var::{R1CSInstanceVar, RelaxedR1CSInstanceVar};
use crate::gadgets::r1cs::{R1CSInstance, RelaxedR1CSInstance};
use ark_crypto_primitives::sponge::constraints::AbsorbGadget;
use ark_ec::short_weierstrass::{Projective, SWCurveConfig};
use ark_ec::{AffineRepr};
use ark_ff::PrimeField;
use ark_r1cs_std::fields::fp::FpVar;

// Functions below are used by Nova, since in Nova we require to hash R1CS instances

impl<G: SWCurveConfig, C: CommitmentScheme<Projective<G>>> R1CSInstance<G, C>
where
    G::BaseField: PrimeField
{
    /// This function should be consistent with the counterpart for R1CSInstanceVar as below
    pub fn to_sponge_field_elements(&self) -> Vec<G::ScalarField> {
        let mut res = Vec::<G::ScalarField>::new();

        // append vector X
        res.extend(self.X.as_slice());

        // convert group into native field elements
        let w = self.commitment_W.into_affine().xy().unwrap();
        let x = convert_field_one_to_field_two::<G::BaseField, G::ScalarField>(w.0);
        let y = convert_field_one_to_field_two::<G::BaseField, G::ScalarField>(w.1);
        res.extend(vec![x, y]);

        res
    }
}

impl<F, G1, C1> R1CSInstanceVar<G1, C1>
where
    G1: SWCurveConfig<ScalarField = F>,
    G1::BaseField: PrimeField,
    F: PrimeField
{
    pub fn to_sponge_field_elements(&self) -> Vec<FpVar<F>> {
        let mut res = Vec::<FpVar<F>>::new();

        // append vector X
        res.extend(self.X.clone());
        // function commitment_W.to_sponge_field_elements() is consistent with the output of
        // convert_field_one_to_field_two::<G::BaseField, G::ScalarField>()
        res.extend(self.commitment_W.to_sponge_field_elements().unwrap());

        res
    }
}

impl<G: SWCurveConfig, C: CommitmentScheme<Projective<G>>> RelaxedR1CSInstance<G, C>
where
    G::BaseField: PrimeField
{
    /// This function should be consistent with the counterpart for RelaxedR1CSInstanceVar as below
    pub fn to_sponge_field_elements(&self) -> Vec<G::ScalarField> {
        let mut res = Vec::<G::ScalarField>::new();

        // append vector X
        res.extend(self.X.as_slice());

        // convert group into native field elements
        let w = self.commitment_W.into_affine().xy().unwrap();
        let x = convert_field_one_to_field_two::<G::BaseField, G::ScalarField>(w.0);
        let y = convert_field_one_to_field_two::<G::BaseField, G::ScalarField>(w.1);
        res.extend(vec![x, y]);

        // convert group into native field elements
        let e = self.commitment_E.into_affine().xy().unwrap();
        let x = convert_field_one_to_field_two::<G::BaseField, G::ScalarField>(e.0);
        let y = convert_field_one_to_field_two::<G::BaseField, G::ScalarField>(e.1);
        res.extend(vec![x, y]);

        res
    }
}

impl<F, G1, C1> RelaxedR1CSInstanceVar<G1, C1>
where
    G1: SWCurveConfig<ScalarField = F>,
    G1::BaseField: PrimeField,
    F: PrimeField
{
    pub fn to_sponge_field_elements(&self) -> Vec<FpVar<F>> {
        let mut res = Vec::<FpVar<F>>::new();

        // append vector X
        res.extend(self.X.clone());
        res.extend(self.commitment_W.to_sponge_field_elements().unwrap());
        res.extend(self.commitment_E.to_sponge_field_elements().unwrap());

        res
    }
}
