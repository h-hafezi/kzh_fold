use crate::commitment::{Commitment, CommitmentScheme};
use crate::gadgets::non_native::util::convert_field_one_to_field_two;
use crate::gadgets::r1cs::r1cs_var::{R1CSInstanceVar, RelaxedR1CSInstanceVar};
use crate::gadgets::r1cs::{R1CSInstance, RelaxedR1CSInstance};
use ark_crypto_primitives::sponge::constraints::AbsorbGadget;
use ark_ec::short_weierstrass::{Projective, SWCurveConfig};
use ark_ec::{AffineRepr, CurveConfig};
use ark_ff::PrimeField;
use ark_r1cs_std::fields::fp::FpVar;

pub fn r1cs_instance_to_sponge_vector<G: SWCurveConfig, C: CommitmentScheme<Projective<G>>>(instance: &R1CSInstance<G, C>) -> Vec<G::ScalarField>
where
    <G as CurveConfig>::BaseField: PrimeField,
{
    let mut res = Vec::<G::ScalarField>::new();

    // append vector X
    res.extend(instance.X.clone());

    // convert group into native field elements
    let w = instance.commitment_W.into_affine().xy().unwrap();
    let x = convert_field_one_to_field_two::<G::BaseField, G::ScalarField>(w.0);
    let y = convert_field_one_to_field_two::<G::BaseField, G::ScalarField>(w.1);
    res.extend(vec![x, y]);

    res
}

pub fn r1cs_instance_var_to_sponge_vector<F, G1, C1>(instance: &R1CSInstanceVar<G1, C1>) -> Vec<FpVar<F>>
where
    G1: SWCurveConfig<ScalarField=F>,
    G1::BaseField: PrimeField,
    F: ark_ff::PrimeField,
{
    let mut res = Vec::<FpVar<F>>::new();

    // append vector X
    res.extend(instance.X.clone());
    res.extend(instance.commitment_W.to_sponge_field_elements().unwrap());

    res
}

pub fn relaxed_r1cs_instance_to_sponge_vector<G: SWCurveConfig, C: CommitmentScheme<Projective<G>>>(instance: &RelaxedR1CSInstance<G, C>) -> Vec<G::ScalarField>
where
    <G as CurveConfig>::BaseField: ark_ff::PrimeField,
{
    let mut res = Vec::<G::ScalarField>::new();

    // append vector X
    res.extend(instance.X.clone());

    // convert group into native field elements
    let w = instance.commitment_W.into_affine().xy().unwrap();
    let x = convert_field_one_to_field_two::<G::BaseField, G::ScalarField>(w.0);
    let y = convert_field_one_to_field_two::<G::BaseField, G::ScalarField>(w.1);
    res.extend(vec![x, y]);

    // convert group into native field elements
    let e = instance.commitment_E.into_affine().xy().unwrap();
    let x = convert_field_one_to_field_two::<G::BaseField, G::ScalarField>(e.0);
    let y = convert_field_one_to_field_two::<G::BaseField, G::ScalarField>(e.1);
    res.extend(vec![x, y]);

    res
}

pub fn relaxed_r1cs_instance_var_to_sponge_vector<F, G1, C1>(instance: &RelaxedR1CSInstanceVar<G1, C1>) -> Vec<FpVar<F>>
where
    G1: SWCurveConfig<ScalarField=F>,
    G1::BaseField: PrimeField,
    F: ark_ff::PrimeField,
{
    let mut res = Vec::<FpVar<F>>::new();

    // append vector X
    res.extend(instance.X.clone());
    res.extend(instance.commitment_W.to_sponge_field_elements().unwrap());
    res.extend(instance.commitment_E.to_sponge_field_elements().unwrap());

    res
}
