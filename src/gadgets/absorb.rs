use crate::commitment::{Commitment, CommitmentScheme};
use crate::gadgets::r1cs::{R1CSInstance, RelaxedR1CSInstance};
use ark_ec::{AffineRepr, CurveConfig, CurveGroup};
use ark_ec::short_weierstrass::{Projective, SWCurveConfig};
use crate::constant_for_curves::E;
use crate::gadgets::non_native::util::{convert_affine_to_scalars, convert_field_one_to_field_two};

pub fn r1cs_instance_to_sponge_vector<G: SWCurveConfig, C: CommitmentScheme<Projective<G>>>(instance: &R1CSInstance<G, C>) -> Vec<G::ScalarField>
where <G as CurveConfig>::BaseField: ark_ff::PrimeField
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

pub fn relaxed_r1cs_instance_to_sponge_vector<G: SWCurveConfig, C: CommitmentScheme<Projective<G>>>(instance: &RelaxedR1CSInstance<G, C>) -> Vec<G::ScalarField>
where <G as CurveConfig>::BaseField: ark_ff::PrimeField
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