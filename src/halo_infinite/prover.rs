/*#![allow(unused)]
use ark_ec::pairing::Pairing;
use ark_std::marker::PhantomData;

/// The polynomial f(x) behind the commiment C, evaluates to f(x) = y
pub struct Instance<E: Pairing> {
    pub C: E::G1Affine,
    pub x: E::ScalarField,
    pub y: E::ScalarField,
}

pub struct PublicAggregationScheme<E: Pairing> {
    _field: PhantomData<E>
}

impl<E: Pairing> PublicAggregationScheme<E> {
    pub fn aggregate(instance_1: Instance<E>, instance_2: Instance<E>) -> () {
        ()
    }
}

#[cfg(test)]
mod tests {
    #[test]
    pub fn aggregate() {
    }
}


 */