////////////////////// borrowed from Haketon

use ark_ff::PrimeField;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::boolean::Boolean;
use ark_r1cs_std::select::CondSelectGadget;
use ark_r1cs_std::uint8::UInt8;
use ark_r1cs_std::{R1CSVar, ToBytesGadget};
use ark_relations::ns;
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use std::borrow::Borrow;

#[derive(Clone)]
pub struct Unsigned32Var<F: PrimeField> {
    pub bits: Vec<Boolean<F>>,
}

#[derive(Clone, Eq, PartialEq, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct Unsigned32 {
    pub bits: Vec<bool>,
}

impl<F: PrimeField> AllocVar<Unsigned32, F> for Unsigned32Var<F> {
    fn new_variable<T: Borrow<Unsigned32>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();
        let res = f();
        let entry = res.as_ref().map(|e| e.borrow()).map_err(|err| *err);
        let mut bits = Vec::new();
        for i in 0..entry.unwrap().bits.len() {
            let b = Boolean::new_variable(ns!(cs, "bit"), || entry.map(|e| e.bits[i]), mode)?;
            bits.push(b);
        }
        Ok(Unsigned32Var { bits })
    }
}

impl<F: PrimeField> R1CSVar<F> for Unsigned32Var<F> {
    type Value = Unsigned32;

    fn cs(&self) -> ConstraintSystemRef<F> {
        let cs = self.bits[0].cs();
        for i in 1..self.bits.len() {
            cs.clone().or(self.bits[i].cs());
        }
        cs.clone()
    }

    fn value(&self) -> Result<Self::Value, SynthesisError> {
        let mut bits = Vec::new();
        for b in &self.bits {
            bits.push(b.value()?);
        }
        Ok(Unsigned32 { bits })
    }
}

impl<F: PrimeField> ToBytesGadget<F> for Unsigned32Var<F> {
    fn to_bytes(&self) -> Result<Vec<UInt8<F>>, SynthesisError> {
        Ok(self.bits.chunks(8).map(UInt8::from_bits_le).collect())
    }
}

impl<F: PrimeField> Unsigned32Var<F> {
    pub fn new(cs: ConstraintSystemRef<F>) -> Self {
        Unsigned32Var {
            bits: vec![Boolean::new_witness(cs.clone(), || Ok(false)).unwrap(); 32],
        }
    }

    pub fn increment_inplace(&mut self) {
        let mut carry = Boolean::new_witness(self.cs().clone(), || Ok(true)).unwrap();
        for index in 0..self.bits.len() {
            let prev_bit = self.bits[index].clone();
            self.bits[index] =
                Boolean::conditionally_select(&carry, &self.bits[index].not(), &self.bits[index])
                    .unwrap();
            carry = Boolean::conditionally_select(&prev_bit, &carry, &Boolean::FALSE).unwrap();
        }
    }
}

impl<F: PrimeField> Unsigned32Var<F> {
    pub fn get_bits_canonical_order(&self, n: usize) -> Vec<Boolean<F>> {
        let res: Vec<Boolean<F>> = self.bits.iter().rev().cloned().collect();
        res[res.len() - n..].to_vec()
    }
}

impl Default for Unsigned32 {
    fn default() -> Self {
        Unsigned32 {
            bits: vec![false; 32],
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::math::Math;
    use crate::nexus_spartan::sparse_polynomial::uint32::{Unsigned32, Unsigned32Var};
    use ark_bls12_381::Fr;
    use ark_ff::PrimeField;
    use ark_r1cs_std::alloc::AllocVar;
    use ark_r1cs_std::R1CSVar;
    use ark_relations::r1cs::ConstraintSystem;
    use rand::random;

    impl Unsigned32 {
        pub fn representation(&self) -> u32 {
            self.bits.iter().enumerate().fold(
                0,
                |acc, (index, bit)| {
                    if *bit {
                        acc | (1 << index)
                    } else {
                        acc
                    }
                },
            )
        }

        // Increment the number by 1
        pub fn increment_inplace(&mut self) {
            let mut carry = true;
            for index in 0..self.bits.len() {
                if carry {
                    carry = self.bits[index];
                    // Flip the bit
                    self.bits[index] = !self.bits[index];
                } else {
                    break;
                }
            }
        }
    }

    impl<F: PrimeField> Unsigned32Var<F> {
        // Converts the U32 boolean vector into a native u32
        pub fn representation(&self) -> u32 {
            self.bits.iter().enumerate().fold(0, |acc, (index, bit)| {
                if bit.value().unwrap() {
                    acc | (1 << index)
                } else {
                    acc
                }
            })
        }
    }

    #[test]
    fn test_u32_increment() {
        let cs = ConstraintSystem::<Fr>::new_ref();
        let mut u = Unsigned32Var::<Fr>::new(cs.clone());
        for i in 1u32..40 {
            // increment u
            u.increment_inplace();
            // make sure u is in fact equal to i
            assert_eq!(u.representation(), i);
            // get the bit
            let res: Vec<bool> = u.get_bits_canonical_order(8)
                .into_iter().
                map(|bit| bit.value().unwrap())
                .collect();
            // make the sure the zk and non-zk versions are consistent
            assert_eq!(res, (i as usize).get_bits_canonical_order(8));
        }
        println!("{}", cs.num_constraints());
    }

    #[test]
    fn test_u32_new_variable() {
        let bits: Vec<bool> = (0..32).map(|_| random::<bool>()).collect();
        let mut u = Unsigned32 { bits };
        let cs = ConstraintSystem::<Fr>::new_ref();
        let mut v = Unsigned32Var::new_witness(cs.clone(), || Ok(u.clone())).unwrap();
        u.increment_inplace();
        v.increment_inplace();
        assert_eq!(u.representation(), v.representation());
    }
}