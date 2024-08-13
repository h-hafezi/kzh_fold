use ark_bn254::{Bn254, Fq, Fr};
use ark_bn254::g1::Config as BNConfig;
use ark_grumpkin::GrumpkinConfig;

pub type E = Bn254;

pub type ScalarField = Fr;

pub type BaseField = Fq;

pub type G1 = BNConfig;

pub type G2 = GrumpkinConfig;
