use ark_bn254::{Bn254, Fq, Fr};
use ark_bn254::g1::Config as BNConfig;
use ark_grumpkin::GrumpkinConfig;
use ark_bn254::g1::{G1Affine as g1};
use ark_bn254::g2::{G2Affine as g2};

/// Since we use the cycle of curves (Bn254, Grumpkin) throughout our tests, we define some types here, so later we can easily use them in out tests

pub type E = Bn254;

pub type ScalarField = Fr;

pub type BaseField = Fq;

pub type G1Affine = g1;

pub type G2Affine = g2;

pub type G1 = BNConfig;

pub type G2 = GrumpkinConfig;
