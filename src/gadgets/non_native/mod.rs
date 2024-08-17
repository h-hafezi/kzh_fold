use std::{borrow::Borrow, marker::PhantomData};

use ark_ff::{Field, PrimeField};
use ark_r1cs_std::{
    alloc::{AllocVar, AllocationMode},
    fields::{
        fp::FpVar,
        nonnative::{params::OptimizationType, AllocatedNonNativeFieldVar, NonNativeFieldVar},
    },
    ToBytesGadget, ToConstraintFieldGadget,
};
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, OptimizationGoal, SynthesisError};

pub mod non_native_affine_var;
pub mod util;
