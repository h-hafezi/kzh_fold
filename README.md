Code Structure
--------------

```
├── accumulatio: Prover/Verifier/Decider for accumulation outside of circuit
│   ├── accumualator.rs
│   ├── eq_tree.rs
│   └── mod.rs
├── accumulation_circuit: Circuit implementation of the accumulation verifier
│   ├── instance_circuit.rs
│   ├── mod.rs
│   ├── prover.rs
│   └── verifier_circuit.rs
├── augmented_circuit: Augmented circuit implementation for IVC and its test, it includes accumulation verifier circuit, Spartan partial verifier and matrix evaluation accumulator verifier
│   ├── augmented_circuit.rs
│   └── mod.rs
├── gadgets
│   ├── non-native
│   │   ├── mod.rs
│   │   ├── non_native_affine_var.rs
│   │   ├── readme.md
│   │   └── util.rs
│   ├── r1cs: It includes R1CS and Relaxed R1CS wrapper implementation, Folding methods for R1CS and ova structure, Functions for ova folding inside and outside of the circuit and R1CSVar instance used for Nova implementation
│   │   ├── conversion.rs
│   │   ├── mod.rs
│   │   ├── ova.rs
│   │   ├── r1cs.rs
│   │   └── r1cs_var.rs
│   ├── mod.rs
│   ├── absorb.rs
│   └── sparse
├── halo_infinite: Halo Infinite Private Aggregation
│   │   ├── errors.rs
│   │   ├── hpi.rs
│   │   ├── mod.rs
│   │   ├── pcs.rs
│   │   ├── private_aggregation.rs
│   │   └── prover.rs
├── hash: Poseidon/PoseidonVar wrappers, Pederson implementation
│   ├── mod.rs
│   ├── pederson.rs
│   └── poseidon.rs
├── nexus_spartan
│   ├── matrix_evaluation_accumulator: Implementation of ABC matrices accumulator + Accumulator verifier inside the circuit
│   │   ├── mod.rs
│   │   ├── prover.rs
│   │   └── verifier_circuit.rs
│   ├── partial_verifier: Spartan partial verifier implementation inside the circuit
│   │   ├── mod.rs
│   │   ├── partial_verifier.rs
│   │   └── partial_verifier_var.rs
│   ├── polycommitment: KZH wrapper used for Spartan
│   │   ├── error.rs
│   │   ├── kzh.rs
│   │   └── mod.rs
│   ├── sparse_polynomial
│   │   ├── mod.rs
│   │   ├── sparse_polynomial.rs
│   │   ├── sparse_polynomial_var.rs
│   │   ├── test.rs
│   │   └── uint32.rs
│   ├── sumcheck_circuit
│   │   ├── mod.rs
│   │   ├── sumcheck_circuit.rs
│   │   ├── sumcheck_circuit_var.rs
│   │   └── test.rs
│   ├── unipoly: Univariate polynomial and its circuit implementation, FFT-based coefficient interpolation used my matrix evaluation accumulator
│   │   ├── mod.rs
│   │   ├── unipoly.rs
│   │   └── unipoly_var.rs
│   ├── commitments.rs
│   ├── committed_relaxed_snark.rs
│   ├── conversion.rs
│   ├── ccr1cs.rs
│   ├── crr1csproof.rs
│   ├── errors.rs
│   ├── mod.rs
│   ├── product_tree.rs
│   ├── r1csinstance.rs
│   ├── sparse_mlpoly.rs
│   ├── sumcheck.rs
│   └── timer.rs
├── nova
│   ├── cycle_fold: Secondary circuit and folding functions for ova instances
│   │   ├── coprocessor.rs
│   │   ├── coprocessor_constraints.rs
│   │   ├── mod.rs
│   │   ├── readme.md
│   │   └── test.rs
│   ├── nova: Nova implementation using our backend, only used for benchmarking
│   │   ├── mod.rs
│   │   ├── prover.rs
│   │   ├── verifier_circuit.rs
│   │   └── verifier_circuit_var.rs
│   ├── mod.rs
│   └── util_test.rs
├── pcs
│   ├── kzh3: KZH3 implementation (comparison purpose only)
│   │   ├── mod.rs
│   │   └── kzh3.rs
│   ├── mod.rs
│   └── multilinear_pcs.rs
├── polynomial
│   ├── eq_poly: EqPolynomial implementation and circuit version, used by Spartan partial verifier circuit
│   │   ├── eq_poly.rs
│   │   ├── eq_poly_var.rs
│   │   └── mod.rs
│   ├── multilinear_poly: Multilinear polynomial implementation (similar to DensePolynomial in arkworks). Circuit version used for sumcheck implementation
│   │   ├── mod.rs
│   │   ├── multilinear_poly.rs
│   │   └── multilinear_poly_var.rs
│   └── univariate: Univariate polynomial and circuit version, circuit version used for sumcheck circuit
│   │   ├── mod.rs
│   │   ├── univariate_poly.rs
│   │   └── univariate_poly_var.rs
│   └── mod
├── signature_aggregation
│   ├── augmented_circuit: Augmented circuit for signature aggregation + Mock test for constraint counting
│   │   ├── ivc.rs
│   │   └── mod.rs
│   ├── verifier_circuit
│   │   ├── mod.rs
│   │   ├── prover.rs
│   │   ├── verifier_circuit.rs
│   │   └── verifier_circuit_var.rs
│   └── mod.rs
│   └── signature_aggregation.rs
├── transcript: Transcript implementation using Poseidon as backend
│   ├── mod.rs
│   ├── transcript.rs
│   └── transcript_var.rs
├── commitment.rs
├── lib.rs
├── kzg.rs ==> KZG implementation (borrowed from arkworks; used for comparison with KZH)
├── math.rs
└── util.rs
```