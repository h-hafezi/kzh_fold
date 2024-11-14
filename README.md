Code Structure
--------------

```
├── accumulation: KZH accumulation implementation (out of circuit)
├── accumulation_circuit: Circuit implementation of the KZH accumulation verifier
├── augmented_circuit: Augmented circuit implementation for KZH-based folding it includes KZH acc verifier circuit, Spartan partial verifier and matrix evaluation accumulator verifier
│   ├── r1cs: It includes R1CS and Relaxed R1CS wrapper implementation, Folding methods for R1CS and ova structure, Functions for ova folding inside and outside of the circuit and R1CSVar instance used for Nova implementation
├── halo_infinite: Halo Infinite Private Aggregation
├── hash: Poseidon implementation (in circuit and out of circuit)
├── nexus_spartan: Spartan implementation by Nexus (plus our own modifications as indicated below)
│   ├── matrix_evaluation_accumulator: Implementation of ABC matrices accumulation
│   ├── partial_verifier: Spartan partial verifier implementation inside the circuit
│   ├── polycommitment: KZH wrapper used for Spartan
│   ├── sumcheck_circuit: Sumcheck verification circuit
├── nova: Nova implementation for benchmarking purposes
│   ├── cycle_fold: Secondary circuit and folding functions for ova instances
├── pcs
│   ├── kzh3: KZH3 implementation
├── signature_aggregation
│   ├── augmented_circuit: Augmented circuit for signature aggregation + Mock test for constraint counting
│   ├── verifier_circuit: Verifier circuit for signature aggregation protocol
├── kzg.rs ==> KZG implementation (borrowed from arkworks; used for comparison with KZH)
```
