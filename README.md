Code Structure
--------------

.
├── accumulation
│   ├── Prover
│   ├── Verifier
│   └── Decider
├── accumulation_circuit
│   └── Accumulation Verifier Circuit
├── augmented_circuit
│   ├── Accumulation Circuit
│   ├── Spartan Partial Verifier Circuit
│   └── Matrix Evaluation Accumulator
├── gadgets
│   ├── non-native
│   │   └── NonNativeAffineVar
│   ├── r1cs
│   │   ├── R1CS Wrapper
│   │   ├── Relaxed R1CS
│   │   ├── Folding Methods
│   │   ├── OVA Structure
│   │   ├── Functions for OVA Folding
│   │   └── R1CSVar (for Nova)
│   ├── absorb
│   │   └── into_sponge_field_function
│   └── sparse
│       └── Sparse Matrix-Vector Evaluation Wrapper
├── halo_infinite
│   └── Halo Infinite Private Aggregation
├── hash
│   ├── Poseidon
│   ├── PoseidonVar
│   └── Pedersen
├── nexus_spartan
│   ├── matrix_evaluation_accumulator
│   │   ├── ABC Matrices Accumulator
│   │   └── Accumulator Verifier (Circuit)
│   ├── partial_verifier
│   │   └── Spartan Partial Verifier (Circuit)
│   ├── polycommitment
│   │   └── KZH PCS Wrapper
│   ├── sparse_polynomial
│   │   ├── Sparse Polynomial
│   │   └── Circuit Version
│   ├── sumcheck_circuit
│   │   └── Sumcheck (Circuit)
│   ├── unipoly
│   │   ├── Univariate Polynomial
│   │   ├── Circuit Functionalities
│   │   └── FFT-based Coefficient Interpolation
│   ├── commitments
│   ├── committed_relaxed_snark
│   ├── conversion
│   ├── ccr1cs
│   ├── crr1csproof
│   ├── errors
│   ├── product_tree
│   ├── r1csinstance
│   ├── sparse_mlpoly
│   └── timer
├── nova
│   ├── cycle_fold
│   │   ├── Secondary Circuit
│   │   └── OVA Folding Functions
│   ├── nova
│   │   ├── Nova Implementation
│   │   └── Benchmarking Utilities
│   └── util_test
├── pcs
│   ├── kzh3
│   │   └── KZH3 Implementation (Comparison)
│   └── multilienear_pcs
│       └── KZH PCS Implementation
├── polynomial
│   ├── eq_poly
│   │   ├── EqPolynomial
│   │   └── Circuit Version
│   ├── multilinear_poly
│   │   ├── Multilinear Polynomial
│   │   └── Circuit Version
│   └── univariate
│       ├── Univariate Polynomial
│       └── Circuit Version
├── signature_aggregation
│   ├── augmented_circuit
│   │   ├── Augmented Circuit
│   │   └── Mock Test (Constraint Counting)
│   ├── verifier_circuit
│   │   └── Signature Aggregation Verifier Circuit
│   └── signature_aggregation
│       └── Signature Aggregation Verifier/Prover
├── transcript
│   ├── transcript
│   │   └── Poseidon Backend
│   └── transcript_var
│       └── PoseidonVar Backend
├── lib
├── commitment
├── kzg
│   └── KZG Implementation (arkworks)
├── math
└── util
