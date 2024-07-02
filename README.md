A $\sqrt{n}$ PCS

## Implementation TODO

What is implemented:
- Bivariate/univariate polynomial backend
- The polynomial commitment scheme
- A skeleton for the entire codebase (including tests and benchmarks)

I would start implementing things in the following order:
- Accumulation verifier
- Bivariate sumcheck
- Accumulation prover

I would start with the accumulation verifier because it's the main part that requires circuit work, and arkworks is
generally considered hard to write circuits. When that's done, the rest should be a breeze.

The codebase is sprinkled with `unimplemented!()` for things that should be implemented.

## Other things

There are a bunch of XXXs in the codebase. Some of them are important, others are really not important. I would suggest
to not spend time on them, except if they are actually blockers to continue.

Please treat the current codebase as a skeleton and feel free to adjust things to your liking.

## Building and running

Tests can be run with
```
cargo test --release
```

And benchmarks can be run with:
```
cargo bench
```

Currently only a PCS commitment test is implemented which can be run with `cargo bench --bench bench_pcs`.


