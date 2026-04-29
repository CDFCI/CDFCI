# Test

CDFCI uses a regression-style test gate driven by doctest and input sets in
`regression_tests/`.

## What gets executed

* `make build_test` builds:
  * `test/cdfci_test` (single-thread)
  * `test/cdfci_test_omp` (OpenMP)
* `make test` runs both executables.

The main test entry is `test/test_system.cpp`, which loads each case from
`regression_tests/<system>/input.json`, runs the solver, and compares the final
energy against `result.txt` tolerance.

## Regression case layout

Each case directory under `regression_tests/` typically contains:

* `input.json`: CDFCI input configuration.
* `result.txt`: reference energy and allowed tolerance.
* `output` (optional): sample output snapshot.
* `README.md` (optional): case-specific notes.

## Notes

* These are regression tests for stability and reproducibility, not strict
  full-convergence proofs for every run length.
* Numerical tolerance is intentionally case-dependent.
