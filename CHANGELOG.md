# Changelog

All notable changes to **CDFCI** will be documented in this file.
This project follows [Semantic Versioning](https://semver.org/) and Git Flow branching model.

---

## \[v0.9.0-pre] - 2025-09-22

### Highlights

* First feature-complete implementation of CDFCI and its variants (OptOrbFCI, xCDFCI, mCDFCI)
* Provides reproducible workflows, regression tests, and CI integration.
* Pre-release version intended for collaborator testing and feedback.

### Added

* Solvers:
    * `cdfci` (single-threaded) and `cdfci_omp` (OpenMP parallel)
    * `optorbfci` (single-threaded)
    * `xcdfci` (single-threaded)
* Auxiliary executable: `cdfci_tools`
* JSON-based configuration interface
* Regression test framework and example systems
* Reproducibility scripts for benchmark datasets
* GitHub Actions workflow for CI testing
* CONTRIBUTING and DEVELOPING guides

### Changed

* Project restructured with dedicated  `examples`, `regression_tests` and `papers` directories
* New `README` and `README_zh` with algorithm overview and usage guide
* Legacy `develop/xxx` branches migrated to `feature/xxx`

### License

* BSD 3-Clause License

### Limitations

* Convergence depends strongly on the choice of initial state
* Python interface not yet available

### Missing Features

* Python interface
---

## [v1.0.0] - 2026-04-29

### Highlights

* First stable public release of the CDFCI C++ solver family.
* Regression tests pass on current `develop` branch release candidate.
* Release pipeline and packaging workflow hardened for reproducible delivery.

### Recap from v0.9.0-pre

The following feature-complete foundations from `v0.9.0-pre` remain part of `v1.0.0`:

* Solver executables: `cdfci`, `cdfci_omp`, `xcdfci`, `optorbfci`, and `cdfci_tools`.
* JSON-based configuration interface for Hamiltonian/solver options.
* Regression test framework with reference systems in `regression_tests/`.
* Reproducibility inputs and scripts in `examples/` and `papers/`.
* Developer and contributor documentation (`README`, `README_zh`, `DEVELOPING`, `CONTRIBUTING`).

### Added

* Build/release targets in `Makefile`:
    * `check` (alias of test gate)
    * `install` / `uninstall` for system deployment
    * `release-check` for end-to-end pre-release validation
* CI trigger expansion:
    * `develop` + `master` branch validation
    * tag trigger `v*` for release builds
* Release artifact workflow in GitHub Actions:
    * build binaries
    * package tarball with `bin/`, `README.md`, `CHANGELOG.md`, `LICENSE`
    * upload versioned release artifact

### Changed

* `clean` target now removes actual generated binaries under `bin/` and debug/test outputs.
* Changelog status aligned with implementation (checkpoint/restart support exists in solver options).
* Release narrative updated from pre-release planning to stable release record.

### Validation

* `make test` passes for regression suite.
* Release targets validated via dry-run (`make -n install`, `make -n uninstall`, `make -n clean`).

### Known Limitations

* Python interface is intentionally not included in v1.0.0.

### License

* BSD 3-Clause License