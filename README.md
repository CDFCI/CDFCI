# CDFCI

[![CDFCI\_Test](https://github.com/CDFCI/CDFCI/actions/workflows/cdfci_test.yml/badge.svg?branch=master)](https://github.com/CDFCI/CDFCI/actions/workflows/cdfci_test.yml)

CDFCI (Coordinate Descent Full Configuration Interaction) is a modern, efficient solver for the ground-state electronic structure problem in quantum chemistry. It is implemented in modern C++17 with support for large-scale sparse wavefunctions and adaptive compression strategies.

📖 中文说明请参见：[README\_zh.md](README_zh.md)

---

## ✨ Overview

CDFCI reformulates the full configuration interaction (FCI) eigenvalue problem as a large-scale unconstrained nonconvex optimization problem. It solves this problem using adaptive coordinate descent with deterministic compression.

Key features:

* Coordinate-wise update of Slater determinants
* Adaptive compression via `z_threshold`
* Deterministic and memory-aware wavefunction storage
* Optional support for OpenMP (multi-threaded)

---

## 📄 Citation

CDFCI implements algorithms from the following papers,

* **CDFCI Core Algorithm**
  Z. Wang, Y. Li, J. Lu, [*JCTC*, 15(6), 2019](https://pubs.acs.org/doi/10.1021/acs.jctc.9b00138)

* **Optimal Orbital Selection (OptOrbFCI)**
  Y. Li, J. Lu, [*JCTC*, 16(10), 2020](https://pubs.acs.org/doi/10.1021/acs.jctc.0c00613)

* **Low-lying Excited States Extension**
  Z. Wang, Z. Zhang, J. Lu, Y.Li, [*JCTC*, 19(21), 2023](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.3c00452)

* **Multicoordinate Parallel Extension**
  Y. Zhang, W. Gao, Y. Li, [*JCTC*, 21(5), 2025](https://pubs.acs.org/doi/10.1021/acs.jctc.4c01530)

with convergence guarantee:

* **Coordinate Descent Convergence Theory**
  Y. Li, J. Lu, Z. Wang, [*SIAM J. Sci. Comput.*, 41(4), 2019](https://doi.org/10.1137/18M1202505)

Please cite the relevant paper(s) if you use this software.

---

## 🚀 Getting Started

### Clone and Build

```bash
git clone --recursive https://github.com/CDFCI/CDFCI.git
cd CDFCI
make
```

The following executables will be generated in the `bin` directory:

* `cdfci`: single-threaded version
* `cdfci_omp`: OpenMP-enabled version
* `xcdfci`: for excited states
* `optorbfci`: for orbital rotation and compression
* `cdfci_tools`: auxiliary utilities

Alternatively, you may use `make cdfci `, `make optorbfci` or  `make xcdfci` to generate specific executable.

Useful release-related targets:

* `make check`: run test gate (`make test` alias)
* `make install PREFIX=/path`: install binaries to `${PREFIX}/bin` (default `/usr/local/bin`)
* `make uninstall PREFIX=/path`: remove installed binaries
* `make release-check`: run clean build + regression tests + examples

### Requirements

* C++17 compiler (GCC ≥ 9, Clang ≥ 10, or ICC)
* [Eigen3](https://eigen.tuxfamily.org) (will be clone recursively)
* OpenMP (optional)
---

## 📥 Input Format (JSON)

The solver accepts a JSON input file. Example:

```json
{
    "hamiltonian": {
        "type": "molecule",
        "molecule": {
            "fcidump_path": "/path/to/FCIDUMP",
            "threshold": 0.0
        }
    },
    "solver":{
        "type": "cdfci",
        "cdfci": {
            "num_iterations": 150000,
            "report_interval": 10000,
            "z_threshold": 0,
            "z_threshold_search": false
        }
    },
    "max_memory": 0.005
}
```

More examples, including [a full documentation for input and output](./examples/README.md), can be found in the `examples/` directory.

---

## 🧪 Testing

To compile and run regression tests:

```bash
make test
```

Example systems and expected energies are provided in `regression_tests`.

---

## 📊 Reproducibility

The `examples/` and `papers/` directories contain all scripts and input files to reproduce published results. Requires Python and [PySCF](https://github.com/pyscf/pyscf) to generate FCIDUMP files.

---

## 👥 Developers

* Zhe Wang (Mathematics Department, Duke University)
* [Jianfeng Lu](https://services.math.duke.edu/~jianfeng/) (Mathematics Department, Duke University)
* [Yingzhou Li](http://yingzhouli.com/) (School of Mathematical Sciences, Fudan University)
* [Yuejia Zhang](https://ninotreve.github.io/) (School of Mathematical Sciences, Fudan University)

---

## 🪙 Acknowledgements

This work was supported in part by the U.S. National Science Foundation (NSF) under Grant Nos. DMS-1454939 and DMS-2012286; by the National Natural Science Foundation of China under Grant Nos. 12271109 and 12526211; by the Shanghai Pilot Program for Basic Research-Fudan University under Grant No. 21TQ1400100 (22TQ017); and by the Scientific Research Innovation Capability Support Project for Young Faculty under Grant No. SRICSPYF-ZY2025159.

---

## ⚖️ License

This project is licensed under the BSD 3-Clause License. See the [LICENSE](./LICENSE) file for details.
