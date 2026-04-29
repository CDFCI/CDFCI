# Developing Guide for CDFCI

This guide provides a technical overview for developers working on the CDFCI codebase.

---

## 🔧 Compiler & Standard

* **C++17** is the required standard.
* Cross-platform compatibility is maintained (Linux, macOS).
* Main compiler targets: `g++ >= 9`, `clang++ >= 10`

---

## 🗂 Directory Structure

```
CDFCI/
├── bin/                    # Compiled executables (ignored by Git)
├── data/                   # Sample FCIDUMP data files
├── doc/                    # Project documentation
├── examples/               # Example input JSONs and configurations
│   ├── demo_input_*.json
│   ├── frozen_core.json
│   └── symm_conn.json
├── papers/                 # Benchmark scripts for publication
│   ├── 2019-cdfci/
│   ├── 2020-optorbfci/
│   ├── 2023-xcdfci/
│   └── 2025-mcdfci/
├── regression_tests/       # Functional tests (small molecules)
├── src/                    # Main source code
│   ├── lib/                # External dependencies (optional)
│   └── main.cpp
├── test/                   # Unit tests
├── Makefile                # Build system (customized)
├── LICENSE
└── README.md
```

---

## 🔄 Build Instructions

### Dependencies:

* Eigen
* OpenMP (optional but recommended)

### Build

```bash
make clean && make -j
```

Optional targets:

* `make check` — test gate alias
* `make test` — run regression tests
* `make install PREFIX=/path` — install binaries to `${PREFIX}/bin`
* `make uninstall PREFIX=/path` — uninstall binaries
* `make release-check` — clean build + tests + examples
* `make format` — run clang-format

---

## 🌱 Coding Style

We follow Visual Studio / LLVM-style formatting. See `.clang-format`:

```json
{
    "BasedOnStyle": "LLVM",
    "IndentWidth": 4,
    "BreakBeforeBraces": "Allman",
    "ColumnLimit": 0
}
```

Please run:

```bash
make format
```

---

## 🧪 Testing Strategy

* Minimal working tests in `regression_tests/`
* Target molecules: H2O, C2, N2, etc.
* Validate correctness by comparing with known energies from literature

---

## 🧭 Version Control Best Practices

* Work on branches from `develop`
* Keep commits small and meaningful
* Sync with upstream regularly

See [CONTRIBUTING.md](./CONTRIBUTING.md) for pull request workflow.

---

For questions, contact [yingzhouli@fudan.edu.cn](mailto:yingzhouli@fudan.edu.cn) or [yuejiazhang21@m.fudan.edu.cn](yuejiazhang21@m.fudan.edu.cn).

