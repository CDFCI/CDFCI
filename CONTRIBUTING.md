# Contributing to CDFCI

Thank you for your interest in contributing to **CDFCI** (Coordinate Descent Full Configuration Interaction)! We welcome all kinds of contributions including code, bug reports, feature requests, and documentation.

---

## 📌 Quick Start

If you're new to the project:

1. **Fork the repository** and clone it to your local machine.
2. **Create a new branch** from `develop`:
   `git checkout -b feature/my-feature develop`
3. Make your changes.
4. **Push to your fork**:
   `git push origin feature/my-feature`
5. **Submit a Pull Request** (PR) to the `develop` branch.

For structural details and codebase organization, see [DEVELOPING.md](./DEVELOPING.md).

---

## 🌿 Branching Strategy

We follow a Git Flow-like model:

* `master`: Stable release-only. Code here is tested and version-tagged.
* `develop`: Main integration branch.
* `feature/<name>`: New features (branched from `develop`).

Please avoid pushing directly to `develop` or `master`.

---

## 🧪 Testing & Validation

* Run unit tests before PR submission.
* If adding a new feature, consider adding a test case under `regression_tests/`.
* For major changes, include example JSON inputs in `examples/`.

CI runs on every PR to ensure build and test correctness.

---

## ✨ Code Style

We use:

* **C++17 standard**
* **4 spaces for indentation**
* Clang-format config consistent with Visual Studio style (see `DEVELOPING.md`)

Auto-format your code with:

```bash
clang-format -i src/**/*.h src/**/*.cpp
```

---

## 📞 Communication

* For major contributions, open an Issue first or email [Prof. Yingzhou Li](mailto:yingzhouli@fudan.edu.cn).
* For general bug reports or feature suggestions, please use [GitHub Issues](https://github.com/CDFCI/CDFCI-private/issues).

We appreciate your time and effort to help improve CDFCI!
