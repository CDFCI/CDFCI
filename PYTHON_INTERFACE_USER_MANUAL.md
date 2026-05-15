# CDFCI Python Interface User Manual

This manual describes how to build and use the Python interface for CDFCI, including:

- Running CDFCI/XCDFCI/OptOrbFCI from Python options
- Running CDFCI directly from RHF integrals in memory (no FCIDUMP intermediate file)
- Integrating CDFCI with PySCF SCF/CASCI workflows
- Accessing RDMs through the PySCF-compatible solver adapter

## 1. Prerequisites

- Linux/macOS with a C++17 compiler
- Python 3.8+
- `numpy`
- Optional for integration: `pyscf`

Recommended environment setup:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install numpy pyscf
```

## 2. Build the Python Extension

From repository root:

```bash
make python-module
```

This builds the extension module under `python/_cdfci<ext-suffix>`.

## 3. Import Setup

When running scripts from repository root, use:

```bash
PYTHONPATH=python python your_script.py
```

Or inside Python:

```python
import sys
sys.path.insert(0, "python")
import cdfci
```

## 4. Public Python API

The package exports the following main symbols:

- `CDFCI`
- `CDFCIResult`
- `run_cdfci(option)`
- `run_xcdfci(option)`
- `run_optorbfci(option)`
- `run_cdfci_rhf_integrals(option, h1e, eri, norb, nelec, core_energy=0.0, ms2=0)`
- `tools_frozen_orbital(option)`
- `tools_symmetry_connection(option)`
- `run_solver(option)`
- `run_from_pyscf_scf(mf, option=None, solver="cdfci", fcidump_path=None, keep_fcidump=False)`
- `CDFCIPySCFFCISolver` and `as_pyscf_fcisolver(...)` (available when PySCF is installed)

### Option Input Types

For APIs that accept `option`, you may provide:

- A Python `dict`
- A JSON string
- A path to a JSON file

## 5. Quick Start: Run CDFCI from FCIDUMP via Options

```python
import cdfci

option = {
    "hamiltonian": {
        "type": "molecule",
        "molecule": {
            "fcidump_path": "data/h2o_sto3g_psi4.FCIDUMP",
            "threshold": 0.0,
        },
    },
    "solver": {
        "type": "cdfci",
        "cdfci": {
            "num_iterations": 150000,
            "report_interval": 10000,
            "z_threshold": 0.0,
            "z_threshold_search": False,
            "stopping_dx_threshold": 1e-10,
        },
    },
    "max_memory": 1.0,
}

res = cdfci.run_cdfci(option)
print(res.energy, res.iterations)
```

## 6. Run XCDFCI (Excited States)

```python
import cdfci

option = {
    "hamiltonian": {
        "type": "molecule",
        "molecule": {
            "fcidump_path": "data/h2o_sto3g_psi4.FCIDUMP",
            "threshold": 0.0,
        },
    },
    "solver": {
        "type": "xcdfci",
        "xcdfci": {
            "num_iterations": 150000,
            "report_interval": 10000,
        },
    },
    "num_states": 2,
    "max_memory": 1.0,
}

res = cdfci.run_xcdfci(option)
print("ground:", res.energy)
print("states:", res.state_energies)
```

## 7. In-Memory RHF Integral Path (No FCIDUMP File)

Use this when you already have MO-basis RHF integrals:

- `h1e` shape: `(norb, norb)`
- `eri` shape: `(norb, norb, norb, norb)`

```python
import cdfci
import numpy as np

norb = 4
nelec = 4
h1e = np.zeros((norb, norb))
eri = np.zeros((norb, norb, norb, norb))

option = {
    "solver": {
        "type": "cdfci",
        "cdfci": {
            "num_iterations": 50000,
            "report_interval": 5000,
        },
    },
    "max_memory": 1.0,
}

res = cdfci.run_cdfci_rhf_integrals(
    option=option,
    h1e=h1e,
    eri=eri,
    norb=norb,
    nelec=nelec,
    core_energy=0.0,
    ms2=0,
)
print(res.energy)
```

## 8. PySCF SCF -> CDFCI Workflow

`run_from_pyscf_scf` supports direct usage from a converged PySCF SCF object.

For RHF + `solver="cdfci"`, it automatically uses the in-memory integral path.
For other cases (for example current XCDFCI flow), it may use FCIDUMP fallback.

```python
import cdfci
from pyscf import gto, scf

mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="sto-3g", spin=0)
mf = scf.RHF(mol)
mf.kernel()

option = {
    "solver": {
        "type": "cdfci",
        "cdfci": {
            "num_iterations": 50000,
            "report_interval": 5000,
        },
    },
    "max_memory": 1.0,
}

res = cdfci.run_from_pyscf_scf(mf, option=option, solver="cdfci")
print(res.energy)
```

## 9. PySCF CASCI/CASSCF Adapter

Use `as_pyscf_fcisolver` to plug CDFCI into PySCF CASCI/CASSCF.

```python
import cdfci
from pyscf import gto, scf, mcscf

mol = gto.M(
    atom="N 0 0 0; N 0 0 1.12079733",
    basis="cc-pvdz",
    spin=0,
)
mf = scf.RHF(mol).run()

option = {
    "solver": {
        "type": "cdfci",
        "cdfci": {
            "num_iterations": 500000,
            "report_interval": 10000,
        },
    },
    "max_memory": 1.0,
}

mc = mcscf.CASCI(mf, ncas=10, nelecas=10)
mc.fcisolver = cdfci.as_pyscf_fcisolver(option=option, solver="cdfci")
e_tot = mc.kernel()[0]
print(e_tot)
```

## 10. RDM Support in Adapter

`CDFCIPySCFFCISolver` provides:

- `make_rdm1(ci, ncas, nelec)`
- `make_rdm12(ci, ncas, nelec)`

Internally, the adapter enables CDFCI RDM output and loads spin-summed `OneRDM.out` and `TwoRDM.out` into NumPy arrays for PySCF consumption.

## 11. Legacy Low-Level Class

The `CDFCI` class is retained for simple direct FCIDUMP-driven workflows:

```python
from cdfci import CDFCI

solver = CDFCI("data/h2o_sto3g_psi4.FCIDUMP")
solver.set_num_iterations(100000)
solver.set_report_interval(5000)
solver.set_max_memory(1.0)
solver.set_max_load_factor(0.79)
res = solver.run()
print(res.energy)
```

## 12. Advanced Solver Options

All solver knobs are passed through the option dictionary and follow the same schema as the C++ JSON input.

Common fields include (depending on solver type):

- `num_iterations`
- `report_interval`
- `z_threshold`
- `z_threshold_search`
- `stopping_dx_threshold`
- `max_wavefunction_size`
- `coordinate_pick`
- `coordinate_update`

For complete option details, see:

- `examples/README.md`
- input JSON files under `examples/` and `papers/`

## 13. Typical Scripting Patterns

### Pattern A: JSON-on-disk configuration

- Keep solver configurations in version-controlled JSON files.
- Pass the JSON file path directly to `run_cdfci`/`run_xcdfci`.

### Pattern B: Programmatic option patching

- Start from a base Python dict.
- Patch runtime-sensitive parameters (`num_iterations`, memory, tolerances) in scripts.

### Pattern C: PySCF-first workflows

- Build molecule and SCF in PySCF.
- Call `run_from_pyscf_scf` for quick CDFCI runs.
- Use `as_pyscf_fcisolver` for CASCI/CASSCF integration.

## 14. Troubleshooting

### ImportError: No module named `_cdfci`

- Rebuild extension: `make python-module`
- Ensure `PYTHONPATH=python` is set when running scripts.

### ImportError: PySCF is required ...

- Install PySCF in the active environment: `pip install pyscf`

### ValueError: invalid h1e/eri shape

- Ensure `h1e.shape == (norb, norb)`
- Ensure `eri.shape == (norb, norb, norb, norb)` for `run_cdfci_rhf_integrals`

### Convergence too slow

- Increase `num_iterations`
- Adjust `report_interval` and stopping thresholds
- Review active space or basis size for PySCF-coupled jobs

## 15. Reproducible Demo Script

A full benchmark-style example is provided in:

- `python/test.py`

Run with:

```bash
PYTHONPATH=python python3 python/test.py
```

---

If you cite CDFCI in publications, please also reference the project README and related papers under `papers/`.
