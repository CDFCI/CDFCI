# CDFCI Examples Documentations

This README focuses on a fast path to try the solvers and understand typical **inputs/outputs**.

## 1) What’s included?

* **Example JSON configs** for:

  * `cdfci` (single-threaded)
  * `cdfci_omp` (OpenMP)
  * block/multi‑coordinate variants (mCDFCI)
  * excited states (xCDFCI)
  * optimal orbitals (`optorbfci`)
  * analysis tools: `cdfci_tools` (frozen core, symmetry connection)
* **One‑click script**: `run_all.sh` to execute representative demos end‑to‑end.

---

## 2) Prerequisites

* A Unix‑like environment (Linux, macOS, WSL).
* Built CDFCI binaries available in `../bin/` from the `examples/` directory:

  * `cdfci`, `cdfci_omp`, `xcdfci`, `optorbfci`, `cdfci_tools`
* FCIDUMP data files
  * FCIDUMP is a file format which contains the integrals among orbitals of a quantum system. These integrals are usually computed by Hartree-Fock calculation and fully define the Hamiltonian of the quantum system.
  * We provide a few example fcidump files in `../data/*.FCIDUMP`. You may also use your own FCIDUMP data file and specify the path.
* Optional: OpenMP runtime for `cdfci_omp`.

> If you don’t have the binaries yet, see the project’s top‑level **README** for Make instructions.

---

## 3) Quickstart: run everything

From `examples/`:

```bash
./run_all.sh
```

This executes (in order):

```bash
../bin/cdfci            demo_input_cdfci.json
../bin/cdfci_omp        demo_input_cdfci.json
../bin/cdfci            demo_input_mcdfci.json
../bin/cdfci_omp        demo_input_mcdfci.json
../bin/xcdfci           demo_input_xcdfci.json
../bin/optorbfci        demo_input_optorbfci.json
../bin/cdfci_tools      frozen_core.json
../bin/cdfci_tools      symm_conn.json
```

If everything works, you should see banners with version/build info, machine details, Hamiltonian info, parsed JSON, and then solver progress tables.

> **Tip (OpenMP):** Control thread count with `OMP_NUM_THREADS` (e.g., `OMP_NUM_THREADS=24 ./run_all.sh`). Banners will print `OpenMP: enabled with N threads` when active.

---

## 4) Minimal: run a single example

```bash
# single-threaded CDFCI on H2O/STO-3G
../bin/cdfci demo_input_cdfci.json

# OpenMP version of the same case
OMP_NUM_THREADS=16 ../bin/cdfci_omp demo_input_cdfci.json

# multi-coordinate (block) update
../bin/cdfci demo_input_mcdfci.json

# excited states (4 roots)
../bin/xcdfci demo_input_xcdfci.json

# optimal-orbital workflow (OptOrbFCI + RDM)
../bin/optorbfci demo_input_optorbfci.json

# tools: freeze first 2 orbitals
../bin/cdfci_tools frozen_core.json

# tools: symmetry/connection graphs
../bin/cdfci_tools symm_conn.json
```

---

## 5) Understanding the JSON inputs

Most example JSONs follow the structure below (keys may differ per solver):

```json
{
    "hamiltonian": {
        "type": "molecule",
        "molecule": {
            "fcidump_path": "../data/h2o_sto3g_psi4.FCIDUMP",
            "threshold": 0.0,
            "verbose": 2
        }
    },
    "solver":{
        "type": "cdfci",
        "cdfci": {
            "coordinate_pick": "gcd_grad",
            "coordinate_update": "ls",
            "num_coordinates": 1,
            "num_iterations": 150000,
            "report_interval": 10000,
            "stopping_dx_threshold": 1e-10,
            "verbose": 2,
            "z_threshold": 0,
            "z_threshold_search": false
        }
    },
    "max_memory": 0.005
}

```

### Hamiltonian block

* `type`: currently only supports `molecule` or `hubbard_k`.

#### Molecule

* `molecule.fcidump_path`: path to **FCIDUMP** integral file. Default `FCIDUMP`.
* `molecule.threshold`: integral screening threshold (often `0.0` for demos). Default `1e-14`.
* `molecule.verbose`: 0–3 output verbosity. Default `2`.

#### Hubbard\_k

* `hubbard_k.hubbard_t`: nearest-neighbor hopping parameter *t*, controls electron kinetic energy between lattice sites. (Hartree units, default `1`).
* `hubbard_k.hubbard_u`: on-site interaction parameter *U*, electron–electron repulsion strength on the same site. (Hartree units, default `0`).
* `hubbard_k.lattice`: lattice dimensions, e.g. `[4,4]` for a 4×4 square lattice (16 sites, 32 spin-orbitals).
* `hubbard_k.nelec`: total electron number to occupy the lattice.
* `hubbard_k.ms2`: spin projection `2Sz`. `0` means spin balanced (e.g. `N↑=N↓=nelec/2`).


### Solver block (`cdfci`, `xcdfci`, …)

* `coordinate_pick`: coordinate selection rule (e.g., `gcd_grad`, `block_gcd_grad`).
* `coordinate_update`: update strategy (e.g., `ls` line‑search, only works with `gcd_grad`; `eig` block eigen solve).
* `num_coordinates`: 1 for classic CDFCI; >1 for block/multi‑coordinate.
* `num_iterations`, `report_interval`: loop control and logging frequency.
* `stopping_*`: convergence gates (currently supports `stopping_dx_threshold` and `stopping_dx_damping_factor`).
* `ref_det_occ` (optional): the occupied spin-orbitals of the explicit reference determinant; if omitted, **HF** from FCIDUMP is used. reference determinant. The spin-orbitals number from 0.
* `z_threshold` The coefficient compression cutoff epsilon. The larger it is, the smaller the wavefunction is, and the larger the error is. If it is 0, the algorithm will converge to the true ground state. Default `0`.
* `z_threshold_search` If it is true, the package will try to find a `z_threshold` such that the wavefunction can be fitted into the `max_memory` GB memory. The package will increase `z_threshold` by 10 times and restart the calculation if the wavefunction gets too large. Warning: it is experimental and may not always work. Default `false`.

### Others

* `max_memory` The maximum memory (GB) that is available for the wavefunction. The package will determine the size of the wavefunction based on it. Note: the `max_memory` GB memory may not be fully used due to the restriction of the hash table implementation. Suggestion: the maximum available memory. Default 0.5.

> **Tip:** The CDFCI package will try to determine the Hartree Fock (HF) state at
first, if `ref_det_occ` is not provided. If there is the energy of each
spin-orbital in the FCIDUMP (e.g., FCIDUMP generated by Psi4), the HF
state is constructed by filling the orbitals with the lowest energy.
Otherwise (e.g., FCIDUMP generated by PySCF), assume the spin-orbitals are
already sorted by the energy. Then, CDFCI runs from the HF state until
```num_iteration``` iterations are reached.

> We suggest to set ```num_iteration``` to a large number and kill the
process manually after CDFCI converges to the precision you need. We
suggest to set ```max_memory``` as large as it can to get the most
accurate result. The only parameter to tune is ```z_threshold```. We
suggest starting from a small threshold such as 0. If the wavefunction
gets too large, increase ```z_threshold``` and run CDFCI again until
truncated wavefunction can be fitted into the memory. The experemental
feature ```z_threshold_search``` tries to do this job automatically. It
can find a rough threshold most of the time, but it may fail sometimes.
This procedure takes a relatively short time compared with the full
calculation.

### xCDFCI specific

* Top‑level `num_states`: number of states desired, including ground state and excited states. For example, if you want three excited states, put 4 here.
* Block `xcdfci`: same structure as `cdfci`, applied to multi‑state optimization.

### OptOrbFCI specific

* `hamiltonian`: choose `type` and molecule/model input.
* `optimal_orbitals`: iterations/optimizer settings and output folder.
* `rdm`: control 1RDM/2RDM evaluation and method.

### Tools (`cdfci_tools`)

#### Input file for frozen core tools

The input of the executable is a JSON file. Two sample input files are as
follows.
```
{
    "tools": {
        "frozen_orbital": {
            "fcidump_path": "path/to/FCIDUMP",
            "fcidump_frozen_path": "path/to/FCIDUMP_FROZEN",
            "frozen_orbital_list": [0, 1]
        }
    }
}
```

```
{
    "tools": {
        "frozen_orbital": {
            "fcidump_path": "path/to/FCIDUMP",
            "fcidump_frozen_path": "path/to/FCIDUMP_FROZEN",
            "num_frozen_orbital": 2
        }
    }
}
```

- ```fcidump_path``` Path of the FCIDUMP file. Required.
- ```fcidump_frozen_path``` Path of the output FCIDUMP file. Required.
- ```frozen_orbital_list``` A list of spin orbitals to be frozen. Default [].
- ```num_frozen_orbital``` The number of spin orbitals to be frozen. Default 0.

**Note**: Allow only one frozen orbital configuration. If
```frozen_orbital_list``` is non-empty, it is used and
```num_frozen_orbital``` option is ignored. If ```frozen_orbital_list```
is empty and ```num_frozen_orbital``` is adopted.

#### Input file for symmetry connection tools

The input of the executable is a JSON file. Two sample input files are as
follows.
```
{
    "tools": {
        "symmetry_connection": {
            "fcidump_path": "path/to/FCIDUMP",
            "one_body_connection": true,
            "gml_one_body_connection_path": "path/to/one_body_conn.gml",
            "two_body_connection": true,
            "gml_two_body_connection_path": "path/to/two_body_conn.gml",
            "threshold": 1e-10
        }
    }
}
```

- ```fcidump_path``` Path of the FCIDUMP file. Required.
- ```one_body_connection``` If true, calculate the symmetry connection
  among one-body integrals. Optional.
- ```gml_one_body_connection_path``` Path of the one-body connection file
  in sparse graph ```gml``` format. Required if ```one_body_connection```
  is true.
- ```two_body_connection``` If true, calculate the symmetry connection
  among two-body integrals. Optional.
- ```gml_two_body_connection_path``` Path of the two-body connection file
  in sparse graph ```gml``` format. Required if ```two_body_connection```
  is true.
- ```threshold``` Discard integrals below the threshold in the FCIDUMP.
  Default 0.0.


---

## 6) Reading the outputs

Each run prints several standard sections:

1. **Banner & Build Info**
   Shows CDFCI version, compile time, GitHub URL, OpenMP status.

2. **Machine information**
   Key types (`size_t` size), determinant packing (how many `size_t` per det), `__float128` support.

3. **Hamiltonian information**
   Electron count, spin‑orbitals, `Ms2`, HF restricted/unrestricted, warnings about optional FCIDUMP keys.

4. **Parsed JSON**
   Echo of the input blocks actually used.

5. **Reference determinant**
   If not provided, HF state is auto‑selected; prints occupied spin‑orbitals and reference energy.

6. **Progress table** (printed every `report_interval` iterations):

| Column      | Meaning                                           |       |                                                                       |
| ----------- | ------------------------------------------------- | ----- | --------------------------------------------------------------------- |
| `Iteration` | iteration counter                                 |       |                                                                       |
| `Energy`    | current variational energy estimate (Hartree)     |       |                                                                       |
| `dx`        | update size / step norm used in stopping criteria |       |                                                                       |
| `\|x\|_0` | number of non‑zero coefficients in wavefunction (determinant count)   |
|`\|z\|_0` | number of non‑zero coefficients in Hx (determinant count)(candidates)                             |
| `\|H_i\|_0` | number of non-zero entries in active Hamiltonian column(s) touched |
| `Time`      | elapsed wall time since start (seconds)           |       |                                                                       |

7. **Final summary**

   * `Final FCI Energy` (or lowest‑state energy for xCDFCI; multiple lines for multiple states)
   * Cumulative timings by phase: `Energy Estimate`, `Wave Function Update`, `Coordinate Update`, `Coordinate Pick`.

### Typical demo results (sanity checks)

* **H2O / STO‑3G** (single‑state CDFCI): final energy around `-75.71610715…` Hartree.
* **Block/multi‑coordinate** on same case: converges in fewer iterations (e.g., 15k with `num_coordinates=10`), reaching the same final energy.
* **xCDFCI (4 states)**: prints a stack of energies at each report; the first root should converge near the ground‑state value above, with higher roots stabilizing slightly above.
* **OptOrbFCI**: prints OptOrb iterations and RDM summary; for tiny demo sizes you’ll see 2RDM computed with the “square method” and very small runtimes.
* **Tools – Frozen core**: lists frozen + active orbital indices and emits `FCIDUMP_FROZEN`.
* **Tools – Symmetry connection**: writes `orb_one_body_conn.gml` / `orb_two_body_conn.gml` (graph files).

> Note: Exact timings depend on CPU/threads; energies should match within numerical tolerance.

---

## 7) Troubleshooting

* **`fcidump_path` not found** → ensure relative paths are correct from where you run the software (or use absolute paths).
* **OpenMP shows “disabled”** → run the `_omp` binary and/or set `OMP_NUM_THREADS>1`; confirm build used OpenMP.
* **Very slow runs** → reduce `num_iterations`, increase `report_interval`, or try block mode (`num_coordinates>1`).
* **Convergence stalls** → provide a better `ref_det_occ` (reference determinant) or relax `z_threshold`.
* **Memory cap hit** → increase `max_memory` or increase `z_threshold`.
* **Quad‑precision missing** → some platforms lack `__float128`; the code will fall back but results/timings can differ slightly.
* **Graph files not created (tools)** → check write permissions to current dir; verify `*connection*_path` values.

