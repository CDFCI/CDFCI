import copy
import json
import os
import tempfile

import numpy as np

from _cdfci import CDFCI
from _cdfci import CDFCIResult
from _cdfci import run_cdfci_json
from _cdfci import run_cdfci_rhf_integrals_json
from _cdfci import run_optorbfci_json
from _cdfci import run_xcdfci_json
from _cdfci import tools_frozen_orbital_json
from _cdfci import tools_symmetry_connection_json


def _load_option(option):
    if isinstance(option, dict):
        return copy.deepcopy(option)

    if isinstance(option, str):
        if os.path.isfile(option):
            with open(option, "r", encoding="utf-8") as f:
                return json.load(f)
        return json.loads(option)

    raise TypeError("option must be dict, JSON string, or path to a JSON file")


def _dump_option(option_dict):
    return json.dumps(option_dict)


def run_cdfci(option):
    """Run CDFCI with a dict / JSON string / JSON file path."""
    opt = _load_option(option)
    return run_cdfci_json(_dump_option(opt))


def run_xcdfci(option):
    """Run XCDFCI with a dict / JSON string / JSON file path."""
    opt = _load_option(option)
    return run_xcdfci_json(_dump_option(opt))


def run_cdfci_rhf_integrals(option, h1e, eri, norb, nelec, core_energy=0.0, ms2=0):
    """Run CDFCI directly from RHF spatial-orbital integrals without FCIDUMP I/O."""
    opt = _load_option(option)

    h1e_array = np.asarray(h1e, dtype=float)
    if h1e_array.shape != (norb, norb):
        raise ValueError("h1e must have shape (norb, norb)")

    eri_array = np.asarray(eri, dtype=float)
    if eri_array.shape != (norb, norb, norb, norb):
        raise ValueError("eri must have shape (norb, norb, norb, norb)")

    return run_cdfci_rhf_integrals_json(
        _dump_option(opt),
        h1e_array.reshape(-1).tolist(),
        eri_array.reshape(-1).tolist(),
        int(norb),
        int(nelec),
        float(core_energy),
        int(ms2),
    )


def run_optorbfci(option):
    """Run OptOrbFCI with a dict / JSON string / JSON file path."""
    opt = _load_option(option)
    return run_optorbfci_json(_dump_option(opt))


def tools_frozen_orbital(option):
    """Run frozen orbital tool with a dict / JSON string / JSON file path."""
    opt = _load_option(option)
    tools_frozen_orbital_json(_dump_option(opt))


def tools_symmetry_connection(option):
    """Run symmetry connection tool with a dict / JSON string / JSON file path."""
    opt = _load_option(option)
    tools_symmetry_connection_json(_dump_option(opt))


def run_solver(option):
    """
    Dispatch solver by option["solver"]["type"].
    Returns CDFCIResult for cdfci/xcdfci.
    """
    opt = _load_option(option)

    solver_type = opt.get("solver", {}).get("type", "cdfci")
    if solver_type == "cdfci":
        return run_cdfci(opt)
    if solver_type == "xcdfci":
        return run_xcdfci(opt)
    raise ValueError("Unsupported solver.type for run_solver: %s" % solver_type)


def run_from_pyscf_scf(
    mf,
    option=None,
    solver="cdfci",
    fcidump_path=None,
    keep_fcidump=False,
):
    """
    Generate FCIDUMP from a PySCF SCF object and run CDFCI family solvers.

    Parameters
    ----------
    mf:
    A converged PySCF SCF object.
    option: dict | None
    CDFCI option dict. If omitted, a minimal default is used.
    solver: str
    One of "cdfci" or "xcdfci".
    fcidump_path: str | None
    Optional path to write FCIDUMP.
    keep_fcidump: bool
        Keep generated temporary FCIDUMP when fcidump_path is None.
    """
    try:
        from pyscf import ao2mo
        from pyscf.tools import fcidump as pyscf_fcidump
    except Exception as exc:
        raise ImportError("PySCF is required for run_from_pyscf_scf") from exc

    user_opt = (
        copy.deepcopy(option)
        if option is not None
        else {
            "hamiltonian": {
                "type": "molecule",
                "molecule": {},
            },
            "solver": {
                "type": solver,
                solver: {
                    "num_iterations": 100000,
                    "report_interval": 10000,
                },
            },
            "max_memory": 0.5,
        }
    )

    if "hamiltonian" not in user_opt:
        user_opt["hamiltonian"] = {"type": "molecule", "molecule": {}}
    if "molecule" not in user_opt["hamiltonian"]:
        user_opt["hamiltonian"]["molecule"] = {}
    if "solver" not in user_opt:
        user_opt["solver"] = {"type": solver, solver: {}}
    if "type" not in user_opt["solver"]:
        user_opt["solver"]["type"] = solver

    solver_type = user_opt.get("solver", {}).get("type", solver)

    if solver_type == "cdfci" and getattr(mf, "mo_coeff", None) is not None:
        mo_coeff = np.asarray(mf.mo_coeff)
        if mo_coeff.ndim == 2:
            h1e = mo_coeff.T @ mf.get_hcore() @ mo_coeff
            eri = ao2mo.restore(1, ao2mo.full(mf.mol, mo_coeff), mo_coeff.shape[1])
            return run_cdfci_rhf_integrals(
                user_opt,
                h1e,
                eri,
                mo_coeff.shape[1],
                mf.mol.nelectron,
                core_energy=mf.mol.energy_nuc(),
                ms2=mf.mol.spin,
            )

    if fcidump_path is None:
        tmp = tempfile.NamedTemporaryFile(prefix="cdfci_", suffix=".FCIDUMP", delete=False)
        tmp.close()
        fcidump_path = tmp.name
        remove_after = not keep_fcidump
    else:
        remove_after = False

    pyscf_fcidump.from_scf(mf, fcidump_path)
    user_opt["hamiltonian"]["molecule"]["fcidump_path"] = fcidump_path

    try:
        if solver_type == "cdfci":
            return run_cdfci(user_opt)
        if solver_type == "xcdfci":
            return run_xcdfci(user_opt)
        raise ValueError("Unsupported solver type in run_from_pyscf_scf: %s" % solver_type)
    finally:
        if remove_after and os.path.exists(fcidump_path):
            os.remove(fcidump_path)


__all__ = [
    "CDFCI",
    "CDFCIResult",
    "run_cdfci",
    "run_cdfci_rhf_integrals",
    "run_xcdfci",
    "run_optorbfci",
    "tools_frozen_orbital",
    "tools_symmetry_connection",
    "run_solver",
    "run_from_pyscf_scf",
]

# Optional PySCF adapter import. Keep it lazy/fail-safe when PySCF is absent.
try:
    from .pyscf_solver import CDFCIPySCFFCISolver, as_pyscf_fcisolver

    __all__.extend([
        "CDFCIPySCFFCISolver",
        "as_pyscf_fcisolver",
    ])
except Exception:
    pass