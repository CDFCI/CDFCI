import copy
import json
import os
import shutil
import tempfile

import numpy as np

from _cdfci import run_cdfci_json
from _cdfci import run_cdfci_rhf_integrals_json
from _cdfci import run_xcdfci_json


def _to_int_nelec(nelec):
    if isinstance(nelec, (tuple, list)):
        return int(nelec[0] + nelec[1])
    return int(nelec)


def _to_ms2(nelec):
    if isinstance(nelec, (tuple, list)):
        return int(nelec[0] - nelec[1])
    return 0


class CDFCIPySCFFCISolver:
    """PySCF-compatible FCI solver adapter backed by CDFCI.

    This adapter implements the `kernel` interface expected by PySCF CASCI/CASSCF.
    Internally, it writes a temporary FCIDUMP and calls CDFCI/XCDFCI.
    """

    def __init__(self, option=None, solver="cdfci", keep_fcidump=False):
        self.option = copy.deepcopy(option) if option is not None else {}
        self.solver = solver
        self.keep_fcidump = keep_fcidump

        # Common PySCF FCI-solver attributes.
        self.nroots = 1
        self.spin = None
        self.verbose = 0

        self.last_result = None
        self.last_fcidump_path = None
        self._last_norb = None
        self._cached_rdm1 = None
        self._cached_rdm2 = None

    def _prepare_rdm_option(self, opt, dump_path):
        # Force generation of spin-summed RDMs that PySCF can consume directly.
        opt["rdm"] = {
            "spin_seperate": False,
            "compute_1RDM": True,
            "compute_2RDM": True,
            "method_2RDM": "auto",
            "dump_path": dump_path,
            "verbose": 0,
        }

    @staticmethod
    def _load_cached_rdm(dump_path, norb):
        one_path = os.path.join(dump_path, "OneRDM.out")
        two_path = os.path.join(dump_path, "TwoRDM.out")

        if not os.path.exists(one_path):
            raise RuntimeError("CDFCI RDM output missing OneRDM.out")
        if not os.path.exists(two_path):
            raise RuntimeError("CDFCI RDM output missing TwoRDM.out")

        with open(one_path, "r", encoding="utf-8") as f:
            dm1_vals = np.fromstring(f.read(), sep=" ", dtype=float)
        with open(two_path, "r", encoding="utf-8") as f:
            dm2_vals = np.fromstring(f.read(), sep=" ", dtype=float)

        dm1_expected = norb * norb
        dm2_expected = norb * norb * norb * norb
        if dm1_vals.size != dm1_expected:
            raise RuntimeError(
                "Invalid OneRDM size: got %d, expected %d" % (dm1_vals.size, dm1_expected)
            )
        if dm2_vals.size != dm2_expected:
            raise RuntimeError(
                "Invalid TwoRDM size: got %d, expected %d" % (dm2_vals.size, dm2_expected)
            )

        dm1 = dm1_vals.reshape((norb, norb))
        dm2 = dm2_vals.reshape((norb, norb, norb, norb))
        return dm1, dm2

    def _build_option(self, fcidump_path):
        opt = copy.deepcopy(self.option)

        if "hamiltonian" not in opt:
            opt["hamiltonian"] = {"type": "molecule", "molecule": {}}
        if "molecule" not in opt["hamiltonian"]:
            opt["hamiltonian"]["molecule"] = {}
        opt["hamiltonian"]["molecule"]["fcidump_path"] = fcidump_path

        if "solver" not in opt:
            opt["solver"] = {
                "type": self.solver,
                self.solver: {
                    "num_iterations": 100000,
                    "report_interval": 10000,
                },
            }
        if "type" not in opt["solver"]:
            opt["solver"]["type"] = self.solver

        return opt

    def _run_cdfci_in_memory(self, opt, h1e, eri, norb, nelec, ecore):
        h1e_array = np.asarray(h1e, dtype=float)
        eri_array = np.asarray(eri, dtype=float)
        if eri_array.size != int(norb) ** 4:
            try:
                from pyscf import ao2mo
            except Exception as exc:
                raise ImportError("PySCF is required to normalize ERI tensors for CDFCI") from exc
            eri_array = np.asarray(ao2mo.restore(1, eri_array, int(norb)), dtype=float)
        return run_cdfci_rhf_integrals_json(
            json.dumps(opt),
            h1e_array.reshape(-1).tolist(),
            eri_array.reshape(-1).tolist(),
            int(norb),
            _to_int_nelec(nelec),
            float(ecore),
            _to_ms2(nelec),
        )

    def kernel(self, h1e, eri, norb, nelec, ci0=None, ecore=0.0, **kwargs):
        """Run CDFCI from CASCI/CASSCF active-space integrals."""
        del ci0, kwargs

        self.last_fcidump_path = None
        self._last_norb = int(norb)
        self._cached_rdm1 = None
        self._cached_rdm2 = None

        rdm_dump_dir = tempfile.mkdtemp(prefix="cdfci_rdm_")

        try:
            opt = copy.deepcopy(self.option)
            if "solver" not in opt:
                opt["solver"] = {
                    "type": self.solver,
                    self.solver: {
                        "num_iterations": 100000,
                        "report_interval": 10000,
                    },
                }
            if "type" not in opt["solver"]:
                opt["solver"]["type"] = self.solver

            self._prepare_rdm_option(opt, rdm_dump_dir)

            solver_type = opt.get("solver", {}).get("type", self.solver)
            if solver_type == "cdfci":
                result = self._run_cdfci_in_memory(opt, h1e, eri, norb, nelec, ecore)
            elif solver_type == "xcdfci":
                try:
                    from pyscf.tools import fcidump as pyscf_fcidump
                except Exception as exc:
                    raise ImportError("PySCF is required to use XCDFCI through the adapter") from exc

                tmp = tempfile.NamedTemporaryFile(prefix="cdfci_pyscf_", suffix=".FCIDUMP", delete=False)
                tmp.close()
                fcidump_path = tmp.name
                self.last_fcidump_path = fcidump_path

                try:
                    pyscf_fcidump.from_integrals(
                        fcidump_path,
                        h1e,
                        eri,
                        int(norb),
                        _to_int_nelec(nelec),
                        nuc=float(ecore),
                        ms=_to_ms2(nelec),
                    )
                except TypeError:
                    pyscf_fcidump.from_integrals(
                        fcidump_path,
                        h1e,
                        eri,
                        int(norb),
                        _to_int_nelec(nelec),
                        float(ecore),
                        _to_ms2(nelec),
                    )

                opt = self._build_option(fcidump_path)
                self._prepare_rdm_option(opt, rdm_dump_dir)
                payload = json.dumps(opt)
                result = run_xcdfci_json(payload)
            else:
                raise ValueError("Unsupported solver type for CDFCIPySCFFCISolver: %s" % solver_type)

            self.last_result = result
            self._cached_rdm1, self._cached_rdm2 = self._load_cached_rdm(rdm_dump_dir, self._last_norb)

            # PySCF expects (energy, ci_vector). CDFCI does not expose CI vector yet.
            return float(result.energy), None
        finally:
            if self.last_fcidump_path is not None and (not self.keep_fcidump) and os.path.exists(self.last_fcidump_path):
                os.remove(self.last_fcidump_path)
            if os.path.isdir(rdm_dump_dir):
                shutil.rmtree(rdm_dump_dir)

    def make_rdm1(self, ci, ncas, nelec):
        del ci, nelec
        if self._cached_rdm1 is None:
            raise RuntimeError("RDM is not available. Please call kernel() first.")
        if int(ncas) != self._last_norb:
            raise ValueError("Requested ncas does not match cached CDFCI RDM dimension.")
        return self._cached_rdm1.copy()

    def make_rdm12(self, ci, ncas, nelec):
        del ci, nelec
        if self._cached_rdm1 is None or self._cached_rdm2 is None:
            raise RuntimeError("RDM is not available. Please call kernel() first.")
        if int(ncas) != self._last_norb:
            raise ValueError("Requested ncas does not match cached CDFCI RDM dimension.")
        return self._cached_rdm1.copy(), self._cached_rdm2.copy()


def as_pyscf_fcisolver(option=None, solver="cdfci", keep_fcidump=False):
    """Factory for a PySCF-compatible CDFCI FCI solver object."""
    return CDFCIPySCFFCISolver(
        option=option,
        solver=solver,
        keep_fcidump=keep_fcidump,
    )
