"""PySCF vs CDFCI benchmark on a moderately harder molecular case.

This example is designed to take noticeably longer than the tiny H2O/STO-3G
smoke test and prints a direct runtime comparison.

Run from repository root:
	PYTHONPATH=python python3 python/test.py
"""

from time import perf_counter

import cdfci
from pyscf import fci
from pyscf import gto
from pyscf import mcscf
from pyscf import scf


def build_n2_sto3g_rhf():
	mol = gto.Mole()
	mol.atom = """
	N 0.0 0.0 0.0
	N 0.0 0.0 1.12079733
	"""
	mol.basis = "sto-3g"
	mol.spin = 0
	mol.charge = 0
	mol.build()

	mf = scf.RHF(mol)
	mf.conv_tol = 1e-12

	t0 = perf_counter()
	mf.kernel()
	rhf_time = perf_counter() - t0

	assert mf.converged
	return mf, rhf_time


def run_cdfci_from_pyscf(mf):
	option = {
		"hamiltonian": {
			"type": "molecule",
			"molecule": {
				"threshold": 0.0,
			},
		},
		"solver": {
			"type": "cdfci",
			"cdfci": {
				"num_iterations": 300000,
				"report_interval": 5000,
				"z_threshold": 0.0,
				"z_threshold_search": False,
				"stopping_dx_threshold": 1e-14,
			},
		},
		"max_memory": 1.0,
	}

	t0 = perf_counter()
	res = cdfci.run_from_pyscf_scf(mf, option=option, solver="cdfci")
	elapsed = perf_counter() - t0
	return res, elapsed


def run_pyscf_fci_reference(mf):
	cisolver = fci.FCI(mf)
	t0 = perf_counter()
	e_fci, _ = cisolver.kernel()
	elapsed = perf_counter() - t0
	return float(e_fci), elapsed


def build_n2_ccpvdz_rhf():
	mol = gto.Mole()
	mol.atom = """
	N 0.0 0.0 0.0
	N 0.0 0.0 1.12079733
	"""
	mol.basis = "cc-pvdz"
	mol.spin = 0
	mol.charge = 0
	mol.build()

	mf = scf.RHF(mol)
	mf.conv_tol = 1e-12

	t0 = perf_counter()
	mf.kernel()
	rhf_time = perf_counter() - t0

	assert mf.converged
	return mf, rhf_time


def run_n2_ccpvdz_casci_benchmark():
	# Full-space FCI for N2/cc-pVDZ is typically too expensive.
	# Use a practical active-space benchmark instead.
	ncas = 10
	nelecas = 10

	mf, rhf_time = build_n2_ccpvdz_rhf()

	cas_pyscf = mcscf.CASCI(mf, ncas, nelecas)
	cas_pyscf.fcisolver.conv_tol = 1e-12
	t0 = perf_counter()
	pyscf_ret = cas_pyscf.kernel()
	pyscf_time = perf_counter() - t0
	pyscf_e = float(pyscf_ret[0])

	option = {
		"solver": {
			"type": "cdfci",
			"cdfci": {
				"num_iterations": 500000,
				"report_interval": 10000,
				"z_threshold": 0.0,
				"z_threshold_search": False,
				"stopping_dx_threshold": 1e-14,
			},
		},
		"max_memory": 1.0,
	}

	cas_cdfci = mcscf.CASCI(mf, ncas, nelecas)
	cas_cdfci.fcisolver = cdfci.as_pyscf_fcisolver(option=option, solver="cdfci")
	t0 = perf_counter()
	cdfci_ret = cas_cdfci.kernel()
	cdfci_time = perf_counter() - t0
	cdfci_e = float(cdfci_ret[0])

	print("\n=== PySCF vs CDFCI Benchmark (N2 / cc-pVDZ, CASCI) ===")
	print(f"Active space: ({nelecas}e, {ncas}o)")
	print(f"RHF time (shared pre-step): {rhf_time:.3f} s")
	print(f"PySCF CASCI energy: {pyscf_e:.12f}")
	print(f"CDFCI CASCI energy: {cdfci_e:.12f}")
	print(f"|Delta E|         : {abs(cdfci_e - pyscf_e):.3e} Ha")
	print(f"PySCF CASCI time  : {pyscf_time:.3f} s")
	print(f"CDFCI CASCI time  : {cdfci_time:.3f} s")

	if cdfci_time > 0:
		print(f"Speed ratio (PySCF/CDFCI): {pyscf_time / cdfci_time:.2f}x")


def main():
	mf, rhf_time = build_n2_sto3g_rhf()

	cdfci_res, cdfci_time = run_cdfci_from_pyscf(mf)
	pyscf_e, pyscf_time = run_pyscf_fci_reference(mf)

	assert isinstance(cdfci_res.energy, float)
	assert abs(cdfci_res.energy - pyscf_e) < 1e-2

	print("\n=== PySCF vs CDFCI Benchmark (N2 / STO-3G) ===")
	print(f"RHF time (shared pre-step): {rhf_time:.3f} s")
	print(f"PySCF FCI energy: {pyscf_e:.12f}")
	print(f"CDFCI energy    : {cdfci_res.energy:.12f}")
	print(f"|Delta E|       : {abs(cdfci_res.energy - pyscf_e):.3e} Ha")
	print(f"PySCF FCI time  : {pyscf_time:.3f} s")
	print(f"CDFCI time      : {cdfci_time:.3f} s")

	if cdfci_time > 0:
		print(f"Speed ratio (PySCF/CDFCI): {pyscf_time / cdfci_time:.2f}x")

	run_n2_ccpvdz_casci_benchmark()


if __name__ == "__main__":
	main()