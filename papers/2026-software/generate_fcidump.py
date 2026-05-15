#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate FCIDUMP (frozen core) for H2O with aug-cc-pVTZ using PySCF.

Geometry in bohr (from CC3(full)/aug-cc-pVTZ).
We run RHF, freeze O 1s (1 spatial orbital, 2 electrons), and write FCIDUMP
for the active space Hamiltonian in MO basis.
"""

from pyscf import gto, scf, ao2mo, tools
from pyscf.mcscf import casci


def compute(atom_config, basis_set, path, ncore = 1):
    # ncore = number of frozen spatial orbitals
    # --- Molecule definition (BOHR) ---
    mol = gto.M(
        atom=atom_config,
        unit="Bohr",
        basis=basis_set,
        charge=0,
        spin=0,
        verbose=4,
    )

    # --- RHF ---
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-10
    mf.kernel()
    if not mf.converged:
        raise RuntimeError("RHF did not converge; adjust conv_tol / max_cycle / initial guess.")

    # --- Frozen core setup ---
    norb_total = mf.mo_coeff.shape[1]
    nele_total = mol.nelectron
    nele_active = nele_total - 2 * ncore
    norb_active = norb_total - ncore

    # Build effective (core-folded) 1e Hamiltonian and core energy in MO basis,
    # and get 2e integrals for the active orbitals.
    # casci.CASCI.get_h1eff returns (h1eff, ecore) with core contributions included.
    mc = casci.CASCI(mf, ncas=norb_active, nelecas=nele_active)
    mc.frozen = ncore
    h1eff, ecore = mc.get_h1eff()

    mo_act = mf.mo_coeff[:, ncore:]  # active MO coefficients
    eri_act = ao2mo.kernel(mol, mo_act, aosym="s4")  # (ij|kl) in active MO basis, packed

    # --- Write FCIDUMP ---
    # Ms = 0
    tools.fcidump.from_integrals(path+"/data.fcidump", h1eff, eri_act, norb_active, nele_active,
                                nuc=ecore, ms=0)

    # Ms = 1
    # tools.fcidump.from_integrals(path+"/FCIDUMP_ms1", h1eff, eri_act, norb_active, nele_active,
    #                             nuc=ecore, ms=1)

    print("\nDone.")
    print(f"Frozen spatial orbitals (ncore): {ncore}  (=> frozen electrons: {2*ncore})")
    print(f"Total electrons: {nele_total}  -> Active electrons: {nele_active}")
    print(f"Total orbitals:  {norb_total}  -> Active orbitals:  {norb_active}")


if __name__ == "__main__":
    h2o_atom_config = """
        O  0.00000000  0.00000000 -0.13209669
        H  0.00000000  1.43152878  0.97970006
        H  0.00000000 -1.43152878  0.97970006
        """
    compute(h2o_atom_config, "aug-cc-pvtz", "table2", 1)

    n2_atom_config = """
        N 0.00000000 0.00000000 1.04008632
        N 0.00000000 0.00000000 -1.04008632
        """
    compute(n2_atom_config, "aug-cc-pvdz", "table3", 2)
    c2h2_atom_config = """
        C 0.00000000 0.00000000 1.14048351
        C 0.00000000 0.00000000 -1.14048351
        H 0.00000000 0.00000000 3.14009043
        H 0.00000000 0.00000000 -3.14009043
    """
    compute(c2h2_atom_config, "aug-cc-pvdz", "table4", 2)