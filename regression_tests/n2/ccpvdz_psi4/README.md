# N<sub>2</sub>

## Note
This is the same system as reported in the Section 4.1 of the paper
- Z. Wang, Y. Li, J. Lu, [Coordinate Descent Full Configuration
  Interaction](https://pubs.acs.org/doi/10.1021/acs.jctc.9b00138), *J.
  Chem. Theory and Comput.*, *15(6)*, 3558-3569 (2019), DOI:
  10.1021/acs.jctc.9b00138.

You can find the converged output in the ```output``` file.

## Properties
  - Basis: cc-pvdz
  - Number of electrons: 14
  - Number of orbitals: 28
  - Dimension: 1.75e11
  - Hartree-Fock: restricted
  - Hartree-Fock energy: -108.9493778784206626
  - CDFCI converged energy: -109.282173

## FCIDUMP
The FCIDUMP file is dumped by [Psi4](http://psicode.org/). The input file
for Psi4 v1.1 with FCIDUMP plugin is as follows.
```
import fcidump

# Equilibruim geometry
# 2.118 a_0 * 0.52917721067 A/a_0 = 1.12079733 A
molecule {
         N 0 0 0
         N 0 0 1.12079733
}

set {
    basis cc-pVDZ
    integrals_file 'FCIDUMP'
}

energy('fcidump')
```
