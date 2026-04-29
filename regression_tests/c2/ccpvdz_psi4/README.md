# C<sub>2</sub>

## Note
This is the same system as reported in the Section 4.1 of the paper
- Z. Wang, Y. Li, J. Lu, [Coordinate Descent Full Configuration
  Interaction](https://pubs.acs.org/doi/10.1021/acs.jctc.9b00138), *J.
  Chem. Theory and Comput.*, *15(6)*, 3558-3569 (2019), DOI:
  10.1021/acs.jctc.9b00138.

You can find the converged output in the ```output``` file.

## Properties
  - Basis: cc-pvdz
  - Number of electrons: 12
  - Number of orbitals: 28
  - Dimension: 1.77e10
  - Hartree-Fock: restricted
  - Hartree-Fock energy: -75.4168819659505800
  - CDFCI converged energy: -75.7319604

## FCIDUMP
The FCIDUMP file is dumped by [Psi4](http://psicode.org/). The input file
for Psi4 v1.1 with FCIDUMP plugin is as follows.
```
import fcidump

molecule {
         C 0 0 0
         C 0 0 1.24253
}

# Reference
# https://aip.scitation.org/doi/pdf/10.1063/1.4928643?class=pdf
# Table 4

set {
    basis cc-pVDZ
    integrals_file 'FCIDUMP'
}

energy('fcidump')
```
