# H<sub>2</sub>O

## Note
This is the same system as reported in the Section 4.1 of the paper
- Z. Wang, Y. Li, J. Lu, [Coordinate Descent Full Configuration
  Interaction](https://pubs.acs.org/doi/10.1021/acs.jctc.9b00138), *J.
  Chem. Theory and Comput.*, *15(6)*, 3558-3569 (2019), DOI:
  10.1021/acs.jctc.9b00138.

You can find the converged output in the ```output``` file.

## Properties
  - Basis: cc-pvdz
  - Number of electrons: 10
  - Number of orbitals: 24
  - Dimension: 4.5e8
  - Hartree-Fock: restricted
  - Hartree-Fock energy: -76.0240385602808999
  - CDFCI converged energy: -76.2418601

## FCIDUMP
The FCIDUMP file is dumped by [Psi4](http://psicode.org/). The input file
for Psi4 v1.1 with FCIDUMP plugin is as follows.
```
import fcidump

molecule {
         O 0 0 -0.004762593
         H 0 0.80184232855 -0.56034446694
         H 0 -0.80184232855 -0.56034446694
}

# Reference
# http://aip.scitation.org/doi/pdf/10.1063/1.471518
# Table 1, 1.0R_e, first line

set {
    basis cc-pVDZ
    integrals_file 'FCIDUMP'
}

energy('fcidump')
```