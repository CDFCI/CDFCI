# Reference
# http://aip.scitation.org/doi/pdf/10.1063/1.471518
# Table 1, 1.0R_e, first line

from pyscf import gto, scf
from pyscf.tools import fcidump

def generate_fcidump(name, configuration, basis):
    mol = gto.Mole()
    mol.atom = configuration
    mol.basis = basis
    mol.spin = 0
    mol.charge = 0
    mol.build()

    mf = scf.RHF(mol)
    mf.kernel()

    filename = "{}_{}.fcidump".format(name, basis)
    fcidump.from_scf(mf, filename)

    lines_to_insert = [
        " {:.16f}    {}  0  0  0\n".format(mf.mo_energy[i], i + 1)
        for i in range(len(mf.mo_energy))
    ]
    with open(filename, 'r') as f:
        lines = f.readlines()
    lines = lines[:-1] + lines_to_insert + lines[-1:]
    with open(filename, 'w') as f:
        f.writelines(lines)

if __name__ == "__main__":
    for basis in ['cc-pvdz', 'cc-pvtz', 'cc-pvqz']:
        configuration = '''
            O 0 0 -0.004762593
            H 0 0.80184232855 -0.56034446694
            H 0 -0.80184232855 -0.56034446694
            '''
        generate_fcidump("h2o", configuration, basis)
        configuration = '''
            N 0 0 0
            N 0 0 1.12079733
            '''
        generate_fcidump("n2", configuration, basis)


