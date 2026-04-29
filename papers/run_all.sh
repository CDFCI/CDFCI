#!/bin/bash
#TO-DO: wrap the experiments in functions

export OMP_NUM_THREADS=128 # change this to the number of threads you want

# ../bin/cdfci_omp 2019-cdfci/h2o_ccpvdz.json
# ../bin/cdfci_omp 2019-cdfci/c2_ccpvdz.json
# ../bin/cdfci_omp 2019-cdfci/n2_ccpvdz.json

# ../bin/optorbfci 2020-optorbfci/h2o_ccpvdz_12.json; date
# ../bin/optorbfci 2020-optorbfci/h2o_ccpvdz_24.json; date
# ../bin/optorbfci 2020-optorbfci/h2o_ccpvtz_12.json; date
# ../bin/optorbfci 2020-optorbfci/h2o_ccpvtz_24.json; date
# ../bin/optorbfci 2020-optorbfci/h2o_ccpvqz_12.json; date
# ../bin/optorbfci 2020-optorbfci/h2o_ccpvqz_24.json; date

# ../bin/xcdfci 2023-xcdfci/h2o_ccpvdz.json
# ../bin/xcdfci 2023-xcdfci/n2_ccpvdz_1e-4.json
# ../bin/xcdfci 2023-xcdfci/n2_ccpvdz_1e-5.json

# ../bin/cdfci_omp 2025-mcdfci/c2_ccpvdz.json
# ../bin/cdfci_omp 2025-mcdfci/n2_ccpvdz.json
../bin/cdfci_omp 2025-mcdfci/cr2_ahlrichs.json
