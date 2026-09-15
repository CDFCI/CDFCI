#!/bin/bash

if [[ -n "${CDFCI_BIN_DIR:-}" ]]; then
    bin_dir="${CDFCI_BIN_DIR}"
elif [[ -x ../build/bin/cdfci ]]; then
    bin_dir=../build/bin
else
    bin_dir=../bin
fi

"${bin_dir}/cdfci" demo_input_cdfci.json
"${bin_dir}/cdfci_omp" demo_input_cdfci.json
"${bin_dir}/cdfci" demo_input_mcdfci.json
"${bin_dir}/cdfci_omp" demo_input_mcdfci.json
"${bin_dir}/xcdfci" demo_input_xcdfci.json
"${bin_dir}/optorbfci" demo_input_optorbfci.json
"${bin_dir}/cdfci_tools" frozen_core.json
"${bin_dir}/cdfci_tools" symm_conn.json
