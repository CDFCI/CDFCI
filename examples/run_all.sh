#!/bin/bash

../bin/cdfci demo_input_cdfci.json
../bin/cdfci_omp demo_input_cdfci.json
../bin/cdfci demo_input_mcdfci.json
../bin/cdfci_omp demo_input_mcdfci.json
../bin/xcdfci demo_input_xcdfci.json
../bin/optorbfci demo_input_optorbfci.json
../bin/cdfci_tools frozen_core.json
../bin/cdfci_tools symm_conn.json