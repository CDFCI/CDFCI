g++ init_hubbard_model.cpp -o init_hubbard_model -I ../../../src/
./init_hubbard_model
../../../bin/cdfci_omp input.json
rm -f hubbard_4x4_half_filled_init.dat
rm -f init_hubbard_model