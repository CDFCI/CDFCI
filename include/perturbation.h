//   ______  _______   _______   ______  __
//  /      ||       \ |   ____| /      ||  |
// |  ,----'|  .--.  ||  |__   |  ,----'|  |
// |  |     |  |  |  ||   __|  |  |     |  |
// |  `----.|  '--'  ||  |     |  `----.|  |
//  \______||_______/ |__|      \______||__|
//
// Coordinate Descent Full Configuration Interaction (CDFCI) package in C++17
// https://github.com/quan-tum/CDFCI
//
// Copyright (c) 2019-2026, CDFCI Developers and Contributors
// All rights reserved.
//
// This source code is licensed under the BSD 3-Clause License found in the
// LICENSE file in the root directory of this source tree.

#ifndef CDFCI_PTB_H
#define CDFCI_PTB_H 1

#include <iostream>

#include "wavefunction.h"
#include "hamiltonian.h"

template <typename H, typename W>
class Perturbation
{
public:
    void compute_perturbation_energy(W &vec_xz, H &ham)
    {
        NumericalType ptb_energy = 0.0;
        NumericalType var_energy = vec_xz.get_variational_energy();
        NumericalType x2norm = vecmath::index(vec_xz.get_xx_double(), 0);
        // Compute perturbation energy here
        auto my_functor = [&](const auto& val) {
            auto det = val.first;
            auto x = vecmath::two_norm(vecmath::slice_first_half(val.second));
            auto z = vecmath::two_norm(vecmath::slice_second_half(val.second));

            if (x == 0.0) {
                NumericalType diag_element = ham.get_diagonal(det);
                ptb_energy += z * z / (var_energy - diag_element);
            }
        };
        vec_xz.loop(my_functor);

        std::cout << "Perturbation Energy: " << std::setw(30)
                    << std::right << std::fixed
                    << std::setprecision(16) << ptb_energy/x2norm << std::endl;
        // Total energy = variational energy + perturbation energy
        NumericalType total_energy = var_energy + ptb_energy/x2norm;
        std::cout << "Total Energy: " << std::setw(30)
                    << std::right << std::fixed
                    << std::setprecision(16) << total_energy << std::endl;
    }
};


#endif