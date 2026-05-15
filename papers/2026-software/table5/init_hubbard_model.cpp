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
// Copyright (c) 2019-2025, CDFCI Developers and Contributors
// All rights reserved.
//
// This source code is licensed under the BSD 3-Clause License found in the
// LICENSE file in the root directory of this source tree.
//
// Build a balanced initialization for 4x4 half-filled Hubbard model.

#include "../../../include/determinant.h"
#include "../../../include/hamiltonian.h"
#include "../../../include/wavefunction.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{
constexpr int    kDetBlocks = 1; // 4x4 Hubbard has 32 spin orbitals, one 64-bit block is enough.
constexpr double kT         = 1.0;
constexpr double kEnergyEps = 1e-8;

using Orbital = int;

void build_combinations(const std::vector<Orbital> &pool, int choose,
                        int start, std::vector<Orbital> &current,
                        std::vector<std::vector<Orbital>> &output)
{
    if (static_cast<int>(current.size()) == choose)
    {
        output.push_back(current);
        return;
    }

    const int need = choose - static_cast<int>(current.size());
    for (int i = start; i <= static_cast<int>(pool.size()) - need; ++i)
    {
        current.push_back(pool[i]);
        build_combinations(pool, choose, i + 1, current, output);
        current.pop_back();
    }
}

std::vector<std::vector<Orbital>> combinations(const std::vector<Orbital> &pool,
                                               int choose)
{
    std::vector<std::vector<Orbital>> output;
    if (choose < 0 || choose > static_cast<int>(pool.size()))
        return output;
    if (choose == 0)
    {
        output.push_back({});
        return output;
    }

    std::vector<Orbital> current;
    current.reserve(choose);
    build_combinations(pool, choose, 0, current, output);
    return output;
}
} // namespace

int main(int argc, char *argv[])
{
    double U = 4.0;
    if (argc > 1)
        U = std::atof(argv[1]);

    std::vector<int> lattice = {4, 4};
    const int        nelec   = 16; // half-filled: number of electrons equals number of lattice sites.
    const int        ms2     = 0;

    HamiltonianHubbardK<kDetBlocks> hamiltonian;
    hamiltonian.init_HamiltonianHubbardK(kT, U, lattice, nelec, ms2);

    const int norb = hamiltonian.orbital_energy.size();
    const int n_alpha_target = (nelec + ms2) / 2;
    const int n_beta_target  = (nelec - ms2) / 2;

    Determinant<kDetBlocks> base_determinant;
    std::vector<Orbital>    alpha_zero_orbitals;
    std::vector<Orbital>    beta_zero_orbitals;
    int                     n_alpha_base = 0;
    int                     n_beta_base  = 0;

    for (int orb = 0; orb < norb; ++orb)
    {
        const double energy = hamiltonian.orbital_energy[orb];
        if (energy < -kEnergyEps)
        {
            base_determinant.set_orbital(orb);
            if (orb % 2 == 0)
                ++n_alpha_base;
            else
                ++n_beta_base;
        }
        else if (std::abs(energy) <= kEnergyEps)
        {
            if (orb % 2 == 0)
                alpha_zero_orbitals.push_back(orb);
            else
                beta_zero_orbitals.push_back(orb);
        }
    }

    const int n_alpha_add = n_alpha_target - n_alpha_base;
    const int n_beta_add  = n_beta_target - n_beta_base;
    if (n_alpha_add < 0 || n_beta_add < 0)
        throw std::runtime_error("Base determinant overfills electrons for target spin sector.");
    if (n_alpha_add > static_cast<int>(alpha_zero_orbitals.size()) ||
        n_beta_add > static_cast<int>(beta_zero_orbitals.size()))
    {
        throw std::runtime_error("Not enough zero-energy orbitals to complete half-filled initialization.");
    }

    const auto alpha_choices = combinations(alpha_zero_orbitals, n_alpha_add);
    const auto beta_choices  = combinations(beta_zero_orbitals, n_beta_add);

    std::vector<Determinant<kDetBlocks>> determinants;
    determinants.reserve(alpha_choices.size() * beta_choices.size());
    for (const auto &alpha_set : alpha_choices)
    {
        for (const auto &beta_set : beta_choices)
        {
            Determinant<kDetBlocks> det = base_determinant;
            for (auto orb : alpha_set)
                det.set_orbital(orb);
            for (auto orb : beta_set)
                det.set_orbital(orb);
            determinants.push_back(det);
        }
    }

    if (determinants.empty())
        throw std::runtime_error("No determinants generated for 4x4 half-filled initialization.");

    std::cout << "4x4 half-filled Hubbard initialization\n"
              << "  U = " << U << "\n"
              << "  norb = " << norb << ", nelec = " << nelec << "\n"
              << "  base occupied (alpha, beta) = (" << n_alpha_base << ", " << n_beta_base << ")\n"
              << "  add from zero-energy shell (alpha, beta) = (" << n_alpha_add << ", " << n_beta_add << ")\n"
              << "  generated determinants = " << determinants.size() << std::endl;

    using ContainerType = ContainerRobinhood<Determinant<kDetBlocks>, std::array<NumericalType, 2>,
                                             DeterminantHash<kDetBlocks>, DeterminantEqual<kDetBlocks>>;
    WaveFunction<ContainerType, 1> wavefunction;

    const NumericalType coeff = 1.0 / std::sqrt(static_cast<NumericalType>(determinants.size()));
    for (const auto &det_const : determinants)
    {
        auto det = det_const;
        wavefunction.update_x(det, coeff);
    }

    wavefunction.dump_wavefunction("hubbard_4x4_half_filled_init.dat");
    return 0;
}

