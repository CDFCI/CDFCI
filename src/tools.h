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

#ifndef CDFCI_TOOLS_H
#define CDFCI_TOOLS_H

#include <iostream>
#include <sstream>
#include <algorithm>
#include <exception>
#include "fcidump.h"

Fcidump fcidump_frozen_core(Fcidump const              &fci,
                            std::vector<Orbital> const &frozen_orbital_list)
{
    int nfrozen = frozen_orbital_list.size();
    if (fci.norb < nfrozen)
        throw std::invalid_argument("The number of spin-orbitals cannot be smaller than the number of frozen orbitals.");
    if (fci.nelec < nfrozen)
        throw std::invalid_argument("The number of electrons cannot be smaller than the number of frozen orbitals.");
    for (const auto &idx : frozen_orbital_list)
        if (idx >= fci.norb)
            throw std::invalid_argument("Frozen orbitals is not in the range of spin-orbitals.");

    Fcidump fci_frozen;
    fci_frozen.norb  = fci.norb - nfrozen;
    fci_frozen.nelec = fci.nelec - nfrozen;
    fci_frozen.ms2   = fci.ms2;
    fci_frozen.uhf   = fci.uhf;
    if (!fci.orbital_symmetries.empty())
        fci_frozen.orbital_symmetries.resize(fci_frozen.norb);
    fci_frozen.core_energy = fci.core_energy;

    std::vector<Orbital> orbital_remap(fci.norb);
    Orbital              orbidx = 0;
    for (Orbital i = 0; i < fci.norb; ++i)
        if (std::find(frozen_orbital_list.begin(),
                      frozen_orbital_list.end(), i) != frozen_orbital_list.end())
            orbital_remap[i] = -1;
        else
        {
            if (!fci.orbital_symmetries.empty())
                fci_frozen.orbital_symmetries[orbidx] = fci.orbital_symmetries[i];
            orbital_remap[i] = orbidx;
            orbidx++;
        }

    // Print frozen orbital information
    std::cout << "Frozen Orbital List:\n";
    for (Orbital i = 0; i < fci.norb; ++i)
        if (orbital_remap[i] < 0)
            std::cout << std::setw(4) << i;
    std::cout << std::endl
              << std::endl;
    // Print active orbital old index
    std::cout << "Active Orbital List (original index):\n";
    for (Orbital i = 0; i < fci.norb; ++i)
        if (orbital_remap[i] >= 0)
            std::cout << std::setw(4) << i;
    std::cout << std::endl
              << std::endl;
    // Print active orbital new index
    std::cout << "Active Orbital List (new index):\n";
    for (Orbital i = 0; i < fci.norb; ++i)
        if (orbital_remap[i] >= 0)
            std::cout << std::setw(4) << orbital_remap[i];
    std::cout << std::endl
              << std::endl;

    // Freeze two-body integral part
    // Step 1: Copy active two-body integral part
    for (const auto &ij_entry : fci.two_body_integral)
    {
        auto ij = ij_entry.first;
        auto i  = ij.first;
        auto j  = ij.second;
        auto ii = orbital_remap[i];
        auto jj = orbital_remap[j];
        if (ii < 0 || jj < 0)
            continue;
        for (const auto &ab_entry : ij_entry.second) // key = (a, b), value = integral
        {
            auto ab       = ab_entry.first;
            auto a        = ab.first;
            auto b        = ab.second;
            auto aa       = orbital_remap[a];
            auto bb       = orbital_remap[b];
            auto integral = ab_entry.second;
            if (aa < 0 || bb < 0)
                continue;
            fci_frozen.two_body_integral[{ii, jj}].insert({{aa, bb}, integral});
        }
    }
    // Step 2: Calculate frozen orbital energy contributed from two boday
    // integral part
    for (auto idx_i = 0; idx_i < frozen_orbital_list.size(); ++idx_i)
    {
        for (auto idx_j = idx_i + 1; idx_j < frozen_orbital_list.size(); ++idx_j)
        {
            auto i = frozen_orbital_list[idx_i];
            auto j = frozen_orbital_list[idx_j];
            if (i < j)
                fci_frozen.core_energy += fci.get_two_body_integral(i, j, i, j) - fci.get_two_body_integral(i, j, j, i);
            else if (i > j)
                fci_frozen.core_energy += fci.get_two_body_integral(j, i, j, i) - fci.get_two_body_integral(j, i, i, j);
        }
    }

    // Initialize one_body_integral
    fci_frozen.one_body_integral.resize(fci_frozen.norb);
    for (auto &v : fci_frozen.one_body_integral)
        v.resize(fci_frozen.norb, 0);

    // Freeze one-body integral part
    for (auto i = 0; i < fci.norb; ++i)
    {
        for (auto a = 0; a < fci.norb; ++a)
        {
            auto ii = orbital_remap[i];
            auto aa = orbital_remap[a];
            if (ii >= 0 && aa >= 0)
            {
                fci_frozen.one_body_integral[ii][aa] = fci.one_body_integral[i][a];
                // Two-body integral contributed to one-body integral
                // after frozen orbital
                for (auto idx_k = 0; idx_k < frozen_orbital_list.size(); ++idx_k)
                {
                    auto k = frozen_orbital_list[idx_k];
                    if (i <= k)
                        fci_frozen.one_body_integral[ii][aa] += fci.get_two_body_integral(i, k, a, k) - fci.get_two_body_integral(i, k, k, a);
                    else
                        fci_frozen.one_body_integral[ii][aa] += fci.get_two_body_integral(k, i, k, a) - fci.get_two_body_integral(k, i, a, k);
                }
            }
            else if (i == a)
                fci_frozen.core_energy += fci.one_body_integral[i][a];
        }
    }

    // Freeze orbital energy
    for (const auto &i_entry : fci.orbital_energy)
    {
        auto i      = i_entry.first;
        auto energy = i_entry.second;
        // Only insert unfrozen orbital energy.
        if (orbital_remap[i] >= 0)
            fci_frozen.orbital_energy.insert({orbital_remap[i], energy});
    }

    return fci_frozen;
}

Fcidump fcidump_frozen_core(Fcidump const &fci, int num_frozen_orbital)
{
    std::vector<std::pair<Orbital, NumericalType>> sorted_energy = get_sorted_orbital_energy(fci.orbital_energy, fci.norb);
    std::vector<Orbital>                           frozen_orbital_list;
    for (auto i = 0; i < num_frozen_orbital; ++i)
        frozen_orbital_list.push_back(sorted_energy[i].first);
    return fcidump_frozen_core(fci, frozen_orbital_list);
}

void connection_from_one_body_integral(Fcidump const              &fci,
                                       NumericalType               threshold,
                                       std::vector<Orbital> const &symmetry_group_list,
                                       std::ostream               &stream)
{
    // Only Alpha Spin Integrals are considered

    std::string tab = "    ";

    // Print Graph Header
    stream << "graph [\n";
    stream << tab << "directed 0\n";

    if (!symmetry_group_list.empty())
        if (fci.orbital_symmetries.empty())
            std::cout << "Warning: symmetry_group_list is ignored since the FCIDUMP does not have orbital symmetry information." << std::endl
                      << std::endl;

    // Print Node Information
    for (Orbital i = 0; i < fci.norb; i += 2)
    {
        if (!symmetry_group_list.empty())
        {
            if (!fci.orbital_symmetries.empty())
            {
                if (std::find(symmetry_group_list.begin(),
                              symmetry_group_list.end(),
                              fci.orbital_symmetries[i]) != symmetry_group_list.end())
                    continue;
            }
        }

        stream << tab << "node [\n";
        stream << tab << tab << "id " << i / 2 << "\n";
        if (!fci.orbital_symmetries.empty())
            stream << tab << tab << "gid " << fci.orbital_symmetries[i] << "\n";
        stream << tab << "]\n";
    }

    std::vector<std::vector<bool>> connection(fci.norb / 2);
    for (Orbital i = 0; i < fci.norb; i += 2)
        connection[i / 2].resize(fci.norb / 2, false);

    // Calculate connectivity from one-body
    for (auto i = 0; i < fci.norb; i += 2)
        for (auto a = 0; a < fci.norb; a += 2)
            connection[i / 2][a / 2] = fabs(fci.one_body_integral[i][a]) > threshold;

    // Print Edge Information
    for (Orbital i = 0; i < fci.norb; i += 2)
        for (Orbital a = i + 2; a < fci.norb; a += 2)
            if (connection[i / 2][a / 2])
            {
                if (!symmetry_group_list.empty())
                {
                    if (!fci.orbital_symmetries.empty())
                    {
                        if (std::find(symmetry_group_list.begin(),
                                      symmetry_group_list.end(),
                                      fci.orbital_symmetries[i]) != symmetry_group_list.end())
                            continue;
                    }
                }
                if (!symmetry_group_list.empty())
                {
                    if (!fci.orbital_symmetries.empty())
                    {
                        if (std::find(symmetry_group_list.begin(),
                                      symmetry_group_list.end(),
                                      fci.orbital_symmetries[a]) != symmetry_group_list.end())
                            continue;
                    }
                }
                stream << tab << "edge [\n";
                stream << tab << tab << "source " << i / 2 << "\n";
                stream << tab << tab << "target " << a / 2 << "\n";
                stream << tab << "]\n";
            }

    stream << "]\n";
}

void connection_from_two_body_integral(Fcidump const              &fci,
                                       NumericalType               threshold,
                                       std::vector<Orbital> const &symmetry_group_list,
                                       std::ostream               &stream)
{
    // Only Alpha Spin Integrals are considered
    std::string tab = "    ";

    // Print Graph Header
    stream << "graph [\n";
    stream << tab << "directed 0\n";

    if (!symmetry_group_list.empty())
        if (fci.orbital_symmetries.empty())
            std::cout << "Warning: symmetry_group_list is ignored since the FCIDUMP does not have orbital symmetry information." << std::endl
                      << std::endl;

    int nsym = 0;
    for (auto i : fci.orbital_symmetries)
        if (i >= nsym)
            nsym = i + 1;

    // Print Node Information
    for (Orbital i = 0; i < fci.norb; i += 2)
        for (Orbital j = i; j < fci.norb; j += 2)
        {
            if (!symmetry_group_list.empty())
            {
                if (!fci.orbital_symmetries.empty())
                {
                    if (std::find(symmetry_group_list.begin(),
                                  symmetry_group_list.end(),
                                  fci.orbital_symmetries[i]) != symmetry_group_list.end())
                        continue;
                    if (std::find(symmetry_group_list.begin(),
                                  symmetry_group_list.end(),
                                  fci.orbital_symmetries[j]) != symmetry_group_list.end())
                        continue;
                }
            }

            stream << tab << "node [\n";
            stream << tab << tab << "id " << (i / 2) * (fci.norb / 2) + (j / 2) << "\n";
            stream << tab << tab << "label \"" << (i / 2) << ", " << (j / 2);
            if (!fci.orbital_symmetries.empty())
            {
                stream << " (" << fci.orbital_symmetries[i] << ", " << fci.orbital_symmetries[j] << ")\"\n";
                stream << tab << tab << "gid " << fci.orbital_symmetries[i] * nsym + fci.orbital_symmetries[j] << "\n";
            }
            else
                stream << "\"\n";
            stream << tab << "]\n";
        }

    // Calculate connectivity from two-body
    for (const auto &ij_entry : fci.two_body_integral)
    {
        auto ij = ij_entry.first;
        auto i  = ij.first;
        auto j  = ij.second;
        if (i % 2 == 1 || j % 2 == 1)
            continue; // Skip Beta spin.
        if (i > j)
            continue;
        for (const auto &ab_entry : ij_entry.second)
        {
            auto ab = ab_entry.first;
            auto a  = ab.first;
            auto b  = ab.second;
            if (a % 2 == 1 || b % 2 == 1)
                continue; // Skip Beta spin.
            if (a > b)
                continue;
            if (fabs(ab_entry.second) > threshold)
            {
                if (!symmetry_group_list.empty())
                {
                    if (!fci.orbital_symmetries.empty())
                    {
                        if (std::find(symmetry_group_list.begin(),
                                      symmetry_group_list.end(),
                                      fci.orbital_symmetries[i]) != symmetry_group_list.end())
                            continue;
                        if (std::find(symmetry_group_list.begin(),
                                      symmetry_group_list.end(),
                                      fci.orbital_symmetries[j]) != symmetry_group_list.end())
                            continue;
                        if (std::find(symmetry_group_list.begin(),
                                      symmetry_group_list.end(),
                                      fci.orbital_symmetries[a]) != symmetry_group_list.end())
                            continue;
                        if (std::find(symmetry_group_list.begin(),
                                      symmetry_group_list.end(),
                                      fci.orbital_symmetries[b]) != symmetry_group_list.end())
                            continue;
                    }
                }
                stream << tab << "edge [\n";
                stream << tab << tab << "source " << i / 2 * (fci.norb / 2) + j / 2 << "\n";
                stream << tab << tab << "target " << a / 2 * (fci.norb / 2) + b / 2 << "\n";
                stream << tab << "]\n";
            }
        }
    }

    stream << "]\n";
}

#endif