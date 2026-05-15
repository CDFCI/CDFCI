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

#ifndef CDFCI_FCIDUMP_H
#define CDFCI_FCIDUMP_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <exception>
#include <type_traits>
#include <algorithm>
#include "lib/json/json.hpp"
#include "determinant.h"

struct pair_hash
{
    template <class T>
    std::size_t operator()(const std::pair<T, T> &p) const
    {
        return std::hash<T>{}((p.first << 10) + p.second);
    }
};

template <typename Orbital, typename ValueType>
std::vector<std::pair<Orbital, ValueType>>
get_sorted_orbital_energy(const std::unordered_map<Orbital, ValueType> &orbital_energy_map, int norb)
{
    std::vector<std::pair<Orbital, ValueType>> sorted_energy;
    if (orbital_energy_map.empty())
    {
        for (auto i = 0; i < norb; ++i)
        {
            sorted_energy.emplace_back(i, i);
        }
        std::cout << "Warning: no orbital energy found. Assume orbitals are sorted."
                  << std::endl
                  << std::endl;
    }
    else
    {
        sorted_energy.assign(orbital_energy_map.begin(), orbital_energy_map.end());
        std::sort(sorted_energy.begin(), sorted_energy.end(),
                  [](const auto &a, const auto &b)
                  {
                      return a.second < b.second;
                  });
    }
    return sorted_energy;
}

class Fcidump
{
public:
    using value_type = double;

    /* FCIDUMP header */

    int              norb;               // Number of spin-orbitals (not orbitals, consistent with UHF).
    int              nelec;              // Number of electrons.
    int              ms2;                // Spin
    std::vector<int> orbital_symmetries; // Optional Orbital symmetries.
    bool             uhf;                // Optional (Psi4/HANDE only). True if Unrestricted Hartree-Fock (UHF) is
                                         // used. False if RHF is used.

    // isym is not used.

    /* FCIDUMP body (integrals) */

    // Store two-body integrals in the format of ((i, j), ((a, b), integral), where i < j.
    std::unordered_map<OrbitalPair, std::unordered_map<OrbitalPair, value_type, pair_hash>, pair_hash> two_body_integral;

    // Store one-body integrals in a norb * norb dense matrix.
    std::vector<std::vector<value_type>> one_body_integral;

    // Store the orbital energy in the format of (i, energy).
    // Note this is optional, which is provided by Psi4/HANDE but not PySCF.
    // Currently, it is only used to find the Hartree-Fock state for initialization.
    std::unordered_map<Orbital, value_type> orbital_energy;

    // The core energy.
    value_type core_energy;

    Fcidump() {}

    Fcidump(int norb_, int nelec_, int ms2_,
            const std::vector<int> &sym  = {},
            bool                    uhf_ = false)
        : norb(norb_), nelec(nelec_), ms2(ms2_),
          orbital_symmetries(sym), uhf(uhf_) {}

    Fcidump(const std::string file_path) { from_file(file_path); }

    void from_file(const std::string file_path)
    {
        try
        {
            std::ifstream f(file_path);
            if (f.fail())
            {
                throw std::invalid_argument("Failed to open the FCIDUMP file " +
                                            file_path);
            }

            // Read FCIDUMP header.
            std::string header, line;
            while (std::getline(f, line) && line.find("&END") == std::string::npos &&
                   line.find("/") == std::string::npos)
            {
                header += " ";
                header += line;
            }
            header += " &END";

            if (f.eof())
            {
                throw std::invalid_argument("FCIDUMP has the wrong header.");
            }

            const std::string int_regex  = R"([ ]*=[ ]*(\d+))";
            const std::string ints_regex = R"([ ]*=[ ]*(([\d]+[ ,\r\n])+))";
            const std::string bool_regex = R"([ ]*=[ .]*(FALSE|TRUE))";

            norb               = read_parameter<int>(header, "NORB", int_regex);
            nelec              = read_parameter<int>(header, "NELEC", int_regex);
            ms2                = read_parameter<int>(header, "MS2", int_regex);
            orbital_symmetries = read_parameters<std::vector<int>>(header, "ORBSYM", ints_regex, std::vector<int>());
            uhf                = read_parameter<bool>(header, "UHF", bool_regex, false);
            if (!uhf)
            {
                norb *= 2; // Convert orbitals to spin-orbitals for RHF.
                if (!orbital_symmetries.empty())
                {
                    orbital_symmetries.resize(norb);
                    for (auto i = norb - 1; i >= 0; i--)
                        orbital_symmetries[i] = orbital_symmetries[i / 2];
                }
            }

            // Initialize one_body_integral
            one_body_integral.resize(norb);
            for (auto &v : one_body_integral)
                v.resize(norb, 0);

            // Read FCIDUMP body (integrals).
            if (uhf)
                read_integral_uhf(f);
            else
                read_integral_rhf(f);
            return;
        }
        catch (const std::exception &e)
        {
            throw std::invalid_argument("Can not read the FCIDUMP file " + file_path +
                                        ". " + e.what());
        }
    }

    static int read_norb(const std::string file_path)
    {
        std::ifstream f(file_path);
        if (f.fail())
        {
            throw std::invalid_argument("FCIDUMP file fails to open: " + file_path);
        }

        // Read FCIDUMP header.
        std::string header, line;
        while (std::getline(f, line) && line.find("&END") == std::string::npos &&
               line.find("/") == std::string::npos)
        {
            header += " ";
            header += line;
        }
        header += " &END";

        if (f.eof())
        {
            throw std::invalid_argument("FCIDUMP has the wrong header.");
        }

        const std::string int_regex  = R"([ ]*=[ ]*(\d+))";
        const std::string bool_regex = R"([ ]*=[ .]*(FALSE|TRUE))";

        auto norb = read_parameter<int>(header, "NORB", int_regex);
        auto uhf  = read_parameter<bool>(header, "UHF", bool_regex, false);
        if (!uhf)
        {
            norb *= 2; // Convert orbitals to spin-orbitals for RHF.
        }
        return norb;
    }

    // Read required parameter
    template <typename T>
    static T read_parameter(const std::string &header, const std::string name,
                            const std::string &regex_string)
    {
        std::regex  r(name + regex_string);
        std::smatch m;
        T           parameter;
        if (std::regex_search(header, m, r))
        {
            std::string parameter_string = m[1];
            // Convert to type T
            if (std::is_same<T, bool>::value)
            {
                // Convert string "TRUE" or "FALSE" (case insensitive) to bool.
                std::transform(parameter_string.begin(), parameter_string.end(),
                               parameter_string.begin(), ::tolower);
                std::istringstream(parameter_string) >> std::boolalpha >> parameter;
            }
            else
            {
                std::istringstream(parameter_string) >> parameter;
            }
        }
        else
        {
            throw std::invalid_argument(name + " is not found.");
        }

        return parameter;
    }

    // Read optional parameter. Default value shoule be provided.
    template <typename T>
    static T read_parameter(const std::string &header, const std::string name,
                            const std::string &regex_string, T default_value)
    {
        T parameter;
        try
        {
            parameter = read_parameter<T>(header, name, regex_string);
        }
        catch (const std::exception &e)
        {
            // Use default value and show a warning. Do not terminate the program.
            parameter = default_value;
            std::cerr << "Warning: FCIDUMP optional parameter " << e.what()
                      << " Use default value " << default_value << "." << std::endl
                      << std::endl;
        }
        return parameter;
    }

    // Read required parameter vector
    // This could be merged with read_parameter function if constexpr in
    // C++17 is supported.
    static std::vector<int> read_parameters(const std::string &header, const std::string name,
                                            const std::string &regex_string)
    {
        std::regex       r(name + regex_string);
        std::smatch      m;
        std::vector<int> parameter;
        if (std::regex_search(header, m, r))
        {
            // Convert to vector
            std::istringstream parameter_ss(m[1]);
            std::string        word;
            while (!parameter_ss.eof())
            {
                std::getline(parameter_ss, word, ',');
                parameter.push_back(stoi(word));
            }
        }
        else
        {
            throw std::invalid_argument(name + " is not found.");
        }

        return parameter;
    }

    // Read optional parameter vector. Default value shoule be provided.
    // This could be merged with read_parameter function if constexpr in
    // C++17 is supported.
    template <typename T>
    static T read_parameters(const std::string &header, const std::string name,
                             const std::string &regex_string, T default_value)
    {
        T parameter;
        try
        {
            parameter = read_parameters(header, name, regex_string);
        }
        catch (const std::exception &e)
        {
            // Use default value and show a warning. Do not terminate the program.
            parameter = default_value;
            std::cerr << "Warning: FCIDUMP optional parameter " << e.what()
                      << " Use default value empty vector." << std::endl
                      << std::endl;
        }
        return parameter;
    }

    void read_integral_uhf(std::ifstream &f)
    {
        Orbital    i, j, a, b;
        value_type integral;
        while (f >> integral >> i >> a >> j >> b)
        {
            // Two-body integrals
            if (i && a && j && b)
            {
                // [symmetry]
                // <ij|g|ab>
                insert_two_body_(i - 1, j - 1, a - 1, b - 1, integral);
                // <ib|g|aj>
                insert_two_body_(i - 1, b - 1, a - 1, j - 1, integral);
                // <aj|g|ib>
                insert_two_body_(a - 1, j - 1, i - 1, b - 1, integral);
                // <ab|g|ij>
                insert_two_body_(a - 1, b - 1, i - 1, j - 1, integral);
            }
            // One-body integrals <i|f|a>
            else if (a)
            {
                one_body_integral[i - 1][a - 1] = integral;
                one_body_integral[a - 1][i - 1] = integral;
            }
            // Orbital energy
            else if (i)
            {
                orbital_energy.insert({i - 1, integral});
            }
            // Core energy
            else
            {
                core_energy = integral;
            }
        }
        return;
    }

    Orbital alpha(const Orbital i) const { return 2 * i - 2; }
    Orbital beta(const Orbital i) const { return 2 * i - 1; }

    void read_integral_rhf(std::ifstream &f)
    {
        Orbital    i, j, a, b;
        value_type integral;
        while (f >> integral >> i >> a >> j >> b)
        {
            // Two-body integrals
            if (i && a && j && b)
            {
                // [symmetry]
                // <ij|g|ab>
                rhf_insert_two_body_(i, j, a, b, integral);
                // <ib|g|aj>
                rhf_insert_two_body_(i, b, a, j, integral);
                // <aj|g|ib>
                rhf_insert_two_body_(a, j, i, b, integral);
                // <ab|g|ij>
                rhf_insert_two_body_(a, b, i, j, integral);
            }
            // One-body integrals <i|f|a>
            else if (a)
            {
                one_body_integral[alpha(i)][alpha(a)] = integral;
                one_body_integral[beta(i)][beta(a)]   = integral;
                one_body_integral[alpha(a)][alpha(i)] = integral;
                one_body_integral[beta(a)][beta(i)]   = integral;
            }
            // Orbital energy
            else if (i)
            {
                orbital_energy.insert({alpha(i), integral});
                orbital_energy.insert({beta(i), integral});
            }
            // Core energy
            else
            {
                core_energy = integral;
            }
        }
        return;
    }

    // Helper function. Insert four combinations of spin-orbitals from orbitals.
    // [symmetry]
    void rhf_insert_two_body_(const Orbital i, const Orbital j, const Orbital a,
                              const Orbital b, const value_type integral)
    {
        // Insert four spin-orbitals combination.
        // alpha-alpha
        insert_two_body_(alpha(i), alpha(j), alpha(a), alpha(b), integral);
        // beta-beta
        insert_two_body_(beta(i), beta(j), beta(a), beta(b), integral);
        // alpha-beta
        insert_two_body_(alpha(i), beta(j), alpha(a), beta(b), integral);
        // beta-alpha
        insert_two_body_(beta(i), alpha(j), beta(a), alpha(b), integral);
        return;
    }

    // Helper function. Keep i <= j.
    // Note: when i = j, <ii|ab> and <ii|ba> are the same and both inserted. [symmetry]
    void insert_two_body_(const Orbital i, const Orbital j, const Orbital a,
                          const Orbital b, const value_type integral)
    {
        if (i <= j)
            two_body_integral[{i, j}].insert({{a, b}, integral});
        // Exchange i, j (and a, b accordingly)
        if (i >= j)
            two_body_integral[{j, i}].insert({{b, a}, integral});
        return;
    }

    void to_file(const std::string file_path)
    {
        try
        {
            std::ifstream existing(file_path);
            if (existing.good())
            {
                std::cerr << "Warning: overwrite existing file: " << file_path << std::endl;
            }
            std::ofstream f;
            f.open(file_path, std::ofstream::out);

            // Write FCIDUMP header.
            f << "&FCI\n";

            // Write required parameters.
            if (uhf)
                f << "NORB=" << norb << ",\n";
            else
                f << "NORB=" << norb / 2 << ",\n";
            f << "NELEC=" << nelec << ",\n";
            if (uhf)
            {
                if (!orbital_symmetries.empty())
                {
                    f << "ORBSYM=";
                    for (auto i : orbital_symmetries)
                        f << i << ",";
                    f << "\n";
                }
            }
            else
            {
                if (!orbital_symmetries.empty())
                {
                    f << "ORBSYM=";
                    // Only write alpha spin orbital symmetrices
                    for (auto i = 0; i < orbital_symmetries.size(); i += 2)
                        f << orbital_symmetries[i] << ",";
                    f << "\n";
                }
            }
            f << "MS2=" << ms2 << ",\n";

            // Write optional parameters.
            if (uhf)
                f << "UHF=.TRUE.,\n";

            // Orbsym is not used and not stored
            // f << "ORBSYM=" << ??? << ",\n";

            f << "&END" << std::endl;

            // Write FCIDUMP body (integrals).
            if (uhf)
                write_integral_uhf(f);
            else
                write_integral_rhf(f);
            return;
        }
        catch (const std::exception &e)
        {
            throw std::invalid_argument("Can not read the FCIDUMP file " + file_path +
                                        ". " + e.what());
        }
    }

    void write_integral_uhf(std::ofstream &f)
    {
        // Write two-body integral values.
        for (const auto &ij_entry : two_body_integral)
        {
            auto ij = ij_entry.first;
            auto i  = ij.first;
            auto j  = ij.second;
            for (const auto &ab_entry : ij_entry.second) // key = (a, b), value = integral
            {
                auto ab            = ab_entry.first;
                auto a             = ab.first;
                auto b             = ab.second;
                auto integral_ijab = ab_entry.second;
                if (i > a || j > b)
                    continue; // Reduce symmetry integral values.
                f << std::scientific << std::setprecision(16) << std::setw(28)
                  << integral_ijab
                  << std::setw(4) << i + 1
                  << std::setw(4) << a + 1
                  << std::setw(4) << j + 1
                  << std::setw(4) << b + 1
                  << "\n";
            }
        }
        f.flush();

        // Write one-body integral values.
        for (auto i = 0; i < norb; ++i)
        {
            for (auto a = 0; a <= i; ++a)
            {
                f << std::scientific << std::setprecision(16) << std::setw(28)
                  << one_body_integral[i][a]
                  << std::setw(4) << i + 1
                  << std::setw(4) << a + 1
                  << std::setw(4) << 0
                  << std::setw(4) << 0
                  << "\n";
            }
        }
        f.flush();

        // Write orbital energy
        for (const auto &i_entry : orbital_energy)
        {
            auto i      = i_entry.first;
            auto energy = i_entry.second;
            f << std::scientific << std::setprecision(16) << std::setw(28)
              << energy
              << std::setw(4) << i + 1
              << std::setw(4) << 0
              << std::setw(4) << 0
              << std::setw(4) << 0
              << "\n";
        }

        // Write core energy
        f << std::scientific << std::setprecision(16) << std::setw(28)
          << core_energy
          << std::setw(4) << 0
          << std::setw(4) << 0
          << std::setw(4) << 0
          << std::setw(4) << 0
          << std::endl;
        return;
    }

    void write_integral_rhf(std::ofstream &f)
    {
        // Write two-body integral values.
        for (const auto &ij_entry : two_body_integral)
        {
            auto ij = ij_entry.first;
            auto i  = ij.first;
            auto j  = ij.second;
            if (i % 2 == 1 || j % 2 == 1)
                continue;                                // Skip Beta spin.
            for (const auto &ab_entry : ij_entry.second) // key = (a, b), value = integral
            {
                auto ab = ab_entry.first;
                auto a  = ab.first;
                auto b  = ab.second;
                if (a % 2 == 1 || b % 2 == 1)
                    continue; // Skip Beta spin.
                auto integral_ijab = ab_entry.second;
                if (i > a || j > b)
                    continue; // Reduce symmetry integral values.
                f << std::scientific << std::setprecision(16) << std::setw(28)
                  << integral_ijab
                  << std::setw(4) << i / 2 + 1
                  << std::setw(4) << a / 2 + 1
                  << std::setw(4) << j / 2 + 1
                  << std::setw(4) << b / 2 + 1
                  << "\n";
            }
        }
        f.flush();

        // Write one-body integral values.
        // Loop over Alpha spin only
        for (auto i = 0; i < norb; i += 2)
        {
            for (auto a = 0; a <= i; a += 2)
            {
                f << std::scientific << std::setprecision(16) << std::setw(28)
                  << one_body_integral[i][a]
                  << std::setw(4) << i / 2 + 1
                  << std::setw(4) << a / 2 + 1
                  << std::setw(4) << 0
                  << std::setw(4) << 0
                  << "\n";
            }
        }
        f.flush();

        // Write orbital energy
        for (const auto &i_entry : orbital_energy)
        {
            auto i      = i_entry.first;
            auto energy = i_entry.second;
            if (i % 2 == 1)
                continue; // Skip Beta spin.
            f << std::scientific << std::setprecision(16) << std::setw(28)
              << energy
              << std::setw(4) << i / 2 + 1
              << std::setw(4) << 0
              << std::setw(4) << 0
              << std::setw(4) << 0
              << "\n";
        }

        // Write core energy
        f << std::scientific << std::setprecision(16) << std::setw(28)
          << core_energy
          << std::setw(4) << 0
          << std::setw(4) << 0
          << std::setw(4) << 0
          << std::setw(4) << 0
          << std::endl;
        return;
    }

    // Return the two_body_integral. Read-only.
    value_type get_two_body_integral(Orbital i, Orbital j, Orbital a, Orbital b) const
    {
        try
        {
            return two_body_integral.at({i, j}).at({a, b});
        }
        catch (const std::out_of_range &e)
        {
            return 0;
        }
    }

    // Return the one_body_integral. Read-only.
    value_type get_one_body_integral(Orbital i, Orbital a) const
    {
        return one_body_integral[i][a];
    }
};

#endif
