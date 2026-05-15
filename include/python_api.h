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

#ifndef PYTHON_API_H
#define PYTHON_API_H

#include <fstream>
#include <string>
#include <vector>
#include "driver.h"
#include "driver_tools.h"

class CDFCIDriverFacade {
public:
    explicit CDFCIDriverFacade(const std::string& fcidump_path);

    void set_num_iterations(int n);
    void set_report_interval(int n);
    void set_max_memory(double gb);
    void set_max_load_factor(double x);

    Result run();

private:
    Option opt_;
};

class CDFCIRunnerFacade {
public:
    static Fcidump build_rhf_fcidump_from_integrals(const std::vector<double>& h1e_flat,
                                                    const std::vector<double>& eri_flat,
                                                    int norb,
                                                    int nelec,
                                                    double core_energy,
                                                    int ms2)
    {
        if (norb <= 0)
            throw std::invalid_argument("norb must be positive.");

        const size_t one_body_size = static_cast<size_t>(norb) * static_cast<size_t>(norb);
        const size_t two_body_size = one_body_size * one_body_size;
        if (h1e_flat.size() != one_body_size)
            throw std::invalid_argument("Invalid h1e size for run_cdfci_rhf_integrals_json.");
        if (eri_flat.size() != two_body_size)
            throw std::invalid_argument("Invalid eri size for run_cdfci_rhf_integrals_json.");

        Fcidump fci(2 * norb, nelec, ms2, {}, false);
        fci.core_energy = core_energy;
        fci.one_body_integral.resize(fci.norb);
        for (auto& row : fci.one_body_integral)
            row.resize(fci.norb, 0.0);

        for (int p = 0; p < norb; ++p)
        {
            for (int q = 0; q < norb; ++q)
            {
                const auto value = h1e_flat[static_cast<size_t>(p) * norb + q];
                fci.one_body_integral[2 * p][2 * q] = value;
                fci.one_body_integral[2 * p + 1][2 * q + 1] = value;
            }
        }

        for (int p = 0; p < norb; ++p)
        {
            for (int q = 0; q < norb; ++q)
            {
                for (int r = 0; r < norb; ++r)
                {
                    for (int s = 0; s < norb; ++s)
                    {
                        const auto idx =
                            ((static_cast<size_t>(p) * norb + q) * norb + r) * norb + s;
                        const auto value = eri_flat[idx];
                        if (value == 0.0)
                            continue;

                        // FCIDUMP RHF body stores (p q r s) while the internal tensor is
                        // queried as get_two_body_integral(p, r, q, s).
                        fci.insert_two_body_(2 * p, 2 * r, 2 * q, 2 * s, value);
                        fci.insert_two_body_(2 * p + 1, 2 * r + 1, 2 * q + 1, 2 * s + 1, value);
                        fci.insert_two_body_(2 * p, 2 * r + 1, 2 * q, 2 * s + 1, value);
                        fci.insert_two_body_(2 * p + 1, 2 * r, 2 * q + 1, 2 * s, value);
                    }
                }
            }
        }

        return fci;
    }

    static Result run_cdfci_with_fcidump(Option opt, Fcidump& fci)
    {
        CDFCIProgramDriver driver;
        driver.opt.merge_patch(opt);

        driver.verbose = driver.opt["verbose"];
        driver.max_memory = driver.opt["max_memory"];
        driver.max_load_factor = driver.opt["max_load_factor"];
        driver.validate_option();

        driver.norb = fci.norb;
        driver.size_t_per_det = 1 + (driver.norb - 1) / driver.bits_per_size_t;
        driver.get_bytes_per_entry();
        driver.get_max_wavefunction_size();
        driver.update_option();

        if (driver.verbose > 1)
        {
            driver.print_header();
            driver.print_machine_header();
        }
        if (driver.verbose > 2)
        {
            driver.print_input_option_header();
        }

        switch (driver.size_t_per_det)
        {
        case 1:
        {
            CDFCIProgram<1> cdfci(driver.opt, fci);
            cdfci.print_header();
            return cdfci.run();
        }
        case 2:
        {
            CDFCIProgram<2> cdfci(driver.opt, fci);
            cdfci.print_header();
            return cdfci.run();
        }
        case 3:
        {
            CDFCIProgram<3> cdfci(driver.opt, fci);
            cdfci.print_header();
            return cdfci.run();
        }
        case 4:
        {
            CDFCIProgram<4> cdfci(driver.opt, fci);
            cdfci.print_header();
            return cdfci.run();
        }
        default:
            throw std::invalid_argument(
                "The CDFCI program only supports 1 - 4 size_t integers for a Slater determinant.");
        }
    }

    static Result run_cdfci_json(const std::string& option_json)
    {
        Option opt = Option::parse(option_json);
        CDFCIProgramDriver driver(opt);
        return driver.run();
    }

    static Result run_xcdfci_json(const std::string& option_json)
    {
        Option opt = Option::parse(option_json);
        XCDFCIProgramDriver driver;
        driver.CDFCIProgramDriver::read_input(opt);
        driver.num_states = opt["num_states"];
        driver.get_size_t_per_det();
        driver.get_num_states_preallocated();
        driver.get_bytes_per_entry();
        driver.get_max_wavefunction_size();
        driver.update_option();
        return driver.run();
    }

    static Result run_cdfci_rhf_integrals_json(const std::string& option_json,
                                               const std::vector<double>& h1e_flat,
                                               const std::vector<double>& eri_flat,
                                               int norb,
                                               int nelec,
                                               double core_energy,
                                               int ms2)
    {
        Option opt = Option::parse(option_json);
        auto fci = build_rhf_fcidump_from_integrals(h1e_flat, eri_flat, norb, nelec, core_energy, ms2);
        return run_cdfci_with_fcidump(opt, fci);
    }

    static double run_optorbfci_json(const std::string& option_json)
    {
        Option opt = Option::parse(option_json);

        OptOrbFCIProgramDriver driver;
        driver.opt.merge_patch(driver.opt_optorb);
        driver.read_input(opt);
        driver.get_size_t_per_det();
        driver.get_max_wavefunction_size();
        driver.update_option();

        return driver.run();
    }

    static void tools_frozen_orbital_json(const std::string& option_json)
    {
        Option option = Option::parse(option_json);

        option_merge_patch(option, default_tools_frozen_orbital_option());
        validate_tools_frozen_orbital_option(option, "tools:frozen_orbital:");

        auto fci_frozen = run_tools_frozen_orbital(option);
        fci_frozen.to_file(option["fcidump_frozen_path"]);
    }

    static void tools_symmetry_connection_json(const std::string& option_json)
    {
        Option option = Option::parse(option_json);

        option_merge_patch(option, default_tools_symmetry_connection_option());
        validate_tools_symmetry_connection_option(option, "tools:symmetry_connection:");

        run_tools_symmetry_connection(option);
    }
};

CDFCIDriverFacade::CDFCIDriverFacade(const std::string& fcidump_path) {
    opt_ = {
        {"hamiltonian",
         {{"type", "molecule"},
          {"molecule",
           {{"fcidump_path", fcidump_path}}}}},
        {"solver",
         {{"type", "cdfci"},
          {"cdfci",
           {{"num_iterations", 100000},
            {"report_interval", 10000}}}}},
        {"num_states", 1},
        {"max_memory", 0.5},
        {"max_load_factor", 0.79}
    };
}

void CDFCIDriverFacade::set_num_iterations(int n) {
    opt_["solver"]["cdfci"]["num_iterations"] = n;
}

void CDFCIDriverFacade::set_report_interval(int n) {
    opt_["solver"]["cdfci"]["report_interval"] = n;
}

void CDFCIDriverFacade::set_max_memory(double gb) {
    opt_["max_memory"] = gb;
}

void CDFCIDriverFacade::set_max_load_factor(double x) {
    opt_["max_load_factor"] = x;
}

Result CDFCIDriverFacade::run() {
    CDFCIProgramDriver cdfci_driver(opt_);
    auto res = cdfci_driver.run();

    Result out;
    out.energy = res.energy;
    return out;
}

#endif