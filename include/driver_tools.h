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

#include <iostream>
#include <string>

#include "config.h"
#include "tools.h"
#include "option.h"

void print_tools_header()
{
    std::cout << "[ CDFCI Tools v" << CDFCI_VERSION_MAJOR
              << '.' << CDFCI_VERSION_MINOR << '.' << CDFCI_VERSION_PATCH
              << " | Compiled on " << __DATE__ << " at " << __TIME__ << " ]" << std::endl;

#ifdef _OPENMP
    std::cout << "[ OpenMP: enabled with " << omp_get_max_threads() << " threads ]" << std::endl;
#else
    std::cout << "[ OpenMP: disabled ]" << std::endl;
#endif

    std::cout << "[ GitHub: https://github.com/quan-tum/CDFCI/ ]" << std::endl;
    std::cout << std::endl;
}

Fcidump run_tools_frozen_orbital(Option &option)
{
    std::cout << "Tools Frozen Orbital Option" << std::endl;
    std::cout << "---------------------------" << std::endl;
    std::cout << option.dump(4) << std::endl
              << std::endl;

    Fcidump fci(option["fcidump_path"]);
    // From default option, only one of frozen_orbital_list and
    // num_frozen_orbital exists
    if (option.contains("frozen_orbital_list"))
    {
        std::vector<Orbital> frozen_orbital_list = option["frozen_orbital_list"];
        auto                 fci_frozen          = fcidump_frozen_core(fci, frozen_orbital_list);
        return fci_frozen;
    }
    else if (option.contains("num_frozen_orbital"))
    {
        int  num_frozen_orbital = option["num_frozen_orbital"];
        auto fci_frozen         = fcidump_frozen_core(fci, num_frozen_orbital);
        return fci_frozen;
    }

    return fci;
}

void run_tools_symmetry_connection(Option &option)
{
    std::cout << "Tools Symmetry Connection Option" << std::endl;
    std::cout << "--------------------------------" << std::endl;
    std::cout << option.dump(4) << std::endl
              << std::endl;

    Fcidump fci(option["fcidump_path"]);
    if (option["one_body_connection"])
    {
        std::ofstream f;
        f.open(option["gml_one_body_connection_path"], std::ofstream::out);
        connection_from_one_body_integral(fci,
                                          option["threshold"],
                                          option["symmetry_group_list"],
                                          f);
    }
    if (option["two_body_connection"])
    {
        std::ofstream f;
        f.open(option["gml_two_body_connection_path"], std::ofstream::out);
        connection_from_two_body_integral(fci,
                                          option["threshold"],
                                          option["symmetry_group_list"],
                                          f);
    }
    return;
}