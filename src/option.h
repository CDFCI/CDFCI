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

#ifndef CDFCI_OPTION_H
#define CDFCI_OPTION_H

#include <iostream>
#include <exception>
#include <fstream>
#include "config.h"
#include "lib/json/json.hpp"

using std::string_literals::operator""s;
using Option = nlohmann::json;

Option option_from_file(const std::string file_path)
{
    Option result;
    try
    {
        std::ifstream f(file_path);
        f >> result;
    }
    catch (...)
    {
        throw std::invalid_argument("Can not read the input file " + file_path);
    }
    return result;
}

void option_merge_patch(Option &option, Option &&default_option)
{
    if (option.is_null())
        option = Option();
    default_option.merge_patch(option);
    option = default_option;
}

namespace optutil
{

    inline int init_iteration(const Option &opt)
    {
        return opt.at("optimal_orbitals").at("init_iteration").get<int>();
    }

    inline int max_iteration(const Option &opt)
    {
        return opt.at("optimal_orbitals").at("max_iteration").get<int>();
    }

    inline int num_orbitals_compressed(const Option &opt)
    {
        return opt.at("optimal_orbitals").at("num_orbitals").get<int>();
    }

    inline std::string fcidump_path(const Option &opt)
    {
        return opt.at("hamiltonian").at("molecule").at("fcidump_path").get<std::string>();
    }

} // namespace option

Option default_tools_frozen_orbital_option()
{
    Option option = {
        // Required Options
        // "fcidump_path"
        // "fcidump_frozen_path"

        // Optional Options
        {"num_frozen_orbital", 0},
        {"frozen_orbital_list", Option::array()}};
    return option;
}

void validate_tools_frozen_orbital_option(Option &option, std::string const &prefix = ""s)
{
    if (option.is_null())
        throw std::invalid_argument(prefix + " is empty. "
                                             "Please define this option.");

    int num_frozen_orbital = option["num_frozen_orbital"];
    if (num_frozen_orbital < 0)
    {
        throw std::invalid_argument(prefix + "num_frozen_orbital is invalid. "
                                             "Please set a positive number.");
    }
    std::vector<Orbital> frozen_orbital_list = option["frozen_orbital_list"];
    // Allow only one frozen orbital configuration in the end
    // If 'frozen_orbital_list' is non-empty, it is used and
    // 'num_frozen_orbital' option is ignored. If 'frozen_orbital_list' is
    // empty (either from default setting or user setting) and
    // 'num_frozen_orbital' is nonzero, 'num_frozen_orbital' is adopted.
    if (frozen_orbital_list.size() > 0)
    {
        if (option["num_frozen_orbital"] > 0)
            std::cout << "Warning: 'frozen_orbital_list' is used and 'num_frozen_orbital' is ignored."
                      << std::endl
                      << std::endl;
        option.erase("num_frozen_orbital");
    }
    else
    {
        if (option["num_frozen_orbital"] > 0)
            option.erase("frozen_orbital_list");
        else
            option.erase("num_frozen_orbital");
    }
}

Option default_tools_symmetry_connection_option()
{
    Option option = {
        // Required Options
        // "fcidump_path"

        // Optional Options
        {"symmetry_group_list", Option::array()},
        {"one_body_connection", false},
        {"two_body_connection", false},
        {"threshold", 0.0}

        // Conditional Required Options
        // "gml_one_body_connection_path" if "one_body_connection" is true
        // "gml_two_body_connection_path" if "two_body_connection" is true
    };
    return option;
}

void validate_tools_symmetry_connection_option(Option &option, std::string const &prefix = ""s)
{
    if (option.is_null())
        throw std::invalid_argument(prefix + " is empty. "
                                             "Please define this option.");
}

#endif