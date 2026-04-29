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

#include "driver.h"
#include "driver_tools.h"

int main(int argc, char *argv[])
{
    try
    {
        CDFCIProgramDriver::print_header();
        if (argc < 2)
        {
            throw std::invalid_argument("No input file specified. Please use the path "
                                        "of the input file as the first argument.");
        }
        std::string input_file = argv[1];

        print_tools_header();

        Option option = option_from_file(input_file);

        if (option["tools"].contains("frozen_orbital"))
        {
            // Use default value for unset parameters in option.
            option_merge_patch(option["tools"]["frozen_orbital"], default_tools_frozen_orbital_option());
            // Validate the input file.
            validate_tools_frozen_orbital_option(option["tools"]["frozen_orbital"], "tools:frozen_orbital:");
            // Run frozen orbital tool.
            auto fci_frozen = run_tools_frozen_orbital(option["tools"]["frozen_orbital"]);
            fci_frozen.to_file(option["tools"]["frozen_orbital"]["fcidump_frozen_path"]);
        }

        if (option["tools"].contains("symmetry_connection"))
        {
            // Use default value for unset parameters in option.
            option_merge_patch(option["tools"]["symmetry_connection"], default_tools_symmetry_connection_option());
            // Validate the input file.
            validate_tools_symmetry_connection_option(option["tools"]["symmetry_connection"], "tools:symmetry_connection:");
            // Run symmetry connection tool.
            run_tools_symmetry_connection(option["tools"]["symmetry_connection"]);
        }

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << std::endl
                  << "ERROR: " << e.what() << std::endl;
        return 1;
    }
}