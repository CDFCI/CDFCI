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

#include "driver.h"

int main(int argc, char *argv[])
{
    try
    {
        if (argc < 2)
        {
            throw std::invalid_argument("No input file specified. Please use the path "
                                        "of the input file as the first argument.");
        }
        std::string input_file = argv[1];

        CDFCIProgramDriver cdfci_driver(input_file);
        cdfci_driver.run();
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << std::endl
                  << "ERROR: " << e.what() << std::endl;
        return 1;
    }
}
