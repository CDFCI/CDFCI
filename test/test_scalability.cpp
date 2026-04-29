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

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cstdlib> // std::stoi
#include "../../src/driver.h"

#include <chrono>

double test_time(std::string path, int num_coordinates)
{
    auto start = std::chrono::high_resolution_clock::now();

    // Read input
    std::string input_file = path + "/input.json";
    CDFCIProgramDriver cdfci_driver(input_file);
    cdfci_driver.opt["solver"]["cdfci"]["num_coordinates"] = num_coordinates;

    // Run CDFCI
    cdfci_driver.run();

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;

    return elapsed.count();
}


int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " [path] [test_count]" << std::endl;
        return 1;
    }

    std::string path = argv[1];
    int test_count = std::stoi(argv[2]);
    std::vector<double> result;

    int coords = 1;
    for (int t = 0; t < test_count; t++)
    {
        double time = test_time(path, coords);
        result.push_back(time);
        coords *= 2;
    }

    std::cout << "Test Results" << std::endl;
    coords = 1;
    for (int t = 0; t < test_count; t++)
    {
        std::cout << "Num of coordinates = " << coords << ", time elapsed = " << std::fixed << std::setprecision(2) << result[t] << " seconds." << std::endl;
        coords *= 2;
    }

    return 0;
}
