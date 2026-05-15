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

#include "test.h"
#include "../include/driver.h"

bool test_system(std::string sys_name, std::string path)
{
    // Read input
    std::string input_file = path + "/" + sys_name + "/input.json";
    CDFCIProgramDriver cdfci_driver(input_file);
    // Run CDFCI
    Result result = cdfci_driver.run();
    double energy = result.energy;
    // Read reference
    std::ifstream f(path + "/" + sys_name + "/" + "result.txt");
    double ref_energy, ref_error;
    f >> ref_energy >> ref_error;
    // Compare
    bool passed = (fabs(ref_energy - energy) < ref_error);

    return passed;
}

TEST_CASE("Test example systems")
{
    const std::string test_path = "../regression_tests";

    SUBCASE("c2_ccpvdz_psi4")
    {
        CHECK(test_system("c2/ccpvdz_psi4", test_path));
    }
    SUBCASE("cr2_ahlrichs_psi4")
    {
        CHECK(test_system("cr2/ahlrichs_psi4", test_path));
    }
    SUBCASE("h2o_ccpvdz_psi4")
    {
        CHECK(test_system("h2o/ccpvdz_psi4", test_path));
    }
    SUBCASE("h2o_ccpvdz_pyscf")
    {
        CHECK(test_system("h2o/ccpvdz_pyscf", test_path));
    }
    SUBCASE("h2o_sto3g_psi4")
    {
        CHECK(test_system("h2o/sto3g_psi4", test_path));
    }
    SUBCASE("hubbard")
    {
        CHECK(test_system("hubbard", test_path));
    }
    SUBCASE("n2_ccpvdz_psi4")
    {
        CHECK(test_system("n2/ccpvdz_psi4", test_path));
    }
    SUBCASE("n2_ccpvdz_psi4_eps1e-2")
    {
        CHECK(test_system("n2/ccpvdz_psi4_eps1e-2", test_path));
    }
    SUBCASE("n2_ccpvdz_psi4_triplet")
    {
        CHECK(test_system("n2/ccpvdz_psi4_triplet", test_path));
    }
}