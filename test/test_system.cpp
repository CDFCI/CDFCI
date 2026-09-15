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

#define CDFCI_REGRESSION_TEST(name, system)                \
    TEST_CASE(name)                                        \
    {                                                      \
        CHECK(test_system(system, "../regression_tests")); \
    }

CDFCI_REGRESSION_TEST("c2_ccpvdz_psi4", "c2/ccpvdz_psi4")
CDFCI_REGRESSION_TEST("cr2_ahlrichs_psi4", "cr2/ahlrichs_psi4")
CDFCI_REGRESSION_TEST("h2o_ccpvdz_psi4", "h2o/ccpvdz_psi4")
CDFCI_REGRESSION_TEST("h2o_ccpvdz_pyscf", "h2o/ccpvdz_pyscf")
CDFCI_REGRESSION_TEST("h2o_sto3g_psi4", "h2o/sto3g_psi4")
CDFCI_REGRESSION_TEST("hubbard", "hubbard")
CDFCI_REGRESSION_TEST("n2_ccpvdz_psi4", "n2/ccpvdz_psi4")
CDFCI_REGRESSION_TEST("n2_ccpvdz_psi4_eps1e-2", "n2/ccpvdz_psi4_eps1e-2")
CDFCI_REGRESSION_TEST("n2_ccpvdz_psi4_triplet", "n2/ccpvdz_psi4_triplet")

#undef CDFCI_REGRESSION_TEST
