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
#ifdef CDFCI_ENABLE_MPI
#include <mpi.h>
#endif

int main(int argc, char *argv[])
{
#ifdef CDFCI_ENABLE_MPI
    int provided = MPI_THREAD_SINGLE;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED)
    {
        std::cerr << "MPI runtime does not support MPI_THREAD_FUNNELED." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
#endif
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
#ifdef CDFCI_ENABLE_MPI
        MPI_Finalize();
#endif
        return 0;
    }
    catch (const std::exception &e)
    {
#ifdef CDFCI_ENABLE_MPI
        if (cdfci_mpi::is_io_root())
        {
#endif
        std::cerr << std::endl
                  << "ERROR: " << e.what() << std::endl;
#ifdef CDFCI_ENABLE_MPI
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
        return 1;
    }
}
