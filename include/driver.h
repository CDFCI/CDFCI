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

#ifndef CDFCI_DRIVER_H
#define CDFCI_DRIVER_H

#include "program.h"
#include "config.h"

#define CDFCI_RUN_CASE(N)           \
    case N:                         \
    {                               \
        CDFCIProgram<N> cdfci(opt); \
        cdfci.print_header();       \
        return cdfci.run();         \
    }

#define XCDFCI_RUN_CASE(N, NSTATES)                         \
    case (N * 100 + NSTATES):                               \
    {                                                       \
        XCDFCIProgram<N, NSTATES> xcdfci(opt, num_states);  \
        xcdfci.print_header();                              \
        return xcdfci.run();                                \
    }

#define OPTORBFCI_RUN_CASE(N)                   \
    case N:                                     \
    {                                           \
        OptOrbFCIProgram<N> optorbfci(opt);     \
        return optorbfci.run();                 \
    }

class CDFCIProgramDriver
{
public:
    // Default option values. They will be updated by the input if provided.
    Option opt = {
        {"hamiltonian",
                 {{"type", "molecule"},
          {"molecule",
           {{"fcidump_path", "FCIDUMP"},
            {"verbose", 2}}},
          {"hubbard_k",
           {{"verbose", 2}}}}},
        {"solver",
                 {{"type", "cdfci"},
          {"cdfci",
           {{"verbose", 2},
            {"max_load_factor", 0.79},
            {"max_wavefunction_size", 1000000}}}}},
        {"num_states", 1},
        {"verbose", 2},
        {"max_memory", 0.5},
        {"max_load_factor", 0.79}};
    int           verbose;
    NumericalType max_memory;
    NumericalType max_load_factor;

    const int bytes_per_size_t = sizeof(CDFCIProgram<1>::det_type);
    const int bits_per_size_t  = bytes_per_size_t * CHAR_BIT;
    const int bytes_per_value  = sizeof(CDFCIProgram<1>::mapped_type);
    int       norb;
    int       size_t_per_det; // The template argument N
    int       bytes_per_entry;

    size_t wavefunction_capacity;
    size_t max_wavefunction_size;

public:
    CDFCIProgramDriver() {}

    CDFCIProgramDriver(std::string input_file)
    {
        read_input(input_file);
    }

    CDFCIProgramDriver(Option option)
    {
        read_input(option);
    }

    ~CDFCIProgramDriver() {}

    void read_input(const std::string input_file)
    {
        Option option = option_from_file(input_file);
        read_input(option);
    }

    void read_input(Option option)
    {
        opt.merge_patch(option);

        verbose         = opt["verbose"];
        max_memory      = opt["max_memory"];
        max_load_factor = opt["max_load_factor"];

        validate_option();

        get_size_t_per_det();
        get_bytes_per_entry();
        get_max_wavefunction_size();
        update_option();

        return;
    }

    void validate_option()
    {
        if (verbose < 0)
        {
            throw std::invalid_argument("verbose (" + std::to_string(verbose) +
                                        ") should be nonnegative.");
        }
        if (max_memory <= 0)
        {
            throw std::invalid_argument("max_memory (" + std::to_string(max_memory) +
                                        ") should be positive.");
        }
        if (max_load_factor <= 0 || max_load_factor >= 1)
        {
            throw std::invalid_argument("max_load_factor (" +
                                        std::to_string(max_load_factor) +
                                        ") should be between 0 and 1.");
        }
        return;
    }

    void get_size_t_per_det()
    {
        norb           = Hamiltonian<1>::read_norb(opt["hamiltonian"]);
        size_t_per_det = 1 + (norb - 1) / bits_per_size_t;
    }

    void get_bytes_per_entry()
    {
        bytes_per_entry = bytes_per_size_t * size_t_per_det + bytes_per_value;
    }

    void get_max_wavefunction_size()
    {
        wavefunction_capacity = 8;
        while (wavefunction_capacity * bytes_per_entry <
               max_memory * 1024ull * 1024 * 1024)
        {
            wavefunction_capacity *= 2;
        }
        wavefunction_capacity /= 2;
        max_wavefunction_size =
            static_cast<size_t>(max_load_factor * wavefunction_capacity);
    }

    void update_option()
    {
        std::string solver = opt["solver"]["type"];
        opt["solver"][solver]["max_load_factor"] = max_load_factor;
        opt["solver"][solver]["max_wavefunction_size"] = max_wavefunction_size;
        return;
    }

    static void print_header()
    {
        std::cout << R"(  ______  _______   _______   ______  __  )" << std::endl;
        std::cout << R"( /      ||       \ |   ____| /      ||  | )" << std::endl;
        std::cout << R"(|  ,----'|  .--.  ||  |__   |  ,----'|  | )" << std::endl;
        std::cout << R"(|  |     |  |  |  ||   __|  |  |     |  | )" << std::endl;
        std::cout << R"(|  `----.|  '--'  ||  |     |  `----.|  | )" << std::endl;
        std::cout << R"( \______||_______/ |__|      \______||__|  )"
                  << CDFCI_VERSION_MAJOR << '.' << CDFCI_VERSION_MINOR << '.'
                  << CDFCI_VERSION_PATCH << std::endl;
        std::cout << std::endl;
        std::cout << "Compiled on " << __DATE__ << " at " << __TIME__ << std::endl;
        std::cout << std::endl;
        std::cout << "Github:" << std::endl;
        std::cout << "    https://github.com/quan-tum/CDFCI/" << std::endl;
        std::cout << "======================================================="
                  << std::endl;
#ifdef _OPENMP
        std::cout << "OpenMP: enabled with " << omp_get_max_threads() << " threads"
                  << std::endl;
#ifdef CDFCI_SOLVER_SERIAL
        std::cout << "OpenMP in CDFCI solver: disabled" << std::endl;
#else
        std::cout << "OpenMP in CDFCI solver: enabled" << std::endl;
#endif
#ifdef CDFCI_RDM_SERIAL
        std::cout << "OpenMP in RDM: disabled" << std::endl;
#else
        std::cout << "OpenMP in RDM: enabled" << std::endl;
#endif
#else
        std::cout << "OpenMP: disabled" << std::endl;
#endif
        std::cout << std::endl;
        return;
    }

    void print_machine_header() const
    {
        std::cout << "Machine information" << std::endl;
        std::cout << "-------------------" << std::endl;
        std::cout << "sizeof(size_t): " << bytes_per_size_t << " bytes or ";
        std::cout << bits_per_size_t << " bits" << std::endl;
        std::cout << "Number of size_t to represent one Slater determinant: "
                  << size_t_per_det << std::endl;
#ifdef CDFCI_USE_LONG_DOUBLE
        std::cout << "__float128: disabled" << std::endl;
        std::cout << "Warning: __float128 may not be supported and long double is "
                     "used instead. Long time iteration (> 1e9) may accumulate "
                     "significant errors. Please use g++/icc or modify "
                     "src/config.h to use the right definition."
                  << std::endl
                  << std::endl;
#else
        std::cout << "__float128: enabled" << std::endl
                  << std::endl;
#endif
    }

    void print_input_option_header() const
    {
        std::cout << "Input option" << std::endl;
        std::cout << "------------" << std::endl;
        std::cout << opt.dump(4) << std::endl;
        std::cout << "Note:" << std::endl;
        std::cout << "In the wavefunction, each determinant needs " << bytes_per_entry
                  << " bytes." << std::endl;
        std::cout << "The size of the wavefunction is 2^"
                  << static_cast<int>(log2(wavefunction_capacity)) << " = "
                  << wavefunction_capacity << ", which can store "
                  << max_wavefunction_size << " determinants at most." << std::endl;
        std::cout << "The wavefunction uses about " << std::fixed << std::setprecision(3)
                  << static_cast<double>(wavefunction_capacity) * bytes_per_entry /
                         1073741824
                  << " GB memory." << std::endl
                  << std::endl;
        return;
    }

    Result run() const
    {
        if (verbose > 1)
        {
            print_header();
            print_machine_header();
        }
        if (verbose > 2)
        {
            print_input_option_header();
        }

        switch (size_t_per_det)
        {
            CDFCI_RUN_CASE(1)
            CDFCI_RUN_CASE(2)
            CDFCI_RUN_CASE(3)
            CDFCI_RUN_CASE(4)
        default:
            throw std::invalid_argument(
                "The CDFCI program only supports 1 - 4 size_t integers for a Slater determinant. \
         Please refer to the documentation, modify, recompile the code and run again.");
            break;
        }
        return Result();
    }
};

class XCDFCIProgramDriver : public CDFCIProgramDriver
{
public:
    int num_states;
    int num_states_preallocated;   // The template argument NSTATES

public:
    XCDFCIProgramDriver() {}

    XCDFCIProgramDriver(std::string input_file)
    {
        read_input(input_file);

        get_size_t_per_det();
        get_num_states_preallocated();
        get_bytes_per_entry();
        get_max_wavefunction_size();
        update_option();
    }

    ~XCDFCIProgramDriver() {}

    void read_input(const std::string input_file)
    {
        CDFCIProgramDriver::read_input(input_file);
        num_states = opt["num_states"];
        // validate
        if (num_states < 1)
        {
            throw std::invalid_argument("num_states (" + std::to_string(num_states) +
                                        ") should be greater than 1.");
        }
        else if (num_states == 1)
        {
            std::cout << "Warning: Only ground-state will be computed. "
                         "In this case, cdfci is more efficient than xcdfci. "
                         "We suggest you to use cdfci instead of xcdfci. "
                         "If this is not intended, specify the number of states "
                         "that you want to compute."
                      << std::endl;
        }
    }

    void get_num_states_preallocated()
    {
        if (num_states <= 2)
            num_states_preallocated = 2;
        else if (num_states <= 4)
            num_states_preallocated = 4;
        else if (num_states <= 8)
            num_states_preallocated = 8;
        else if (num_states <= 16)
            num_states_preallocated = 16;
        else
            throw std::invalid_argument("num_states (" + std::to_string(num_states) +
                                        ") should be smaller than 16.");
    }

    void get_bytes_per_entry()
    {
        bytes_per_entry = bytes_per_size_t * size_t_per_det + bytes_per_value * num_states_preallocated;
    }

    Result run() const
    {
        if (verbose > 1)
        {
            print_header();
            print_machine_header();
        }
        if (verbose > 2)
        {
            print_input_option_header();
        }

        switch (size_t_per_det * 100 + num_states_preallocated)
        {
            XCDFCI_RUN_CASE(1, 2)
            XCDFCI_RUN_CASE(1, 4)
            XCDFCI_RUN_CASE(1, 8)
            XCDFCI_RUN_CASE(1, 16)
            XCDFCI_RUN_CASE(2, 2)
            XCDFCI_RUN_CASE(2, 4)
            XCDFCI_RUN_CASE(2, 8)
            XCDFCI_RUN_CASE(2, 16)
            XCDFCI_RUN_CASE(3, 2)
            XCDFCI_RUN_CASE(3, 4)
            XCDFCI_RUN_CASE(3, 8)
            XCDFCI_RUN_CASE(3, 16)
            XCDFCI_RUN_CASE(4, 2)
            XCDFCI_RUN_CASE(4, 4)
            XCDFCI_RUN_CASE(4, 8)
            XCDFCI_RUN_CASE(4, 16)
        default:
            throw std::invalid_argument(
                "The XCDFCI program only supports 1 - 4 size_t integers for a Slater determinant. \
         Please refer to the documentation, modify, recompile the code and run again.");
            break;
        }
        return Result();
    }
};

class OptOrbFCIProgramDriver : public CDFCIProgramDriver
{
public:
    Option opt_optorb = {
        {"optimal_orbitals",
         {{"num_orbitals", 12},
          {"output_path", "./output"},
          {"init_iteration", 0},
          {"max_iteration", 3},
          {"tol", 1e-4},
          {"optimizer",
           {{"max_rand_iter", 1},
            {"max_iterations", 1000},
            {"report_freq", 10},
            {"tolerance", 1e-7},
            {"initial_stepsize", 1e-4}}}}}};

        OptOrbFCIProgramDriver()
        {
            opt.merge_patch(opt_optorb);
        }

    OptOrbFCIProgramDriver(std::string input_file)
    {
        opt.merge_patch(opt_optorb);
        read_input(input_file);

        get_size_t_per_det();
        get_max_wavefunction_size();
        update_option();
    }

    void get_size_t_per_det()
    {
        norb           = optutil::num_orbitals_compressed(opt) * 2;
        size_t_per_det = 1 + (norb - 1) / bits_per_size_t;

        bytes_per_entry = bytes_per_size_t * size_t_per_det + bytes_per_value;
    }

    NumericalType run() const
    {
        if (verbose > 1)
        {
            print_header();
            print_machine_header();
        }
        if (verbose > 2)
        {
            print_input_option_header();
        }

        switch (size_t_per_det)
        {
            OPTORBFCI_RUN_CASE(1)
            OPTORBFCI_RUN_CASE(2)
            OPTORBFCI_RUN_CASE(3)
            OPTORBFCI_RUN_CASE(4)
        default:
            throw std::invalid_argument(
                "The OptOrbFCI program only supports 1 - 4 size_t integers for a Slater determinant. \
         Please refer to the documentation, modify, recompile the code and run again.");
            break;
        }
        return 0.0;
    }
};

#endif
