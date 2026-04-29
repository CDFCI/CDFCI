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

#ifndef CDFCI_SOLVER_H
#define CDFCI_SOLVER_H 1

#include "coordinate_strategy.h"
#include "hamiltonian.h"
#include "wavefunction.h"
#include <chrono>
#include <iomanip>
#include <memory>

/**
 * @brief Coordinate Descent Full Configuration Interaction (CDFCI) Solver Framework.
 *
 * This class provides a general-purpose solver framework for performing full
 * configuration interaction (FCI) via adaptive coordinate descent. The core logic
 * is implemented in the solve() functions, which follow a standard iterative
 * coordinate descent procedure:
 *
 * 1. Pick coordinates (determinants) to update.
 * 2. Compute update directions and magnitudes.
 * 3. Apply updates to the wave function.
 * 4. Check stopping criteria and output results.
 *
 * This class is meant to be subclassed. The member functions that are
 * declared as `virtual` must be implemented in the derived class
 * to define problem-specific behavior, for example:
 *
 *   - `initialize(...)`: Initialize wavefunction before the main loop.
 *   - `check_convergence(...)`: Define custom stopping criteria.
 *   - `print_solver_header()`: (Optional) Customize solver output headers.
 *   - ...
 *
 * The high-level interface allows `solve(h)` or `solve(h, vec_xz)` to run
 * a complete CDFCI solve, handling internal initialization and threshold
 * adjustment logic.
 *
 * Template Parameters:
 *   - H: Hamiltonian type (must implement required interface).
 *   - W: WaveFunction type (must implement coordinate access and updates).
 */
template <typename H, typename W>
class Solver
{
public:
    using det_type = typename W::key_type;
    using wff_type = typename W::wff_type;
    using data_type = typename W::data_type;

    /* Options */
    Option opt = {{"num_iterations", 10},
                  {"report_interval", 1},
                  {"num_coordinates", 1},
                  {"ref_det_occ", Option::array()},
                  {"z_threshold", 0.0},
                  {"z_threshold_search", false},
                  {"max_wavefunction_size", 1000000},
                  {"max_load_factor", 0.79},
                  {"stopping_dx_threshold", 1e-10},
                  {"stopping_dx_damping_factor", 0.0},
                  {"load_configuration", false},
                  {"load_configuration_location", "./wavefunction.in"},
                  {"recalculate_z_on_load", false},
                  {"save_configuration", false},
                  {"save_configuration_location", "./wavefunction.out"},
                  {"verbose", 2}};

    size_t num_iter;                         ///< Total number of coordinate descent iterations to perform.
    size_t report_interval;                  ///< Interval for reporting progress and checking convergence.
    std::string coord_pick_str;              ///< Strategy name for coordinate (determinant) selection.
    std::string coord_update_str;            ///< Strategy name for coordinate update rule.
    size_t num_coordinates;                  ///< Number of coordinates to descend per iteration.

    std::vector<Orbital> ref_det_occ;        ///< Orbital occupation of reference determinant (e.g. Hartree–Fock).

    NumericalType z_threshold;               ///< Threshold for pruning small z components during update.
    bool z_threshold_search;                 ///< Whether to automatically adjust z_threshold when wavefunction overflows.

    size_t max_wavefunction_size;            ///< Maximum allowed size of wavefunction storage
    NumericalType max_load_factor;           ///< Load factor for internal wavefunction storage (< 1).

    NumericalType stopping_dx_damping_factor;  ///< Damping factor to control convergence rate of dx accumulation.
    NumericalType stopping_dx_threshold;       ///< Threshold for convergence of accumulated coordinate updates.

    int verbose = 2;

    bool shifted = false; ///< Whether the Hamiltonian has been shifted.
    NumericalType shift_value = 0.0; ///< Value of the energy shift applied to the Hamiltonian.

    static std::unique_ptr<Solver<H, W>> init(Option &option);

    Solver() {}
    virtual ~Solver() {}

    // helper function
    void init_Solver()
    {
        num_iter                   = opt["num_iterations"];
        report_interval            = opt["report_interval"];
        ref_det_occ                = opt.at("ref_det_occ").get<std::vector<int>>();
        coord_pick_str             = opt["coordinate_pick"];
        coord_update_str           = opt["coordinate_update"];
        num_coordinates            = opt["num_coordinates"];
        z_threshold                = opt["z_threshold"];
        z_threshold_search         = opt["z_threshold_search"];
        max_wavefunction_size      = opt["max_wavefunction_size"];
        max_load_factor            = opt["max_load_factor"];
        stopping_dx_damping_factor = opt["stopping_dx_damping_factor"];
        stopping_dx_threshold      = opt["stopping_dx_threshold"];
        verbose                    = opt["verbose"];
    }

    void validate_option()
    {
        if (report_interval < 1)
        {
            throw std::invalid_argument("solver.report_interval (" +
                                        std::to_string(report_interval) +
                                        ") should be at least 1.");
        }
        if (num_iter < report_interval)
        {
            throw std::invalid_argument(
                "solver.num_iterations (" + std::to_string(num_iter) +
                ") should be no smaller than solver.report_interval.");
        }
        if (num_coordinates < 1)
        {
            throw std::invalid_argument("solver.num_coordinates (" +
                                        std::to_string(num_coordinates) +
                                        ") should be at least 1.");
        }
        if (num_coordinates > 1 && coord_pick_str == "gcd_grad")
        {
            throw std::invalid_argument("solver.coord_pick_str (" +
                                        coord_pick_str +
                                        ") cannot be gcd_grad when " +
                                        "solver.num_coordinates > 1.");
        }
        if (z_threshold < 0)
        {
            throw std::invalid_argument("solver.z_threshold (" +
                                        std::to_string(z_threshold) +
                                        ") should be nonnegative.");
        }
        if (stopping_dx_threshold < 0)
        {
            throw std::invalid_argument("solver.stopping_dx_threshold (" +
                                        std::to_string(stopping_dx_threshold) +
                                        ") should be nonnegative.");
        }
        if (stopping_dx_damping_factor < 0 || stopping_dx_damping_factor > 1)
        {
            throw std::invalid_argument("solver.stopping_dx_damping_factor (" +
                                        std::to_string(stopping_dx_damping_factor) +
                                        ") should be between 0 and 1.");
        }
        if (max_load_factor <= 0 || max_load_factor >= 1)
        {
            throw std::invalid_argument("solver.max_load_factor (" +
                                        std::to_string(max_load_factor) +
                                        ") should be between 0 and 1.");
        }
        // Warnings
        if (z_threshold_search)
        {
            std::cout << "Warning: z cut-off threshold auto search (z_threshold_search) "
                         "is enabled. ";
            std::cout << "This feature is experimental and may not always work."
                      << std::endl
                      << std::endl;
            const double start_z_threshold = 1e-10;
            if (z_threshold < start_z_threshold)
            {
                z_threshold = start_z_threshold;
                std::cout << "Warning: The initial z_threshold may be too small. Set "
                             "z_threshold = "
                          << std::scientific << std::setprecision(2) << start_z_threshold
                          << std::endl
                          << std::endl;
            }
        }
        if (z_threshold_search && (report_interval < 1000))
        {
            throw std::invalid_argument(
                "Please set report_interval >= 1000 when z_threshold_search = true.");
        }
    }

    virtual void print_solver_header() const
    {
        std::cout << "Solver information" << std::endl;
        std::cout << "-------------------" << std::endl;
    }

    void output_header() const
    {
        if (verbose > 0)
        {
            std::cout << std::setw(13) << std::left << "Iteration";
            std::cout << std::setw(18) << std::right << "Energy";
            std::cout << std::setw(18) << "dx";
            std::cout << std::setw(15) << "|x|_0";
            std::cout << std::setw(15) << "|z|_0";
            std::cout << std::setw(9) << "|H_i|_0";
            std::cout << std::setw(10) << "Time";
            std::cout << std::endl;
        }
    }

    virtual void output_line(const size_t iteration,
                             const data_type energy,
                             const NumericalType dx,
                             const size_t x_size,
                             const size_t z_size,
                             const size_t H_i_size,
                             const double time
                             ) const = 0;

    virtual void initialize(H &h, W &vec_xz, wff_type &sub_xz) = 0;

    virtual NumericalType get_default_dx_accumulated() const = 0;

    virtual bool check_convergence(const W &vec_xz,
                                   NumericalType &dx_accumulated,
                                   const wff_type &det_picked) = 0;

    virtual void output_final(const W &vec_xz,
                              const size_t &iterations) = 0;

    virtual NumericalType get_ground_state_energy() const = 0;

    int solve(H                      &h,
              CoordinatePick<W>      &coord_pick,
              CoordinateUpdate<H, W> &coord_update,
              W                      &vec_xz,
              NumericalType           z_threshold_)
    {
        size_t    iterations;
        data_type energy;

        typename H::Column column;
        size_t             last_xz_size = 0;

        Timer timer;
        auto  time_start = std::chrono::high_resolution_clock::now();

        // Initialize sub_xz
        wff_type sub_xz;
        if (vec_xz.size() == 0)
            initialize(h, vec_xz, sub_xz);
        else {
            // Copy everything in vec_xz to sub_xz
            auto move_wff_functor = [&](const auto& val) {
                sub_xz.push_back(val.first, val.second);
            };
            vec_xz.loop(move_wff_functor);
        }

        // Checking the inital energy is negative, see Appendix A of paper MCDFCI

        NumericalType curr_var_energy = vecmath::index(vec_xz.get_variational_energy(), 0);
        if (curr_var_energy > 0) {
            std::cout << "Warning: The initial variational energy is positive: "
                      << std::scientific << std::setprecision(6)
                      << curr_var_energy << std::endl;
            std::cout << "Doing a shift to the Hamiltonian (invisible to user)."
                      << std::endl;
            shifted = true;
            shift_value = curr_var_energy + 1.0;
            h.apply_shift(shift_value);
        }

        // Initialize stopping criterion
        NumericalType dx_accumulated = get_default_dx_accumulated();
        bool early_stop = false;

        output_header();

        for (auto i = 0; i < num_iter / report_interval; ++i)
        {
            for (auto j = 0; j < report_interval; ++j)
            {
                iterations = i * report_interval + j + 1;

                // Coordinate Pick
                // det_picked : wff<{det, {x, z}}>
                wff_type det_picked;
                timer.measure("Coordinate Pick", [&]()
                              { det_picked = coord_pick(vec_xz, sub_xz, num_coordinates); });

                // Coordinate Update
                ScaleFactorType scale_factor;
                timer.measure("Coordinate Update", [&]()
                              { scale_factor = coord_update(h, vec_xz, det_picked); });

                // Update Wavefunction
                timer.measure("Wave Function Update", [&]()
                              {
                    vec_xz.set_scale(scale_factor);
                    vec_xz.update_coordinate(det_picked, h, sub_xz, z_threshold_); });

                // Energy Estimate
                timer.measure("Energy Estimate", [&]()
                              { vec_xz.compute_variational_energy(); });

                // Update the stopping criterion
                if (check_convergence(vec_xz, dx_accumulated, det_picked)) {
                    early_stop = true;
                    break;
                }

                // Check overflow to prevent the container to expand automatically
                vec_xz.check_overflow();
            }
            energy       = vec_xz.get_variational_energy();

            auto xz_size = vec_xz.size_z();
            auto x_size  = vec_xz.size_x();

            auto time_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_elpse = time_end - time_start;

            // Output
            output_line((i + 1) * report_interval,
                        energy,
                        dx_accumulated,
                        x_size,
                        xz_size,
                        sub_xz.size(),
                        time_elpse.count());

            if (early_stop) break;

            // Auto adjust z_threshold
            if (z_threshold_search)
            {
                auto inc_ratio =
                    1000.0 * (xz_size - last_xz_size) / sub_xz.size() / report_interval;
                if ((inc_ratio >= 1) &&
                    (xz_size > 0.9 * vec_xz.max_load_factor() * vec_xz.capacity()))
                {
                    return -1;
                }
            }
            last_xz_size = vec_xz.size_z();
        }
        output_final(vec_xz, iterations);
        timer.print_results();

        return 0;
    }

    int solve(H                &h,
              CoordinatePick<W>      &coord_pick,
              CoordinateUpdate<H, W> &coord_update,
              W                      &vec_xz)
    {
        // Loop to search z_threshold
        while (true)
        {
            auto ierr = solve(h, coord_pick, coord_update, vec_xz,
                              z_threshold);
            if (z_threshold_search && (ierr == -1))
            {
                z_threshold *= 2;
                opt["z_threshold"] = z_threshold;
                std::cout << "z_threshold_search: z_threshold is too small. Increase "
                             "z_threhold to ";
                std::cout << std::scientific << std::setprecision(2) << z_threshold
                          << " and restart." << std::endl;

                if (opt["load_configuration"]) {
                    vec_xz.load_wavefunction(opt["load_configuration_location"]);
                } else
                    vec_xz.clear();
            }
            else
            {
                return ierr;
            }
        }
        return 0;
    }

    int solve(H &h, W &vec_xz)
    {
        auto coord_pick   = CoordinatePick<W>::init(coord_pick_str);
        auto coord_update = CoordinateUpdate<H, W>::init(coord_update_str);
        return solve(h, *coord_pick, *coord_update, vec_xz);
    }

    W solve(H &h, int num_states = 1)
    {
        auto vec_xz = W(max_wavefunction_size, max_load_factor);
        vec_xz.num_states = num_states;
        if (opt["load_configuration"]) {
            vec_xz.load_wavefunction(opt["load_configuration_location"]);
            if (opt["recalculate_z_on_load"]) {
                wff_type my_wff, sub_xz;
                auto move_wff_functor = [&](const auto& val) {
                    my_wff.push_back(val.first, val.second);
                };
                vec_xz.loop(move_wff_functor);
                vec_xz.clear();
                vec_xz.update_coordinate(my_wff, h, sub_xz,
                                         opt["z_threshold"]);
                std::cout << "Number of determinants in current wave function: "
                          << vec_xz.size_x() << " | " << vec_xz.size_z() << std::endl;
                vec_xz.compute_variational_energy();
                std::cout << "Current variational energy: " << vecmath::index(vec_xz.get_variational_energy(), 0) << std::endl;

            }
        }
        solve(h, vec_xz);
        if (opt["save_configuration"]) {
            vec_xz.dump_wavefunction(opt["save_configuration_location"]);
        }
        return vec_xz;
    }
};

template <typename H, typename W>
class CDFCISolver : public Solver<H, W>
{
public:
    using det_type = typename W::key_type;
    using wff_type = typename W::wff_type;
    using data_type = typename W::data_type;

    using Solver<H, W>::opt;

    NumericalType ground_state_energy;
    NumericalType hard_stop_energy_level;      ///< Optional hard stop condition: terminate when energy drops below this level.

    /* Constructors */
    CDFCISolver() { }

    CDFCISolver(Option &option)
    {
        opt.merge_patch(option);
        if (opt["num_coordinates"] == 1)
            opt.merge_patch({{"coordinate_pick", "gcd_grad"},
                             {"coordinate_update", "ls"}});
        else
            opt.merge_patch({{"coordinate_pick", "block_gcd_grad"},
                             {"coordinate_update", "eig"}});
        this->init_Solver();
        this->validate_option();
    }

    virtual ~CDFCISolver() {}

    void print_solver_header() const
    {
        if (this->verbose < 2)
            return;

        this->Solver<H, W>::print_solver_header();
        Option opt_   = {{"cdfci", ""}};
        opt_["cdfci"] = this->opt;
        std::cout << opt_.dump(4) << std::endl << std::endl;
        return;
    }

    NumericalType get_default_dx_accumulated() const { return 1.0; }

    void initialize(H &h, W &vec_xz, wff_type &sub_xz)
    {
        // Get reference determinant
        det_type hf = get_ref_det(h);

        // Print information
        if (this->verbose > 0)
        {
            std::cout << "CDFCI calculation" << std::endl;
            std::cout << "-----------------" << std::endl;
            std::cout << "Reference determinant occupied spin-orbitals:" << std::endl;
            for (auto &orb : hf.get_occupied_orbitals())
            {
                std::cout << orb << "  ";
            }
            std::cout << std::endl;
            std::cout << "Reference energy: " << std::fixed << std::setprecision(10)
                    << h.get_diagonal(hf) << std::endl
                    << std::endl;
        }

        // Initialize sub_xz
        sub_xz.push_back(hf, {0.0, 0.0});
    }

    det_type get_ref_det(H &h) const
    {
        std::vector<Orbital> occ_list = this->ref_det_occ;
        det_type             hf;
        if (occ_list.empty())
        {
            if (this->verbose > 1)
            {
                std::cout << "Reference determinant is not provided. Use Hartree Fock "
                             "state from FCIDUMP."
                          << std::endl
                          << std::endl;
            }
            return h.get_hartree_fock();
        }
        else
        {
            for (auto orb : occ_list)
            {
                if ((orb >= 0) && (orb < h.norb))
                {
                    hf.set_orbital(orb);
                }
            }
            if (hf.get_occupied_orbitals().size() != h.nelec)
            {
                throw std::invalid_argument(
                    "Reference determinant is invalid. The number of valid orbitals does "
                    "not match the number of electrons. Please provide " +
                    std::to_string(h.nelec) + " orbitals from 0 to " +
                    std::to_string(h.norb - 1) + " or remove the argument.");
            }
        }
        return hf;
    }

    bool check_convergence(const W &vec_xz, NumericalType &dx_accumulated,
                           const wff_type &det_picked)
    {
        // Hard Stop

        NumericalType dx_norm = 0.0;
        for (const auto &val : det_picked)
        {
            auto dx = vecmath::two_norm(vecmath::slice_first_half(val.second));
            dx_norm += dx * dx;
        }
        dx_norm = sqrt(dx_norm) * fabs(vec_xz.get_scale());
        dx_accumulated = this->stopping_dx_damping_factor * dx_accumulated
                       + dx_norm;

        // Check stopping criterion
        return (dx_accumulated < this->stopping_dx_threshold);
    }

    void output_line(const size_t iteration,
                     const data_type energy,
                     const NumericalType dx,
                     const size_t x_size,
                     const size_t z_size,
                     const size_t H_i_size,
                     const double time
                    ) const
    {
        if (this->verbose == 0) return;

        std::cout << std::setw(13) << std::left << iteration;
        if (this->shifted) {
            std::cout << std::setw(18) << std::right << std::fixed
                      << std::setprecision(10) << vecmath::index(energy, 0) + this->shift_value;
        } else {
            std::cout << std::setw(18) << std::right << std::fixed
                      << std::setprecision(10) << vecmath::index(energy, 0);
        }
        std::cout << std::setw(18) << std::scientific
                  << std::setprecision(4)  << dx;
        std::cout << std::setw(15) << x_size;
        std::cout << std::setw(15) << z_size;
        std::cout << std::setw(9)  << H_i_size;
        std::cout << std::setw(10) << std::fixed
                    << std::setprecision(2) << time;
        std::cout << std::endl;
    }

    void output_final(const W &vec_xz, const size_t &iterations)
    {
        // Store the final computed energy in the Solver class.
        ground_state_energy = vecmath::index(vec_xz.get_variational_energy(), 0);
        if (this->verbose == 0) return;

        if (this->shifted) {
            ground_state_energy += this->shift_value;
        }
        std::cout << "Final FCI Energy: " << std::setw(30)
                    << std::right << std::fixed
                    << std::setprecision(16) << ground_state_energy << std::endl;
        std::cout << "Iterations: " << iterations << std::endl;
    }

    NumericalType get_ground_state_energy() const
    {
        return ground_state_energy;
    }
};

template<typename H, typename W>
class XCDFCISolver : public CDFCISolver<H, W>
{
public:
    using det_type = typename W::key_type;
    using mapped_type = typename W::mapped_type;
    using wff_type = typename W::wff_type;
    using data_type = typename W::data_type;

    using Solver<H, W>::opt;

    int num_states;

    /* Constructors. CDFCI Constructor is called by default */
    XCDFCISolver() { }

    XCDFCISolver(Option &option, int n_states)
    {
        opt.merge_patch(option);
        opt.merge_patch({{"coordinate_pick", "gcd_grad"},
                         {"coordinate_update", "ls"}});
        num_states = n_states;
        this->init_Solver();
    }

    virtual ~XCDFCISolver() {}

    void print_solver_header() const
    {
        if (this->verbose < 2)
            return;

        this->Solver<H, W>::print_solver_header();
        Option opt_ = {{"num_states", num_states}, {"xcdfci", ""}};
        opt_["xcdfci"] = this->opt;
        std::cout << opt_.dump(4) << std::endl << std::endl;
        return;
    }

    void initialize(H &h, W &vec_xz, wff_type &sub_xz)
    {
        // Get initial determinant vector
        auto ref_det_vec = get_ref_det_vec(h, num_states);
        if (this->verbose > 0)
        {
            // Print information
            std::cout << "XCDFCI calculation" << std::endl;
            std::cout << "-----------------" << std::endl;
            std::cout << "Initial determinant occupied spin-orbitals | Energy" << std::endl;
            for (det_type &det : ref_det_vec)
            {
                for (auto &orb : det.get_occupied_orbitals())
                {
                    std::cout << orb << "  ";
                }
                std::cout << " | ";
                std::cout << std::fixed << std::setprecision(10)
                        << h.get_diagonal(det) << std::endl;
            }
        }

        // Initialize sub_xz
        wff_type ref_det_wff;
        for (auto it = 0; it < num_states; ++it)
        {
            mapped_type dx = {};
            dx[it] = 1.0;
            ref_det_wff.push_back(ref_det_vec[it], dx);
        }

        vec_xz.update_coordinate(ref_det_wff, h, sub_xz, this->z_threshold);
        vec_xz.print_xx_xz();
    }

    std::vector<det_type> get_ref_det_vec(H &h, const int num_states) const
    {
        // return num_states determinants with lowest diagonal energies
        // among all determinants that are connected to ref_det.
        det_type hf = CDFCISolver<H, W>::get_ref_det(h);
        auto column = h.get_column(hf);
        for (auto& entry: column)
        {
            auto& det = entry.first;
            entry.second = h.get_diagonal(det);
        }
        std::sort (column.begin(), column.end(),
                [](const std::pair<det_type, double> & a,
                    const std::pair<det_type, double> & b) -> bool
                    {
                        return a.second < b.second;
                    });
        std::vector<det_type> ref_vec;
        for (int it = 0; it < num_states; ++it)
            ref_vec.push_back(column[it].first);
        return ref_vec;
    }

    void output_line(const size_t iteration,
                     const data_type energy,
                     const NumericalType dx,
                     const size_t x_size,
                     const size_t z_size,
                     const size_t H_i_size,
                     const double time
                    ) const
    {
        if (this->verbose == 0) return;

        std::cout << std::setw(13) << std::left << iteration;
        std::cout << std::setw(18) << std::right << std::fixed
                  << std::setprecision(10) << vecmath::index(energy, 0);
        std::cout << std::setw(18) << std::scientific
                  << std::setprecision(4)  << dx;
        std::cout << std::setw(15) << x_size;
        std::cout << std::setw(15) << z_size;
        std::cout << std::setw(9)  << H_i_size;
        std::cout << std::setw(10) << std::fixed
                  << std::setprecision(2) << time;
        std::cout << std::endl;

        for (int i = 1; i < num_states; i++)
        {
            std::cout << std::setw(13) << std::left << " ";
            std::cout << std::setw(18) << std::right << std::fixed
                      << std::setprecision(10) << vecmath::index(energy, i);
            std::cout << std::endl;
        }

    }

    void output_final(const W &vec_xz, const size_t &iterations)
    {
        auto energy = vec_xz.get_variational_energy();
        this->ground_state_energy = vecmath::index(energy, 0);

        if (this->shifted)
        {
            this->ground_state_energy += this->shift_value;
        }

        if (this->verbose == 0) return;

        std::cout << "Final FCI Energies:" << std::endl;
        for (int i = 0; i < num_states; ++i)
        {
            auto state_energy = vecmath::index(energy, i);
            if (this->shifted)
            {
                state_energy += this->shift_value;
            }
            std::cout << "  State " << i << ": " << std::setw(30)
                      << std::right << std::fixed << std::setprecision(16)
                      << state_energy << std::endl;
        }
        std::cout << "Iterations: " << iterations << std::endl;
    }

    NumericalType get_ground_state_energy() const
    {
        return this->ground_state_energy;
    }
};

template<typename H, typename W>
std::unique_ptr<Solver<H,W>> Solver<H,W>::init(Option &option)
{
    if (!option.contains("type") || !option["type"].is_string()) {
        throw std::invalid_argument("Missing or invalid 'type' field in solver option.");
    }

    std::string type = option["type"];

    if (type == "cdfci")
    {
        if (!option.contains("cdfci") || !option["cdfci"].is_object()) {
            throw std::invalid_argument("Missing or invalid 'cdfci' configuration in solver option.");
        }
        return std::make_unique<CDFCISolver<H, W>>(option["cdfci"]);
    }

    if (type == "xcdfci")
    {
        if (!option.contains("xcdfci") || !option["xcdfci"].is_object()) {
            throw std::invalid_argument("Missing or invalid 'xcdfci' configuration in solver option.");
        }
        return std::make_unique<XCDFCISolver<H, W>>(option["xcdfci"], option["num_states"]);
    }

    throw std::invalid_argument("Invalid solver type: " + type);
}
#endif
