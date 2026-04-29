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

#ifndef PROGRAM_H
#define PROGRAM_H

#include <chrono>
#include <iostream>

#include "solver.h"
#include "rdm.h"
#include "perturbation.h"
#include "optimizer.h"

template <typename Func>
void time_block(const std::string& tag, Func&& f)
{
    auto t0 = std::chrono::high_resolution_clock::now();
    f();
    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "[Timer] " << tag << " took " << dt << " s\n";
}

template <int N>
class CDFCIProgram
{
public:
    using ham_type    = Hamiltonian<N>;
    using det_type    = typename ham_type::det_type;
    using key_type    = det_type;
    using mapped_type = std::array<NumericalType, 2>;
    using equal_type  = DeterminantEqual<N>;
#ifdef CDFCI_SOLVER_SERIAL
    using hash_type      = DeterminantHash<N>;
    using container_type = ContainerRobinhood<key_type, mapped_type, hash_type, equal_type>;
#else
    using hash_type      = DeterminantHashRobinhood<N>;
    using container_type = ContainerCuckoo<key_type, mapped_type, hash_type, equal_type>;
#endif
    using wf_type     = WaveFunction<container_type>;
    using solver_type = Solver<ham_type, wf_type>;

    Option opt;

    std::unique_ptr<ham_type> ham_ptr;

    std::unique_ptr<solver_type> solver_ptr;

    // wf_type vec_xz;

    std::unique_ptr<RDM<wf_type>> rdm_ptr;
    std::unique_ptr<Perturbation<ham_type,wf_type>> ptb_ptr;

    CDFCIProgram() { }

    CDFCIProgram(Option option)
    {
        opt = option;
        ham_ptr = Hamiltonian<N>::init(opt["hamiltonian"]);
        solver_ptr = Solver<ham_type, wf_type>::init(opt["solver"]);
    }

    CDFCIProgram(Option option, Fcidump &fci)
    {
        opt = option;

        ham_ptr = std::make_unique<HamiltonianMolecule<N>>(fci);

        solver_ptr = Solver<ham_type,wf_type>::init(opt["solver"]);
    }

    ~CDFCIProgram() {}

    NumericalType run()
    {
        det_type::constuct_masks();
        wf_type vec_xz;
        time_block("Solver", [&]{
            vec_xz = solver_ptr->solve(*ham_ptr);
        });

        if (opt.contains("perturbation"))
        {
            time_block("Perturbation energy", [&]{
                using wff_type = typename wf_type::wff_type;
                wff_type my_wff, sub_xz;
                auto move_wff_functor = [&](const auto& val) {
                    my_wff.push_back(val.first, val.second);
                };
                vec_xz.loop(move_wff_functor);
                vec_xz.clear();
                vec_xz.update_coordinate(my_wff, *ham_ptr, sub_xz,
                                         opt["perturbation"].value("z_threshold", 0.0));
                std::cout << "Number of determinants in perturbation wave function: "
                          << vec_xz.size_x() << " | " << vec_xz.size_z() << std::endl;
                ptb_ptr = std::make_unique<Perturbation<ham_type,wf_type>>();
                ptb_ptr->compute_perturbation_energy(vec_xz, *ham_ptr);
            });
        }

        if (opt.contains("rdm"))
        {
            time_block("RDM computation", [&]{
                rdm_ptr = std::make_unique<RDM<wf_type>>(
                    opt["rdm"], static_cast<int>(ham_ptr->norb / 2), ham_ptr->nelec);
                rdm_ptr->computeRDM(vec_xz);
                rdm_ptr->dumpFile();
            });
        }

        return solver_ptr->get_ground_state_energy();
    }

    void print_header()
    {
        ham_ptr->print_header();
        solver_ptr->print_solver_header();
    }
};

template <int N, int NSTATES>
class XCDFCIProgram
{
public:
    using ham_type    = Hamiltonian<N>;
    using det_type    = typename ham_type::det_type;
    using key_type    = det_type;
    using mapped_type = std::array<NumericalType, 2*NSTATES>;
    using equal_type  = DeterminantEqual<N>;
#ifdef CDFCI_SOLVER_SERIAL
    using hash_type      = DeterminantHash<N>;
    using container_type = ContainerRobinhood<key_type, mapped_type, hash_type, equal_type>;
#else
    using hash_type      = DeterminantHashRobinhood<N>;
    using container_type = ContainerCuckoo<key_type, mapped_type, hash_type, equal_type>;
#endif
    using wf_type = WaveFunction<container_type, NSTATES>;
    using solver_type = Solver<ham_type, wf_type>;

    Option opt;

    std::unique_ptr<ham_type> ham_ptr;

    std::unique_ptr<solver_type> solver_ptr;

    XCDFCIProgram(Option option, int num_states)
    {
        opt = option;
        opt["solver"]["num_states"] = num_states;
        ham_ptr = Hamiltonian<N>::init(opt["hamiltonian"]);
        solver_ptr = Solver<ham_type, wf_type>::init(opt["solver"]);
    }

    ~XCDFCIProgram() {}

    void run()
    {
        det_type::constuct_masks();
        solver_ptr->solve(*ham_ptr, opt["solver"]["num_states"]);
    }

    void print_header()
    {
        ham_ptr->print_header();
        solver_ptr->print_solver_header();
    }
};


template <int N>
class OptOrbFCIProgram
{
public:
    Option  opt;
    Fcidump fci;
    int     init_iter;
    int     max_iter;
    int     norb_orig; // norb. spatial orbitals.
    int     norb_comp;

    Eigen::VectorXd orb_energies;
    NumericalType   zero_int;
    Matrix          one_int; // One-electron integrals (norb × norb) index: (i, a)
    Tensor4         two_int; // Two-electron integrals (norb⁴)       index: (i, a, j, b)

    inline Orbital alpha(const int i) const { return 2 * i; }
    inline Orbital beta(const int i) const { return 2 * i + 1; }

public:
    OptOrbFCIProgram(Option option)
    {
        opt = option;
        try
        {
            std::string fcidump_path = optutil::fcidump_path(opt);
            fci                      = Fcidump(fcidump_path);
        }
        catch (const Option::out_of_range &e)
        {
            throw std::runtime_error("Error: Missing 'hamiltonian.molecule.fcidump_path' in input JSON.");
        }
        catch (const Option::type_error &e)
        {
            throw std::runtime_error("Error: 'fcidump_path' should be a string.");
        }

        init_iter = optutil::init_iteration(opt);
        max_iter  = optutil::max_iteration(opt);
        norb_orig = fci.norb / 2;                          // spatial orbitals before compression
        norb_comp = optutil::num_orbitals_compressed(opt); // after compression
        std::cout << "norb_orig = " << norb_orig << " norb_comp = " << norb_comp << std::endl;

        // transform orb_energies → Eigen::VectorXd
        orb_energies.resize(norb_orig);
        for (int i = 0; i < norb_orig; ++i)
            orb_energies(i) = fci.orbital_energy.count(alpha(i)) ? fci.orbital_energy.at(alpha(i)) : 0.0;

        std::cout << "orb_energies" << std::endl;

        zero_int = fci.core_energy;

        // transform one_body_integral → Matrix
        one_int.resize(norb_orig, norb_orig);
        for (int p = 0; p < norb_orig; ++p)
            for (int q = 0; q < norb_orig; ++q)
            {
                auto x        = fci.one_body_integral[alpha(p)][alpha(q)];
                one_int(p, q) = x;
                one_int(q, p) = x;
            }

        // transform two_body_integral → Eigen::Tensor<value_type, 4>
        two_int.resize(norb_orig, norb_orig, norb_orig, norb_orig);
        two_int.setZero();
        for (const auto &[pq, inner] : fci.two_body_integral)
        {
            int p = pq.first;
            int q = pq.second;
            if (p % 2 == 1 || q % 2 == 1)
                continue;
            p /= 2;
            q /= 2;
            for (const auto &[rs, val] : inner)
            {
                int r = rs.first;
                int s = rs.second;
                if (r % 2 == 1 || s % 2 == 1)
                    continue;
                r /= 2;
                s /= 2;
                two_int(p, r, q, s) = val;
                two_int(q, s, p, r) = val;
            }
        }
    }

    ~OptOrbFCIProgram() {}

    NumericalType run()
    {
        NumericalType energy;

        Matrix        U(norb_orig, norb_comp);
        NumericalType zero_rdm;
        Matrix        one_rdm(norb_comp, norb_comp);
        Tensor4       two_rdm(norb_comp, norb_comp, norb_comp, norb_comp);

        U.setZero();

        if (init_iter == 0)
        {
            std::vector<std::pair<Orbital, NumericalType>> Eidx =
                get_sorted_orbital_energy(fci.orbital_energy, norb_orig);

            for (int i = 0; i < norb_comp; i++)
                U(Eidx[alpha(i)].first / 2, i) = 1.0;

            fci = rotate_basis(fci, U);
            // run cdfci
            CDFCIProgram<N> cdfci(opt, fci);
            energy = cdfci.run();
            cdfci.rdm_ptr->output_rdm(zero_rdm, one_rdm, two_rdm);

            init_iter = 1;
        }
        else
        {
            // load history file
        }

        Optimizer myOpt(opt["optimal_orbitals"]["optimizer"], zero_int, one_int, two_int);
        for (int iter = init_iter; iter < max_iter; iter++)
        {
            myOpt.update_rdm(zero_rdm, one_rdm, two_rdm);

            // try optimize U with random initial guess
            bool rand_iter_success_flag = false;
            for (int randiter = 0; randiter < myOpt.opt["max_rand_iter"]; randiter++)
            {
                Matrix V = myOpt.add_Gaussian_noise(U, 0.1);
                V        = myOpt.project_to_manifold(V);
                V        = myOpt.optimize_U(V);
                if (myOpt.final_value < energy)
                {
                    std::cout << "Rand initialization success!" << std::endl;
                    U                      = V;
                    rand_iter_success_flag = true;
                    break;
                }
            }
            if (!rand_iter_success_flag)
                U = myOpt.optimize_U(U);

            fci = rotate_basis(fci, U);
            // run cdfci
            CDFCIProgram<N> cdfci(opt, fci);
            NumericalType   new_energy = cdfci.run();

            // check convergence for optorbfci
            if (fabs(new_energy - energy) < opt["optimal_orbitals"]["tol"])
                return new_energy;

            energy = new_energy;
            cdfci.rdm_ptr->output_rdm(zero_rdm, one_rdm, two_rdm);
        }

        return energy;
    }

    Fcidump rotate_basis(const Fcidump &orig_fci, const Matrix &U)
    {
        Fcidump new_fci(norb_comp * 2, orig_fci.nelec, orig_fci.ms2,
                        std::vector<int>(norb_comp * 2, 1), orig_fci.uhf);

        new_fci.core_energy = orig_fci.core_energy;

        // Rotate orbitals
        Eigen::VectorXd new_orb_energies = (U.transpose() * orb_energies.asDiagonal() * U).diagonal();
        Matrix          new_one_int      = contracter::rotate_one_int(one_int, U);
        Tensor4         new_two_int      = contracter::rotate_two_int(two_int, U);

        // Write to FCI
        new_fci.one_body_integral.resize(norb_comp * 2);
        for (auto &v : new_fci.one_body_integral)
            v.resize(norb_comp * 2, 0.0);

        for (int i = 0; i < norb_comp; ++i)
        {
            new_fci.orbital_energy.insert({alpha(i), new_orb_energies(i)});
            new_fci.orbital_energy.insert({beta(i), new_orb_energies(i)});

            for (int a = 0; a <= i; ++a)
            {
                auto integral                                 = new_one_int(i, a);
                new_fci.one_body_integral[alpha(i)][alpha(a)] = integral;
                new_fci.one_body_integral[beta(i)][beta(a)]   = integral;
                new_fci.one_body_integral[alpha(a)][alpha(i)] = integral;
                new_fci.one_body_integral[beta(a)][beta(i)]   = integral;

                for (int j = 0; j < norb_comp; ++j)
                    for (int b = 0; b <= j; ++b)
                    {
                        auto x = new_two_int(i, a, j, b);
                        new_fci.rhf_insert_two_body_(i + 1, j + 1, a + 1, b + 1, x);
                        new_fci.rhf_insert_two_body_(i + 1, b + 1, a + 1, j + 1, x);
                        new_fci.rhf_insert_two_body_(a + 1, j + 1, i + 1, b + 1, x);
                        new_fci.rhf_insert_two_body_(a + 1, b + 1, i + 1, j + 1, x);
                    }
            }
        }

        return new_fci;
    }
};

#endif