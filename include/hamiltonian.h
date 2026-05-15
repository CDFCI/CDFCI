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

#ifndef CDFCI_HAMILTONIAN_H
#define CDFCI_HAMILTONIAN_H 1

#include <algorithm>
#include <memory>
#include <cassert>

#include "option.h"
#include "fcidump.h"
#include "determinant.h"

template <int N = 1>
class Hamiltonian
{
public:
    using det_type   = Determinant<N>;
    using value_type = NumericalType;
    using Column     = std::vector<std::pair<Determinant<N>, NumericalType>>;

    Option opt = {
        {"type", ""}, // required, no default value
    };

    int  norb;
    int  nelec;
    int  ms2;
    bool uhf;

    double diagonal_shift = 0.0; // Energy shift for the Hamiltonian, positive (H_new = H - shift).

    void apply_shift(double shift) { diagonal_shift = shift; }
    double get_shift() const { return diagonal_shift; }

    Hamiltonian() : norb{0}, nelec{0}, ms2{0}, uhf{false} {}

    // Defined outside
    static std::unique_ptr<Hamiltonian<N>> init(Option option);
    static int                             read_norb(Option option);

    virtual void print_header() const
    {
        std::cout << "Hamiltonian information" << std::endl;
        std::cout << "-------------------" << std::endl;
        std::cout << "Number of electrons: " << nelec << std::endl;
        std::cout << "Number of spin-orbitals: " << norb << std::endl;
        std::cout << "Ms2: " << ms2 << std::endl;
        std::string hf_type = uhf ? "unrestricted" : "restricted";
        std::cout << "Hartree Fock: " + hf_type << std::endl
                  << std::endl;
        return;
    }

    virtual Column        get_column(Determinant<N> &det) const                                                                             = 0;
    virtual Column        get_column_det(Determinant<N> &det) const                                                                         = 0;
    virtual void          get_column_val(Determinant<N> &det, typename Column::iterator col_begin, typename Column::iterator col_end) const = 0;
    virtual Column        get_column_serial_part(Determinant<N> &det) const                                                                 = 0;
    virtual Column        get_column_parallel_part(Determinant<N> &det, int tid,
                                                   int nthreads) const                                                                      = 0;
    virtual NumericalType get_diagonal(Determinant<N> &det) const                                                                           = 0;
    virtual NumericalType get_entry(Determinant<N> &det1, Determinant<N> &det2) const                                                       = 0;

    // virtual Column get_column(DeterminantDecoded<N>& det) const = 0;
    virtual Determinant<N> get_hartree_fock() const                                     = 0;
    virtual Column         get_column_diagonal_single(DeterminantDecoded<N> &det) const = 0;
    virtual void           get_double_excitation(DeterminantDecoded<N> &det, Column &result,
                                                 int idx_i, int idx_j) const            = 0;

    virtual ~Hamiltonian() {}
};

template <int N = 1>
class HamiltonianMolecule : public Hamiltonian<N>
{
    using OrbitalList = std::vector<Orbital>;
    using typename Hamiltonian<N>::det_type;
    using typename Hamiltonian<N>::Column;
    using Hamiltonian<N>::norb;
    using Hamiltonian<N>::nelec;
    using Hamiltonian<N>::ms2;
    using Hamiltonian<N>::uhf;

private:
    // Double excitation data structure.
    using DoubleExcitationEntry = std::pair<OrbitalPair, NumericalType>;
    std::vector<std::vector<DoubleExcitationEntry>> double_excitation; // norb * norb

    // Single excitation data structure.
    struct SingleExcitationEntry
    {
        Orbital                    a;
        NumericalType              one_body_int; // <i|f|a>
        std::vector<NumericalType> two_body_int; // <ik||ak>

        SingleExcitationEntry(int norb) : a{0}, one_body_int{0}
        {
            two_body_int.resize(norb);
        }
    };
    std::vector<std::vector<SingleExcitationEntry>> single_excitation;

    // Diagonal data structure.
    struct DiagonalData
    {
        std::vector<NumericalType> one_body_integral; // <i|f|i>, size = norb
        std::vector<NumericalType> two_body_integral; // <ij||ij>, size = norb * norb
    };
    DiagonalData diagonal;

    // Other members
    NumericalType                              core_energy;
    std::unordered_map<Orbital, NumericalType> orbital_energy;

    int verbose = 2;

public:
    Option opt = {{"fcidump_path", "FCIDUMP"}, {"threshold", 1e-14}, {"verbose", 2}};

    std::string   fcidump_path;
    NumericalType threshold;

    using Hamiltonian<N>::Hamiltonian;

    HamiltonianMolecule() { init_HamiltonianMolecule(); }

    HamiltonianMolecule(Option option)
    {
        opt.merge_patch(option);
        init_HamiltonianMolecule();
    }

    HamiltonianMolecule(Fcidump &fci, NumericalType threshold_ = 0.0)
    {
        init_HamiltonianMolecule(fci, threshold_);
    }

    void init_HamiltonianMolecule()
    {
        fcidump_path = opt["fcidump_path"];
        threshold    = opt["threshold"];
        verbose      = opt["verbose"];

        Fcidump fci(fcidump_path);
        init_HamiltonianMolecule(fci, threshold);
    }

    void init_HamiltonianMolecule(Fcidump &fci, NumericalType threshold_ = 0.0)
    {
        norb        = fci.norb;
        nelec       = fci.nelec;
        ms2         = fci.ms2;
        uhf         = fci.uhf;
        core_energy = fci.core_energy;

        validate_option();

        construct_double_excitation(fci, threshold_);
        construct_single_excitation(fci, threshold_);
        construct_diagonal(fci);
        orbital_energy = fci.orbital_energy;
    }

    void validate_option()
    {
        if (norb <= 0)
        {
            throw std::invalid_argument("Number of spin-orbitals norb (" +
                                        std::to_string(norb) + ") shoule be positive.");
        }
        if (nelec <= 0)
        {
            throw std::invalid_argument("Number of electrons nelec (" +
                                        std::to_string(norb) + ") shoule be positive.");
        }
        if (norb < nelec)
        {
            throw std::invalid_argument(
                "Number of spin-orbitals norb (" + std::to_string(norb) +
                ") shoule be no less than the number of electrons nelec (" +
                std::to_string(nelec) + ")");
        }
        if ((nelec + ms2) % 2 != 0)
        {
            throw std::invalid_argument("Required nelec (" + std::to_string(nelec) +
                                        ") and ms2 (" + std::to_string(ms2) +
                                        ") are not compatible.");
        }
    }

    static int read_norb(Option &option)
    {
        return Fcidump::read_norb(option["fcidump_path"]);
    }

    void print_header() const
    {
        if (verbose < 2)
            return;

        this->Hamiltonian<N>::print_header();
        Option opt_      = {{"molecule", ""}};
        opt_["molecule"] = opt;
        std::cout << opt_.dump(4) << std::endl
                  << std::endl;
        return;
    }

    int index(const Orbital i, const Orbital j) const { return i * norb + j; }

    // Construct double excitation <ij||ab> = <ij|g|ab> - <ij|g|ba>.
    // Only store i < j. (In Fcidump, i <= j).
    void construct_double_excitation(const Fcidump &fci, double threshold)
    {
        // Initialize
        double_excitation.resize(norb * norb);

        for (const auto &ij_entry : fci.two_body_integral)
        {
            auto ij = ij_entry.first;
            auto i  = ij.first;
            auto j  = ij.second;
            if (i == j)
            {
                continue; // <ii||ab> is always zero.
            }

            for (const auto &ab_entry : ij_entry.second) // key = (a, b), value = integral
            {
                auto ab            = ab_entry.first;
                auto a             = ab.first;
                auto b             = ab.second;
                auto integral_ijab = ab_entry.second;

                // opposite string, simply <ij|g|ab>. Note a > b and a < b are both
                // allowed. i % 2 = 0 if i is alpha. i % 2 = 1 if i is beta.
                if ((i % 2) ^ (j % 2))
                {
                    if (fabs(integral_ijab) > threshold)
                    {
                        double_excitation[index(i, j)].push_back(ab_entry);
                    }
                }
                // same string, <ij||ab> = <ij|g|ab> - <ij|g|ba>
                // only store a < b
                else if (a < b)
                {
                    auto integral =
                        integral_ijab -
                        fci.get_two_body_integral(
                            i, j, b, a); // Q: Why can not const fci? A: fixed. Cannot use
                                         // [] for unordered_map.
                    std::pair<OrbitalPair, NumericalType> ab_int{ab, integral};
                    if (fabs(integral) > threshold)
                    {
                        double_excitation[index(i, j)].push_back(ab_int);
                    }
                }
                else if ((b < a) && (fci.two_body_integral.at(ij).find({b, a}) ==
                                     fci.two_body_integral.at(ij).end()))
                {
                    auto integral = -integral_ijab; // <ij|g|min(a,b) max(a,b)> = 0 but
                                                    // <ij|g|max(a,b) min(a,b)> != 0
                    std::pair<OrbitalPair, NumericalType> ab_int{{b, a}, integral};
                    if (fabs(integral) > threshold)
                    {
                        double_excitation[index(i, j)].push_back(ab_int);
                    }
                }
            }
        }
        // Sort
        const auto comp = [](auto &x, auto &y)
        {
            return fabs(x.second) > fabs(y.second);
        };
        for (auto &ij_column : double_excitation)
            std::sort(ij_column.begin(), ij_column.end(), comp);

        return;
    }

    // Construct single excitation
    void construct_single_excitation(const Fcidump &fci, double threshold)
    {
        // Initialization
        single_excitation.resize(norb);

        // for (const auto& i_entry : fci.one_body_integral)
        for (auto i = 0; i < norb; ++i)
        {
            // auto i = i_entry.first;
            // for (const auto& a_entry : i_entry.second)
            for (auto a = 0; a < norb; ++a)
            {
                SingleExcitationEntry entry(norb);
                entry.a            = a;
                entry.one_body_int = fci.get_one_body_integral(i, a);

                auto max_val = fabs(entry.one_body_int);

                // Find all double integrals <ik|g|ak> and <ik||ak> = <ik|g|ak> -
                // <ik|g|ka>
                for (auto k = 0; k < norb; ++k)
                {
                    NumericalType integral;
                    if (i <= k)
                    {
                        integral = fci.get_two_body_integral(i, k, a, k) -
                                   fci.get_two_body_integral(i, k, k, a);
                    }
                    else
                    {
                        integral = fci.get_two_body_integral(k, i, k, a) -
                                   fci.get_two_body_integral(k, i, a, k);
                    }

                    entry.two_body_int[k] = integral;
                    max_val               = std::max(max_val, fabs(integral));
                }

                // Store the excitation if large enough.
                if (max_val > threshold)
                {
                    single_excitation[i].push_back(entry);
                }
            }
        }
        // Sort
        const auto comp = [](auto &x, auto &y)
        {
            return fabs(x.one_body_int) > fabs(y.one_body_int);
        };
        for (auto &i_column : single_excitation)
            std::sort(i_column.begin(), i_column.end(), comp);

        return;
    }

    // Size of one_body_integral = norb.
    // one_body_integral[i] = <i|f|i>
    // Size of two_body_integral = norb * norb.
    // two_body_integral{{i, j}} = <ij||ij>, i < j.
    // The lower diagonal is not used.
    void construct_diagonal(const Fcidump &fci)
    {
        // Reserve space
        diagonal.one_body_integral.resize(norb);
        diagonal.two_body_integral.resize(norb * norb);

        // One-body integral <i|f|i>, diagonal
        for (auto i = 0; i < norb; ++i)
        {
            diagonal.one_body_integral[i] = fci.get_one_body_integral(i, i);
        }

        // Two-body integral <ij||ij> (same string) and <ij|g|ij> (opposite string)
        for (auto i = 0; i < norb; ++i)
            for (auto j = i + 1; j < norb; ++j)
            {
                diagonal.two_body_integral[index(i, j)] =
                    fci.get_two_body_integral(i, j, i, j) -
                    fci.get_two_body_integral(i, j, j, i);
            }

        return;
    }

    // int same_index(Orbital i, Orbital j) {return j * norb + i;}
    // int oppo_index(Orbital i, Orbital j) {return i * norb + j;}

    void get_diagonal(DeterminantDecoded<N> &det, Column &result) const
    {
        auto value = get_diagonal(det);
        result.push_back({det, value});
    }

    NumericalType get_diagonal(DeterminantDecoded<N> &det) const
    {
        auto occupied_orbitals = det.get_occupied_orbitals();
        return get_diagonal(occupied_orbitals);
    }

    NumericalType get_diagonal(Determinant<N> &det) const
    {
        DeterminantDecoded<N> det_decoded(det);
        return get_diagonal(det_decoded);
    }

    NumericalType get_diagonal(OrbitalList &occupied_orbitals) const
    {
        NumericalType value = 0;

        for (auto idx_i = 0; idx_i < occupied_orbitals.size(); ++idx_i)
        {
            auto i = occupied_orbitals[idx_i];
            // One body integral
            value += diagonal.one_body_integral[i];
            // Two body integral, i < j
            for (auto idx_j = idx_i + 1; idx_j < occupied_orbitals.size(); ++idx_j)
            {
                auto j = occupied_orbitals[idx_j];
                value += diagonal.two_body_integral[index(i, j)];
            }
        }
        value += core_energy;
        value -= this->diagonal_shift;
        return value;
    }

    NumericalType get_entry(Determinant<N> &det1, Determinant<N> &det2) const
    {
        NumericalType  value    = 0.0;
        Determinant<N> xor_det  = det1 ^ det2;
        int            popcount = xor_det.popcount();
        if (popcount == 0)
            return get_diagonal(det1);
        else if (popcount == 2)
        { // single excitation from i to a
            OrbitalList i_orbitals = (xor_det & det1).get_occupied_orbitals();
            OrbitalList a_orbitals = (xor_det & det2).get_occupied_orbitals();
            Orbital     i          = i_orbitals[0];
            for (auto &entry_a : single_excitation[i])
            {
                Orbital a = entry_a.a;
                if (a == a_orbitals[0])
                {
                    // <i|f|a>
                    value = entry_a.one_body_int;
                    // <ik||ak>
                    for (auto k : det1.get_occupied_orbitals())
                        value += entry_a.two_body_int[k];
                    // Calculate sign
                    Determinant<N> new_det(det1);
                    new_det.clear_orbital(i);
                    if (new_det.parity(i, a))
                        value = -value;
                    break;
                }
            }
        }
        else if (popcount == 4)
        { // double excitation from i, j to a, b
            OrbitalList ij_orbitals = (xor_det & det1).get_occupied_orbitals();
            OrbitalList ab_orbitals = (xor_det & det2).get_occupied_orbitals();
            Orbital     i           = ij_orbitals[0];
            Orbital     j           = ij_orbitals[1];
            for (auto &entry : double_excitation[index(i, j)])
            {
                Orbital a = entry.first.first;
                Orbital b = entry.first.second;
                if ((a == ab_orbitals[1] && b == ab_orbitals[0]) ||
                    (a == ab_orbitals[0] && b == ab_orbitals[1]))
                {
                    value = entry.second;
                    Determinant<N> new_det(det1);
                    // Calculate sign
                    new_det.clear_orbital(i);
                    auto sign1 = new_det.parity(i, a);
                    new_det.set_orbital(a);
                    new_det.clear_orbital(j);
                    auto sign2 = new_det.parity(j, b);
                    new_det.set_orbital(b);
                    auto sign = sign1 ^ sign2;
                    if (sign)
                        value = -value;
                    break;
                }
            }
        }
        return value;
    }

    Column get_single_excitation(DeterminantDecoded<N> &det) const
    {
        Column result;
        get_single_excitation(det, result);
        return result;
    }

    void get_single_excitation(DeterminantDecoded<N> &det, Column &result) const
    {
        // Single excitation from orbital i
        for (auto i : det.occupied_orbitals)
        {
            // Loop over all possible excitation orbital a
            for (auto &entry_a : single_excitation[i])
            {
                auto a = entry_a.a;
                // The excitation exists if a is not occupied.
                if (!det.is_occupied(a))
                {
                    // <i|f|a>
                    NumericalType value = entry_a.one_body_int;
                    // <ik||ak>
                    for (auto k : det.occupied_orbitals)
                        value += entry_a.two_body_int[k];
                    // value is the Hamiltonian. Construct the new Determinant.
                    Determinant<N> new_det(det);
                    new_det.clear_orbital(i);
                    // Calculate sign
                    auto sign = new_det.parity(i, a);
                    if (sign)
                        value = -value;

                    new_det.set_orbital(a);
                    result.push_back({new_det, value});
                }
            }
        }
        return;
    }

    Column get_double_excitation(DeterminantDecoded<N> &det) const
    {
        Column result;
        get_double_excitation(det, result);
        return result;
    }

    void get_double_excitation(DeterminantDecoded<N> &det, Column &result) const
    {
        // Double excitation from orbitals (i, j).
        for (auto idx_i = 0; idx_i < det.occupied_orbitals.size(); ++idx_i)
        {
            // auto i = det.occupied_orbitals[idx_i];
            for (auto idx_j = idx_i + 1; idx_j < det.occupied_orbitals.size(); ++idx_j)
            {
                // auto j = det.occupied_orbitals[idx_j];
                get_double_excitation(det, result, idx_i, idx_j);
            }
        }

        return;
    }

    // Get all double excitations from occupied orbital index (idx_i, idx_j).
    void get_double_excitation(DeterminantDecoded<N> &det, Column &result, int idx_i,
                               int idx_j) const
    {
        auto i = det.occupied_orbitals[idx_i];
        auto j = det.occupied_orbitals[idx_j];
        // Loop over all possible orbitals a, b
        for (auto &entry : double_excitation[index(i, j)])
        {
            auto a = entry.first.first;
            auto b = entry.first.second;
            // The excitation exists if a, b are not occupied
            if (!det.is_occupied(a) & !det.is_occupied(b))
            {
                NumericalType  value = entry.second;
                Determinant<N> new_det(det);
                // Calculate sign
                new_det.clear_orbital(i);
                auto sign1 = new_det.parity(i, a);
                new_det.set_orbital(a);
                new_det.clear_orbital(j);
                auto sign2 = new_det.parity(j, b);
                new_det.set_orbital(b);
                auto sign = sign1 ^ sign2;
                if (sign)
                    value = -value;

                result.push_back({new_det, value});
            }
        }
        return;
    }

    Column get_column(Determinant<N> &det) const
    {
        DeterminantDecoded<N> det_decoded(det);
        Column                result;
        get_diagonal(det_decoded, result);
        get_single_excitation(det_decoded, result);
        get_double_excitation(det_decoded, result);
        return result;
    }

    /* This function only returns the determinants, without the values.*/
    Column get_column_det(Determinant<N> &det) const
    {
        DeterminantDecoded<N> det_decoded(det);
        Column                result;
        result.push_back({det, 0.0});
        // Single excitation from orbital i
        for (auto i : det_decoded.occupied_orbitals)
        {
            // Loop over all possible excitation orbital a
            for (auto &entry_a : single_excitation[i])
            {
                auto a = entry_a.a;
                // The excitation exists if a is not occupied.
                if (!det.is_occupied(a))
                {
                    Determinant<N> new_det(det);
                    new_det.clear_orbital(i);
                    new_det.set_orbital(a);
                    result.push_back({new_det, 0.0});
                }
            }
        }
        // Double excitation from orbitals (i, j).
        for (auto idx_i = 0; idx_i < det_decoded.occupied_orbitals.size(); ++idx_i)
        {
            auto i = det_decoded.occupied_orbitals[idx_i];
            for (auto idx_j = idx_i + 1; idx_j < det_decoded.occupied_orbitals.size(); ++idx_j)
            {
                auto j = det_decoded.occupied_orbitals[idx_j];
                // Loop over all possible orbitals a, b
                for (auto &entry : double_excitation[index(i, j)])
                {
                    auto a = entry.first.first;
                    auto b = entry.first.second;
                    // The excitation exists if a, b are not occupied
                    if (!det.is_occupied(a) & !det.is_occupied(b))
                    {
                        Determinant<N> new_det(det);
                        new_det.clear_orbital(i);
                        new_det.clear_orbital(j);
                        new_det.set_orbital(a);
                        new_det.set_orbital(b);
                        result.push_back({new_det, 0.0});
                    }
                }
            }
        }
        return result;
    }

    /* This function computes the value of given determinants in one column.*/
    void get_column_val(Determinant<N> &det, typename Column::iterator begin, typename Column::iterator end) const
    {
        assert(std::distance(begin, end) >= 0);
        for (auto it = begin; it != end; it++)
            it->second = get_entry(det, it->first);
    }

    Column get_column_diagonal_single(DeterminantDecoded<N> &det) const
    {
        Column result;
        get_diagonal(det, result);
        get_single_excitation(det, result);
        return result;
    }

    Column get_column_serial_part(Determinant<N> &det) const
    {
        DeterminantDecoded<N> det_decoded(det);
        Column                result;
        get_diagonal(det_decoded, result);
        get_single_excitation(det_decoded, result);

        return result;
    }

    Column get_column_parallel_part(Determinant<N> &det, int tid, int nthreads) const
    {
        DeterminantDecoded<N> det_decoded(det);
        Column                result;

        int n_excitations = nelec * (nelec - 1) / 2;
        // Compute the range the thread needs to compute
        // 0 <= tid < nthreads
        int idx_start = tid * n_excitations / nthreads;
        int idx_end   = (tid + 1) * n_excitations / nthreads;

        for (auto idx = idx_start; idx < idx_end; ++idx)
        {
            int idx_j = 0.5 + sqrt(0.25 + 2 * idx);
            int idx_i = idx - (idx_j - 1) * idx_j / 2;
            get_double_excitation(det_decoded, result, idx_i, idx_j);
        }

        return result;
    }

    std::vector<std::pair<Orbital, NumericalType>> get_sorted_orbital() const
    {
        return get_sorted_orbital_energy(orbital_energy, norb);
    }

    Determinant<N> get_hartree_fock() const
    {
        Determinant<N> hf;
        auto           sorted_orbital = get_sorted_orbital();

        // Number of occupied alpha orbials and occupied beta orbitals.
        auto n_alpha = (nelec + ms2) / 2;
        auto n_beta  = (nelec - ms2) / 2;

        for (auto i_alpha = 0, i_beta = 0, i = 0; i < nelec; ++i)
        {
            auto orb = sorted_orbital[i].first;
            // Alpha orbital
            if ((orb % 2 == 0) && i_alpha < n_alpha)
            {
                hf.set_orbital(orb);
                ++i_alpha;
            }
            // Beta orbital
            else if ((orb % 2 == 1) && i_beta < n_beta)
            {
                hf.set_orbital(orb);
                ++i_beta;
            }
        }

        return hf;
    }
};

template <int N = 1>
class HamiltonianHubbardK : public Hamiltonian<N>
{
    using Hamiltonian<N>::Hamiltonian;
    using OrbitalList = std::vector<Orbital>;
    using typename Hamiltonian<N>::Column;
    using Hamiltonian<N>::norb;
    using Hamiltonian<N>::nelec;
    using Hamiltonian<N>::ms2;

public:
    NumericalType              hubbard_t; // nearest neighbor hopping strength
    NumericalType              hubbard_u; // onsite interaction strength
    NumericalType              coulomb_k; // hubbard_u / (norb / 2)
    int                        dim;       // dimension of the lattice = lattice.size()
    std::vector<Orbital>       lattice;   // orthogonal lattice
    Orbital                    n_alpha;
    Orbital                    n_beta;
    std::vector<NumericalType> orbital_energy;

    int verbose = 2;

    Option opt = {{"hubbard_t", 1.0}, {"hubbard_u", 0.0}, {"lattice", {4}}, {"nelec", 4}, {"ms2", 0}, {"verbose", 2}};

    HamiltonianHubbardK() { init_HamiltonianHubbardK(); }

    HamiltonianHubbardK(Option &option)
    {
        opt.merge_patch(option);
        init_HamiltonianHubbardK();
    }

    HamiltonianHubbardK(NumericalType t, NumericalType u, std::vector<Orbital> &lattice_,
                        int nelec_, int ms2_)
    {
        init_HamiltonianHubbardK(t, u, lattice_, nelec_, ms2_);
    }

    void init_HamiltonianHubbardK()
    {
        NumericalType        hubbard_t_ = opt["hubbard_t"];
        NumericalType        hubbard_u_ = opt["hubbard_u"];
        std::vector<Orbital> lattice_   = opt["lattice"];
        int                  nelec_     = opt["nelec"];
        int                  ms2_       = opt["ms2"];
        verbose                         = opt["verbose"];

        init_HamiltonianHubbardK(hubbard_t_, hubbard_u_, lattice_, nelec_, ms2_);
    }

    void init_HamiltonianHubbardK(NumericalType t, NumericalType u,
                                  std::vector<Orbital> &lattice_, int nelec_, int ms2_)
    {
        nelec     = nelec_;
        ms2       = ms2_;
        hubbard_t = t;
        hubbard_u = u;
        lattice   = lattice_;
        // number of alpha and beta electrons
        n_alpha = (nelec + ms2) / 2;
        n_beta  = (nelec - ms2) / 2;
        // size of lattice
        dim = lattice.size();

        // number of spin-orbitals
        norb = 2;
        for (auto i : lattice)
        {
            if (i < 1)
            {
                throw std::invalid_argument("Length " + std::to_string(i) +
                                            " in the lattice shoule be at least 2.");
            }
            norb *= i;
        }

        validate_option();

        // coulomb interaction = U / #(sites)
        coulomb_k = hubbard_u / (norb / 2);
        // compute orbital energy
        orbital_energy.clear();
        orbital_energy.resize(norb, 0);
        for (auto i = 0; i < norb; ++i)
        {
            auto orb_i = i / 2; // index of orbital
            for (auto d = 0; d < dim; ++d)
            {
                auto wave_number_i = orb_i % lattice[d];
                orbital_energy[i] += cos(2.0 * M_PI * wave_number_i / lattice[d]);
                orb_i /= lattice[d];
            }
            orbital_energy[i] *= (-2 * hubbard_t);
            // std::cout << "Orbital " << i << " energy: " << orbital_energy[i] << std::endl;
        }
    }

    void validate_option()
    {
        if (nelec <= 0)
        {
            throw std::invalid_argument("Number of electrons nelec (" +
                                        std::to_string(norb) + ") shoule be positive.");
        }
        if (norb < nelec)
        {
            throw std::invalid_argument(
                "Number of spin-orbitals norb (" + std::to_string(norb) +
                ") shoule be no less than the number of electrons nelec (" +
                std::to_string(nelec) + "),");
        }
        if ((nelec + ms2) % 2 != 0)
        {
            throw std::invalid_argument("Required nelec (" + std::to_string(nelec) +
                                        ") and ms2 (" + std::to_string(ms2) +
                                        ") are not compatible.");
        }
        if (dim == 0)
        {
            throw std::invalid_argument("Empty lattice.");
        }
        return;
    }

    static int read_norb(Option &option)
    {
        std::vector<Orbital> lattice = option["lattice"];
        int                  norb    = 2;
        for (auto i : lattice)
        {
            if (i <= 0)
            {
                throw std::invalid_argument("Non-positive length " + std::to_string(i) +
                                            " in the lattice.");
            }
            norb *= i;
        }
        return norb;
    }

    void print_header() const
    {
        if (verbose < 2)
            return;

        this->Hamiltonian<N>::print_header();
        Option opt_       = {{"hubbard_k", ""}};
        opt_["hubbard_k"] = opt;
        std::cout << opt_.dump(4) << std::endl
                  << std::endl;
        return;
    }

    void get_diagonal(DeterminantDecoded<N> &det, Column &result) const
    {
        auto value = get_diagonal(det);
        result.push_back({det, value});
    }

    NumericalType get_diagonal(DeterminantDecoded<N> &det) const
    {
        auto occupied_orbitals = det.get_occupied_orbitals();
        return get_diagonal(occupied_orbitals);
    }

    NumericalType get_diagonal(Determinant<N> &det) const
    {
        DeterminantDecoded<N> det_decoded(det);
        return get_diagonal(det_decoded);
    }

    NumericalType get_diagonal(OrbitalList &occupied_orbitals) const
    {
        NumericalType value = 0;

        for (auto idx_i = 0; idx_i < occupied_orbitals.size(); ++idx_i)
        {
            auto i = occupied_orbitals[idx_i];
            // One body integral
            value += orbital_energy[i];
        }
        // Two body integral
        value += coulomb_k * n_alpha * n_beta;
        value -= this->diagonal_shift;
        return value;
    }

    NumericalType get_entry(Determinant<N> &det1, Determinant<N> &det2) const
    {
        NumericalType  value    = 0.0;
        Determinant<N> xor_det  = det1 ^ det2;
        int            popcount = xor_det.popcount();
        if (popcount == 0)
            return get_diagonal(det1);
        else if (popcount == 4)
        { // double excitation from i, j to a, b
            OrbitalList ij_orbitals = (xor_det & det1).get_occupied_orbitals();
            OrbitalList ab_orbitals = (xor_det & det2).get_occupied_orbitals();
            Orbital     i           = ij_orbitals[0];
            Orbital     j           = ij_orbitals[1];
            Orbital     a           = ab_orbitals[0];
            Orbital     b           = ab_orbitals[1];
            // The spin of i and j are different
            if (i % 2 == j % 2)
                return value;
            // The spin of i and a are the same
            if (i % 2 != a % 2)
                std::swap(a, b);

            // momentum conservation
            auto orb_i = i / 2;
            auto orb_j = j / 2;
            auto orb_a = a / 2;
            auto orb_b = determine_b(orb_i, orb_j, orb_a);
            if (b != 2 * orb_b + (j % 2))
                return value;

            value = coulomb_k;
            Determinant<N> new_det(det1);
            // Calculate sign
            new_det.clear_orbital(i);
            auto sign1 = new_det.parity(i, a);
            new_det.set_orbital(a);
            new_det.clear_orbital(j);
            auto sign2 = new_det.parity(j, b);
            new_det.set_orbital(b);
            auto sign = sign1 ^ sign2;
            if (sign)
                value = -value;
        }
        return value;
    }

    // helper function. determine b from momentum conservation
    Orbital determine_b(Orbital orb_i, Orbital orb_j, Orbital orb_a) const
    {
        Orbital orb_b      = 0;
        int     multiplier = 1;
        // compute wave number
        for (auto d = 0; d < dim; ++d)
        {
            auto wave_number_i = orb_i % lattice[d];
            auto wave_number_j = orb_j % lattice[d];
            auto wave_number_a = orb_a % lattice[d];
            // momentum conservation
            auto wave_number_b =
                (wave_number_i + wave_number_j - wave_number_a) % lattice[d];
            wave_number_b =
                (wave_number_b + lattice[d]) %
                lattice[d]; // make sure its nonnegative. strange c++ behavior...
            orb_b += wave_number_b * multiplier;
            // prepare for next dimension
            multiplier *= lattice[d];
            orb_i /= lattice[d];
            orb_j /= lattice[d];
            orb_a /= lattice[d];
        }
        return orb_b;
    }

    // Get all double excitations from occupied orbital index (idx_i, idx_j).
    void get_double_excitation(DeterminantDecoded<N> &det, Column &result, int idx_i,
                               int idx_j) const
    {
        auto i = det.occupied_orbitals[idx_i];
        auto j = det.occupied_orbitals[idx_j];

        // The spin of i and j are different
        if (i % 2 == j % 2)
            return;

        // Loop over all possible orbitals a
        // Note b is determined by momentum conservation
        for (auto a = i % 2; a < norb; a += 2)
        {
            // transform spin-orbital to orbital
            auto orb_i = i / 2;
            auto orb_j = j / 2;
            auto orb_a = a / 2;
            if (!det.is_occupied(a))
            {
                auto orb_b = determine_b(orb_i, orb_j, orb_a);
                // transform orbital to spin-orbital
                auto b = 2 * orb_b + (j % 2);
                // double excitation (i, j) -> (a, b) exists if (a, b) are not occupied
                if (!det.is_occupied(b))
                {
                    NumericalType  value = coulomb_k;
                    Determinant<N> new_det(det);
                    // Calculate sign
                    new_det.clear_orbital(i);
                    auto sign1 = new_det.parity(i, a);
                    new_det.set_orbital(a);
                    new_det.clear_orbital(j);
                    auto sign2 = new_det.parity(j, b);
                    new_det.set_orbital(b);
                    auto sign = sign1 ^ sign2;
                    if (sign)
                        value = -value;

                    // std::cout << i << ' ' << j << ' ' << a << ' ' << b << ' ' << value
                    // << std::endl;
                    result.push_back({new_det, value});
                }
            }
        }
        return;
    }

    Column get_double_excitation(DeterminantDecoded<N> &det) const
    {
        Column result;
        get_double_excitation(det, result);
        return result;
    }

    void get_double_excitation(DeterminantDecoded<N> &det, Column &result) const
    {
        // Double excitation from orbitals (i, j).
        for (auto idx_i = 0; idx_i < det.occupied_orbitals.size(); ++idx_i)
        {
            // auto i = det.occupied_orbitals[idx_i];
            for (auto idx_j = idx_i + 1; idx_j < det.occupied_orbitals.size(); ++idx_j)
            {
                // auto j = det.occupied_orbitals[idx_j];
                get_double_excitation(det, result, idx_i, idx_j);
            }
        }

        return;
    }

    Column get_column(Determinant<N> &det) const
    {
        DeterminantDecoded<N> det_decoded(det);
        Column                result;
        get_diagonal(det_decoded, result);
        // get_single_excitation(det, result);
        get_double_excitation(det_decoded, result);
        return result;
    }

    /* This function only returns the determinants, without the values.*/
    Column get_column_det(Determinant<N> &det) const
    {
        DeterminantDecoded<N> det_decoded(det);
        Column                result;
        result.push_back({det, 0.0});
        // Double excitation from orbitals (i, j).
        for (auto idx_i = 0; idx_i < det_decoded.occupied_orbitals.size(); ++idx_i)
        {
            auto i = det_decoded.occupied_orbitals[idx_i];
            for (auto idx_j = idx_i + 1; idx_j < det_decoded.occupied_orbitals.size(); ++idx_j)
            {
                auto j = det_decoded.occupied_orbitals[idx_j];

                // The spin of i and j are different
                if (i % 2 == j % 2)
                    continue;

                // Loop over all possible orbitals a
                // Note b is determined by momentum conservation
                for (auto a = i % 2; a < norb; a += 2)
                {
                    // transform spin-orbital to orbital
                    auto    orb_i      = i / 2;
                    auto    orb_j      = j / 2;
                    auto    orb_a      = a / 2;
                    Orbital orb_b      = 0;
                    int     multiplier = 1;
                    if (!det.is_occupied(a))
                    {
                        auto orb_b = determine_b(orb_i, orb_j, orb_a);
                        // transform orbital to spin-orbital
                        auto b = 2 * orb_b + (j % 2);
                        // double excitation (i, j) -> (a, b) exists if (a, b) are not occupied
                        if (!det.is_occupied(b))
                        {
                            Determinant<N> new_det(det);
                            // Calculate sign
                            new_det.clear_orbital(i);
                            new_det.set_orbital(a);
                            new_det.clear_orbital(j);
                            new_det.set_orbital(b);
                            result.push_back({new_det, 0.0});
                        }
                    }
                }
            }
        }

        return result;
    }

    /* This function computes the value of given determinants in one column.*/
    void get_column_val(Determinant<N> &det, typename Column::iterator begin, typename Column::iterator end) const
    {
        assert(std::distance(begin, end) >= 0);
        for (auto it = begin; it != end; it++)
            it->second = get_entry(det, it->first);
    }

    Column get_column_diagonal_single(DeterminantDecoded<N> &det) const
    {
        Column result;
        get_diagonal(det, result);
        // get_single_excitation(det, result);
        return result;
    }

    Column get_column_serial_part(Determinant<N> &det) const
    {
        DeterminantDecoded<N> det_decoded(det);
        Column                result;
        get_diagonal(det_decoded, result);

        return result;
    }

    Column get_column_parallel_part(Determinant<N> &det, int tid, int nthreads) const
    {
        DeterminantDecoded<N> det_decoded(det);
        Column                result;

        int n_excitations = nelec * (nelec - 1) / 2;
        // Compute the range the thread needs to compute
        // 0 <= tid < nthreads
        int idx_start = tid * n_excitations / nthreads;
        int idx_end   = (tid + 1) * n_excitations / nthreads;

        for (auto idx = idx_start; idx < idx_end; ++idx)
        {
            int idx_j = 0.5 + sqrt(0.25 + 2 * idx);
            int idx_i = idx - (idx_j - 1) * idx_j / 2;
            get_double_excitation(det_decoded, result, idx_i, idx_j);
        }

        return result;
    }

    std::vector<std::pair<Orbital, NumericalType>> get_sorted_orbital() const
    {
        std::vector<std::pair<Orbital, NumericalType>> result;
        for (auto i = 0; i < norb; ++i)
        {
            result.push_back({i, orbital_energy[i]});
        }
        const auto comp = [](auto &x, auto &y)
        { return x.second < y.second; };
        std::sort(result.begin(), result.end(), comp);
        return result;
    }

    Determinant<N> get_hartree_fock() const
    {
        Determinant<N> hf;
        auto           sorted_orbital = get_sorted_orbital();

        for (auto i_alpha = 0, i_beta = 0, i = 0; i < nelec; ++i)
        {
            auto orb = sorted_orbital[i].first;
            // Alpha orbital
            if ((orb % 2 == 0) && i_alpha < n_alpha)
            {
                hf.set_orbital(orb);
                ++i_alpha;
            }
            // Beta orbital
            else if ((orb % 2 == 1) && i_beta < n_beta)
            {
                hf.set_orbital(orb);
                ++i_beta;
            }
        }

        return hf;
    }
};

template <int N>
int Hamiltonian<N>::read_norb(Option option)
{
    std::string type = option["type"];
    if (type == "molecule")
    {
        return HamiltonianMolecule<N>::read_norb(option["molecule"]);
    }
    if (type == "hubbard_k")
    {
        return HamiltonianHubbardK<N>::read_norb(option["hubbard_k"]);
    }
    throw std::invalid_argument(
        "hamiltonian.type " + type +
        R"( is invalid in the input file. Currently "molecule" and "hubbard_k" are supported.)");
    return 0;
}

template <int N>
std::unique_ptr<Hamiltonian<N>>
Hamiltonian<N>::init(Option option)
{
    std::string type = option["type"];
    if (type == "molecule")
    {
        return std::make_unique<HamiltonianMolecule<N>>(option["molecule"]);
    }
    if (type == "hubbard_k")
    {
        return std::make_unique<HamiltonianHubbardK<N>>(option["hubbard_k"]);
    }
    throw std::invalid_argument(
        "hamiltonian.type " + type +
        R"( is invalid in the input file. Currently "molecule" and "hubbard_k" are supported.)");
    return nullptr;
}

#endif
