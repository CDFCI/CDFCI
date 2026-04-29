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

#ifndef CDFCI_RDM_H
#define CDFCI_RDM_H 1

#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>

#include "wavefunction.h"

// This class defines RDM
template <typename W>
class RDM
{
    using det_type = typename W::key_type;

private:
    /* Options */
    Option      opt = {{"spin_seperate", false},
                       {"compute_1RDM", false},
                       {"compute_2RDM", false},
                       {"method_2RDM", "auto"},
                       {"dump_path", ""},
                       {"verbose", 1}};
    int         norb_;
    int         norb2_; // norb^2
    int         norb3_; // norb^3
    int         nelec_;
    bool        spinsep_; // Spin seperate
    bool        flag1RDM_;
    bool        flag2RDM_;
    std::string method2RDM_;
    std::string dump_path_;
    int         verbose_;
    // data are viewed as a matrix
    NumericalType              data0RDM_;
    std::vector<NumericalType> data1RDM_;           // orbital 1rdm
    std::vector<NumericalType> data2RDM_;           // orbital 2rdm
    std::vector<NumericalType> alpha_data1RDM_;     // orbital 1rdm
    std::vector<NumericalType> alpha_data2RDM_;     // orbital 2rdm
    std::vector<NumericalType> beta_data1RDM_;      // orbital 1rdm
    std::vector<NumericalType> beta_data2RDM_;      // orbital 2rdm
    std::vector<NumericalType> alphabeta_data2RDM_; // alpha/beta spin, aabb

public:
    RDM(Option &option, int norb, int nelec)
    {
        opt.merge_patch(option);

        norb_  = norb;
        norb2_ = norb * norb;
        norb3_ = norb2_ * norb;
        nelec_ = nelec;

        spinsep_    = opt["spin_seperate"];
        flag1RDM_   = opt["compute_1RDM"];
        flag2RDM_   = opt["compute_2RDM"];
        method2RDM_ = opt["method_2RDM"];
        dump_path_  = opt["dump_path"];
        verbose_    = opt["verbose"];
    }

    void compute1RDM(W &vec_xz);

    void compute1RDM();

    NumericalType compute2RDMLin(W &vec_xz, int ntest = 0);

    void symmetrize2RDMLin();

    NumericalType compute2RDMSqr(W &vec_xz, int ntest = 0);

    void symmetrize2RDMSqr();

    void compute2RDM(W &vec_xz)
    {
        int methodflag = 0;

        if (method2RDM_.compare("auto") == 0)
        {
            auto siz_x = vec_xz.size_x();
            if (siz_x > 1e5)
            {
                auto tim1 = compute2RDMLin(vec_xz, 1000);
                auto tim2 = compute2RDMSqr(vec_xz, 50000);
                auto est1 = tim1 * siz_x / 1000;
                auto est2 = tim2 * siz_x * siz_x / 50000 / 50000;
                if (verbose_ > 0)
                {
                    std::cout << "Estimated 2RDM runtime: "
                              << std::setw(10) << std::fixed
                              << std::setprecision(2) << est1
                              << "sec (linear), "
                              << std::setw(10) << std::fixed
                              << std::setprecision(2) << est2
                              << "sec (square). "
                              << std::endl;
                }
                if (est1 < est2)
                    methodflag = 1;
            }
            else if (verbose_ > 0)
            {
                std::cout << "x-vector size is "
                          << std::setw(8) << std::fixed
                          << siz_x << " (<1e5). "
                          << std::endl;
            }

            if (verbose_ > 0)
            {
                if (methodflag == 1)
                    std::cout << "Linear method is used for 2RDM."
                              << std::endl;
                else
                    std::cout << "Square method is used for 2RDM."
                              << std::endl;
            }
        }
        else if (method2RDM_.compare("linear") == 0)
            methodflag = 1;

        NumericalType tim;
        if (methodflag == 1)
        {
            tim = compute2RDMLin(vec_xz);
            symmetrize2RDMLin();
        }
        else
        {
            tim = compute2RDMSqr(vec_xz);
            symmetrize2RDMSqr();
        }

        if (verbose_ > 0)
        {
            std::cout << "Two RDM runtime: "
                      << std::setw(10) << std::fixed
                      << std::setprecision(2) << tim
                      << std::endl;
        }
    }

    void computeRDM(W &vec_xz)
    {
        if (flag1RDM_)
        {
            if (spinsep_)
            {
                alpha_data1RDM_.resize(norb2_, 0);
                beta_data1RDM_.resize(norb2_, 0);
            }
            else
            {
                data1RDM_.resize(norb2_, 0);
            }
        }

        if (flag2RDM_)
        {
            alpha_data2RDM_.resize(norb3_ * norb_, 0);
            beta_data2RDM_.resize(norb3_ * norb_, 0);
            alphabeta_data2RDM_.resize(norb3_ * norb_, 0);
            if (!spinsep_)
                data2RDM_.resize(norb3_ * norb_, 0);
        }

        data0RDM_ = (NumericalType)vec_xz.get_xx();
        if (flag2RDM_)
        {
            compute2RDM(vec_xz);
            if (flag2RDM_)
                compute1RDM();
        }
        else
        {
            if (flag1RDM_)
                compute1RDM(vec_xz);
        }

        if (flag2RDM_ && !spinsep_)
        {
            std::vector<NumericalType>().swap(alpha_data2RDM_);
            std::vector<NumericalType>().swap(beta_data2RDM_);
            std::vector<NumericalType>().swap(alphabeta_data2RDM_);
        }
    }

    void dumpFile()
    {
        if (dump_path_.empty())
            return;

        try
        {
            if (std::filesystem::create_directories(dump_path_))
            {
                std::cout << "RDM dir created successfully: " << dump_path_ << std::endl;
            }
            else
            {
                std::cout << "RDM dir already exist: " << dump_path_ << std::endl;
            }
        }
        catch (const std::filesystem::filesystem_error &e)
        {
            std::cerr << "Error creating directories: " << e.what() << std::endl;
        }

        std::ofstream fzeroRDM(dump_path_ + "/ZeroRDM.out");
        fzeroRDM << std::setprecision(16)
                 << data0RDM_ << std::endl;
        fzeroRDM.close();

        if (flag1RDM_)
        {
            std::ofstream foneRDM(dump_path_ + "/OneRDM.out");
            print1RDM(foneRDM);
            foneRDM.close();
        }

        if (flag2RDM_)
        {
            std::ofstream ftwoRDM(dump_path_ + "/TwoRDM.out");
            print2RDM(ftwoRDM);
            ftwoRDM.close();
        }
    }

    void print1RDM(std::ostream &stream) const;
    void print2RDM(std::ostream &stream) const;

    void output_rdm(NumericalType &zero_rdm, Matrix &one_rdm, Tensor4 &two_rdm);

    ~RDM(){};

    // Get info
    int num_elec() const { return nelec_; }

    int num_orb() const { return norb_; }

    bool spin_seperate() const { return spinsep_; }

    // Get 1RDM
    NumericalType get1RDM(int i, int a) const
    {
        return data1RDM_[i * norb_ + a];
    }

    NumericalType aget1RDM(int i, int a) const
    {
        return alpha_data1RDM_[i * norb_ + a];
    }

    NumericalType bget1RDM(int i, int a) const
    {
        return beta_data1RDM_[i * norb_ + a];
    }

    // Get 2RDM
    NumericalType get2RDM(int i, int a, int j, int b) const
    {
        return data2RDM_[i * norb3_ + a * norb2_ + j * norb_ + b];
    }

    NumericalType aget2RDM(int i, int a, int j, int b) const
    {
        return alpha_data2RDM_[i * norb3_ + a * norb2_ + j * norb_ + b];
    }

    NumericalType bget2RDM(int i, int a, int j, int b) const
    {
        return beta_data2RDM_[i * norb3_ + a * norb2_ + j * norb_ + b];
    }

    NumericalType abget2RDM(int i, int a, int j, int b) const
    {
        return alphabeta_data2RDM_[i * norb3_ + a * norb2_ + j * norb_ + b];
    }

private:
    // Idx 1RDM
    int idx1RDM(int i, int a)
    {
        return i * norb_ + a;
    }

    // Set 1RDM
    void set1RDM(int i, int a, NumericalType value)
    {
        data1RDM_[i * norb_ + a] = value;
    }

    void aset1RDM(int i, int a, NumericalType value)
    {
        alpha_data1RDM_[i * norb_ + a] = value;
    }

    void bset1RDM(int i, int a, NumericalType value)
    {
        beta_data1RDM_[i * norb_ + a] = value;
    }

    // Update 1RDM
    void update1RDM(int i, int a, NumericalType value)
    {
        data1RDM_[i * norb_ + a] += value;
    }

    void aupdate1RDM(int i, int a, NumericalType value)
    {
        alpha_data1RDM_[i * norb_ + a] += value;
    }

    void bupdate1RDM(int i, int a, NumericalType value)
    {
        beta_data1RDM_[i * norb_ + a] += value;
    }

    // Idx 2RDM
    int idx2RDM(int i, int a, int j, int b)
    {
        return i * norb3_ + a * norb2_ + j * norb_ + b;
    }

    // Set 2RDM
    void set2RDM(int i, int a, int j, int b, NumericalType value)
    {
        data2RDM_[i * norb3_ + a * norb2_ + j * norb_ + b] = value;
    }

    void aset2RDM(int i, int a, int j, int b, NumericalType value)
    {
        alpha_data2RDM_[i * norb3_ + a * norb2_ + j * norb_ + b] = value;
    }

    void bset2RDM(int i, int a, int j, int b, NumericalType value)
    {
        beta_data2RDM_[i * norb3_ + a * norb2_ + j * norb_ + b] = value;
    }

    void abset2RDM(int i, int a, int j, int b, NumericalType value)
    {
        alphabeta_data2RDM_[i * norb3_ + a * norb2_ + j * norb_ + b] = value;
    }

    // Update 2RDM
    void update2RDM(int i, int a, int j, int b, NumericalType value)
    {
        data2RDM_[i * norb3_ + a * norb2_ + j * norb_ + b] += value;
    }

    void aupdate2RDM(int i, int a, int j, int b, NumericalType value)
    {
        alpha_data2RDM_[i * norb3_ + a * norb2_ + j * norb_ + b] += value;
    }

    void bupdate2RDM(int i, int a, int j, int b, NumericalType value)
    {
        beta_data2RDM_[i * norb3_ + a * norb2_ + j * norb_ + b] += value;
    }

    void abupdate2RDM(int i, int a, int j, int b, NumericalType value)
    {
        alphabeta_data2RDM_[i * norb3_ + a * norb2_ + j * norb_ + b] += value;
    }
};

/************************************
  Compute 1RDM
*************************************/
// Compute 1RDM from wavefunction
template <typename W>
void RDM<W>::compute1RDM(W &vec_xz)
{
    if (spinsep_)
    {
        alpha_data1RDM_.resize(norb_ * norb_, 0);
        beta_data1RDM_.resize(norb_ * norb_, 0);
    }
    else
        data1RDM_.resize(norb_ * norb_, 0);

    auto time_start = std::chrono::high_resolution_clock::now();
    for (auto iter : vec_xz)
    {
        auto itervalue = iter.second[0];
        if (itervalue == 0)
            continue;
        auto det            = iter.first;
        auto occ_orbitals   = det.get_occupied_orbitals();
        auto unocc_orbitals = det.get_unoccupied_orbitals(2 * norb_);
        // Single annihilation from orbital a
        for (auto a : occ_orbitals)
        {
            int aa = a / 2;
            // Single excitation from orbital i
            for (auto i : unocc_orbitals)
            {
                if ((i % 2) ^ (a % 2))
                    continue;

                int ii = i / 2;

                det_type det_ia(det);
                det_ia.clear_orbital(a);
                auto sign = det_ia.parity(i, a);
                det_ia.set_orbital(i);
                auto value = itervalue * vec_xz.get_x(det_ia);
                if (sign)
                    value = -value;

                if (spinsep_)
                {
                    if (a % 2)
                        bupdate1RDM(ii, aa, value);
                    else
                        aupdate1RDM(ii, aa, value);
                }
                else
                    update1RDM(ii, aa, value);
                det_ia.clear_orbital(i);
            }

            // Single excitation from orbital a
            if (spinsep_)
            {
                if (a % 2)
                    bupdate1RDM(aa, aa, itervalue * itervalue);
                else
                    aupdate1RDM(aa, aa, itervalue * itervalue);
            }
            else
                update1RDM(aa, aa, itervalue * itervalue);
        }
    }
    auto                                 time_end   = std::chrono::high_resolution_clock::now();
    std::chrono::duration<NumericalType> time_elpse = time_end - time_start;
    if (verbose_ > 0)
    {
        std::cout << "One RDM runtime: " << std::setw(10) << std::fixed
                  << std::setprecision(2) << time_elpse.count()
                  << std::endl;
    }
}

// Compute 1RDM from 2RDM
template <typename W>
void RDM<W>::compute1RDM()
{
    if (spinsep_)
    {
        alpha_data1RDM_.resize(norb_ * norb_, 0);
        beta_data1RDM_.resize(norb_ * norb_, 0);
    }
    else
        data1RDM_.resize(norb_ * norb_, 0);

    for (auto a = 0; a < norb_; ++a)
        for (auto i = 0; i < norb_; ++i)
            for (auto k = 0; k < norb_; ++k)
            {
                if (spinsep_)
                {
                    aupdate1RDM(i, a, aget2RDM(i, a, k, k));
                    aupdate1RDM(i, a, abget2RDM(i, a, k, k));
                    bupdate1RDM(i, a, bget2RDM(i, a, k, k));
                    bupdate1RDM(i, a, abget2RDM(i, a, k, k));
                }
                else
                    update1RDM(i, a, get2RDM(i, a, k, k));
            }
    for (auto k = 0; k < norb_ * norb_; ++k)
        if (spinsep_)
        {
            alpha_data1RDM_[k] /= (nelec_ - 1);
            beta_data1RDM_[k] /= (nelec_ - 1);
        }
        else
            data1RDM_[k] /= (nelec_ - 1);
    // Very fast, no need to report runtime
}

/************************************
  Compute 2RDM
*************************************/
template <typename W>
NumericalType RDM<W>::compute2RDMLin(W &vec_xz, int ntest)
{
    using det_type = typename W::key_type;

    auto time_start = std::chrono::high_resolution_clock::now();

    std::vector<NumericalType> alpha_data2RDM;
    std::vector<NumericalType> beta_data2RDM;
    std::vector<NumericalType> alphabeta_data2RDM;
    alpha_data2RDM.resize(norb3_ * norb_, 0);
    beta_data2RDM.resize(norb3_ * norb_, 0);
    alphabeta_data2RDM.resize(norb3_ * norb_, 0);

    int cnt_x = 0;

    for (auto iter : vec_xz)
    {
        auto itervalue = iter.second[0];
        if (itervalue == 0)
            continue;
        cnt_x++;
        if (ntest > 0 && cnt_x > ntest)
        {
            std::chrono::duration<NumericalType> time_elpse =
                std::chrono::high_resolution_clock::now() - time_start;
            return time_elpse.count();
        }
        auto det = iter.first;
        auto occ_alpha_orbitals =
            det.get_occupied_alpha_orbitals();
        auto occ_beta_orbitals =
            det.get_occupied_beta_orbitals();

        // aaaa
        for (auto aiter = occ_alpha_orbitals.begin();
             aiter != occ_alpha_orbitals.end(); ++aiter)
            for (auto biter = aiter + 1;
                 biter != occ_alpha_orbitals.end(); ++biter)
            {
                int      a  = *aiter;
                int      b  = *biter;
                int      aa = a / 2;
                int      bb = b / 2;
                det_type det_ba(det);
                det_ba.clear_orbital(a);
                auto sign1 = det_ba.parity(a, b);
                det_ba.clear_orbital(b);
                auto unocc_alpha_orbitals =
                    det_ba.get_unoccupied_alpha_orbitals(2 * norb_);
                for (auto iiter = unocc_alpha_orbitals.begin();
                     iiter != unocc_alpha_orbitals.end(); ++iiter)
                    for (auto jiter = iiter + 1;
                         jiter != unocc_alpha_orbitals.end(); ++jiter)
                    {
                        int      i  = *iiter;
                        int      j  = *jiter;
                        int      ii = i / 2;
                        int      jj = j / 2;
                        det_type det_iajb(det_ba);
                        det_iajb.set_orbital(j);
                        auto sign2 = det_iajb.parity(i, j);
                        det_iajb.set_orbital(i);
                        auto value = itervalue * vec_xz.get_x(det_iajb);
                        if (sign1 ^ sign2)
                            value = -value;
                        alpha_data2RDM[idx2RDM(ii, aa, jj, bb)] += value;
                    }
            }

        // bbbb
        for (auto aiter = occ_beta_orbitals.begin();
             aiter != occ_beta_orbitals.end(); ++aiter)
            for (auto biter = aiter + 1;
                 biter != occ_beta_orbitals.end(); ++biter)
            {
                int      a  = *aiter;
                int      b  = *biter;
                int      aa = a / 2;
                int      bb = b / 2;
                det_type det_ba(det);
                det_ba.clear_orbital(a);
                auto sign1 = det_ba.parity(a, b);
                det_ba.clear_orbital(b);
                auto unocc_beta_orbitals =
                    det_ba.get_unoccupied_beta_orbitals(2 * norb_);
                for (auto iiter = unocc_beta_orbitals.begin();
                     iiter != unocc_beta_orbitals.end(); ++iiter)
                    for (auto jiter = iiter + 1;
                         jiter != unocc_beta_orbitals.end(); ++jiter)
                    {
                        int      i  = *iiter;
                        int      j  = *jiter;
                        int      ii = i / 2;
                        int      jj = j / 2;
                        det_type det_iajb(det_ba);
                        det_iajb.set_orbital(j);
                        auto sign2 = det_iajb.parity(i, j);
                        det_iajb.set_orbital(i);
                        auto value = itervalue * vec_xz.get_x(det_iajb);
                        if (sign1 ^ sign2)
                            value = -value;
                        beta_data2RDM[idx2RDM(ii, aa, jj, bb)] += value;
                    }
            }

        // aabb
        for (auto a : occ_alpha_orbitals)
            for (auto b : occ_beta_orbitals)
            {
                det_type det_ba(det);
                det_ba.clear_orbital(a);
                auto sign1 = det_ba.parity(a, b);
                det_ba.clear_orbital(b);
                int  aa = a / 2;
                int  bb = b / 2;
                auto unocc_alpha_orbitals =
                    det_ba.get_unoccupied_alpha_orbitals(2 * norb_);
                auto unocc_beta_orbitals =
                    det_ba.get_unoccupied_beta_orbitals(2 * norb_);
                for (auto i : unocc_alpha_orbitals)
                    for (auto j : unocc_beta_orbitals)
                    {
                        int      ii = i / 2;
                        int      jj = j / 2;
                        det_type det_iajb(det_ba);
                        det_iajb.set_orbital(j);
                        auto sign2 = det_iajb.parity(i, j);
                        det_iajb.set_orbital(i);
                        auto value = itervalue * vec_xz.get_x(det_iajb);
                        if (sign1 ^ sign2)
                            value = -value;
                        alphabeta_data2RDM[idx2RDM(ii, aa, jj, bb)] += value;
                    }
            }
    }
    for (auto it = 0; it < alpha_data2RDM.size(); ++it)
        alpha_data2RDM_[it] += alpha_data2RDM[it];
    for (auto it = 0; it < beta_data2RDM.size(); ++it)
        beta_data2RDM_[it] += beta_data2RDM[it];
    for (auto it = 0; it < alphabeta_data2RDM.size(); ++it)
        alphabeta_data2RDM_[it] += alphabeta_data2RDM[it];

    std::chrono::duration<NumericalType> time_elpse =
        std::chrono::high_resolution_clock::now() - time_start;
    return time_elpse.count();
}

template <typename W>
void RDM<W>::symmetrize2RDMLin()
{
    for (auto ii = 0; ii < norb_; ++ii)
        for (auto jj = ii + 1; jj < norb_; ++jj)
            for (auto aa = 0; aa < norb_; ++aa)
                for (auto bb = aa + 1; bb < norb_; ++bb)
                {
                    alpha_data2RDM_[idx2RDM(jj, aa, ii, bb)] =
                        -alpha_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    alpha_data2RDM_[idx2RDM(ii, bb, jj, aa)] =
                        -alpha_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    alpha_data2RDM_[idx2RDM(jj, bb, ii, aa)] =
                        alpha_data2RDM_[idx2RDM(ii, aa, jj, bb)];

                    beta_data2RDM_[idx2RDM(jj, aa, ii, bb)] =
                        -beta_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    beta_data2RDM_[idx2RDM(ii, bb, jj, aa)] =
                        -beta_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    beta_data2RDM_[idx2RDM(jj, bb, ii, aa)] =
                        beta_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                }

    if (spinsep_)
        return;

    for (auto ii = 0; ii < norb_; ++ii)
        for (auto aa = 0; aa < norb_; ++aa)
            for (auto jj = 0; jj < norb_; ++jj)
                for (auto bb = 0; bb < norb_; ++bb)
                {
                    data2RDM_[idx2RDM(ii, aa, jj, bb)] +=
                        alpha_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    data2RDM_[idx2RDM(ii, aa, jj, bb)] +=
                        beta_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    data2RDM_[idx2RDM(ii, aa, jj, bb)] +=
                        alphabeta_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    data2RDM_[idx2RDM(ii, aa, jj, bb)] +=
                        alphabeta_data2RDM_[idx2RDM(jj, bb, ii, aa)];
                }
}

template <typename W>
NumericalType RDM<W>::compute2RDMSqr(W &vec_xz, int ntest)
{
    using wff_type = typename W::wff_type;
    wff_type sub_xz;

    auto time_start = std::chrono::high_resolution_clock::now();

    for (auto iter : vec_xz)
    {
        auto itervalue = iter.second[0];
        if (itervalue == 0)
            continue;
        sub_xz.push_back(iter.first, iter.second);
    }

    if (ntest > 0)
        time_start = std::chrono::high_resolution_clock::now();

#ifndef CDFCI_RDM_SERIAL
#pragma omp parallel
    //shared(sub_xz, spinsep_, norb_, norb2_, norb3_, \
        //alpha_data2RDM_, beta_data2RDM_, alphabeta_data2RDM_, data2RDM_)
    {
        int num_threads = omp_get_num_threads();
        int tid         = omp_get_thread_num();
#else
    int num_threads = 1;
    int tid         = 0;
#endif

        std::vector<NumericalType> alpha_data2RDM;     // orbital 2rdm
        std::vector<NumericalType> beta_data2RDM;      // orbital 2rdm
        std::vector<NumericalType> alphabeta_data2RDM; // alpha/beta spin, aabb
        alpha_data2RDM.resize(norb3_ * norb_, 0);
        beta_data2RDM.resize(norb3_ * norb_, 0);
        alphabeta_data2RDM.resize(norb3_ * norb_, 0);

        auto it_start = sub_xz.begin();
        auto it_end   = sub_xz.begin();
        if (ntest == 0)
        {
            it_start = it_start + (int)(sqrt((NumericalType)tid / num_threads) * sub_xz.size() + 0.5);
            it_end   = it_end + (int)(sqrt((tid + 1.0) / num_threads) * sub_xz.size() + 0.5);
        }
        else
            it_end = it_end + ntest;

        for (auto lit = it_start; lit != it_end; ++lit)
        {
            auto ldet               = (*lit).first;
            auto lval               = (*lit).second[0];
            auto occ_alpha_orbitals = ldet.get_occupied_alpha_orbitals();
            auto occ_beta_orbitals  = ldet.get_occupied_beta_orbitals();
            // aaaa
            for (auto aiter = occ_alpha_orbitals.begin();
                 aiter != occ_alpha_orbitals.end(); ++aiter)
                for (auto biter = aiter + 1;
                     biter != occ_alpha_orbitals.end(); ++biter)
                {
                    int  aa    = *aiter / 2;
                    int  bb    = *biter / 2;
                    auto value = lval * lval;
                    alpha_data2RDM[idx2RDM(aa, aa, bb, bb)] += value;
                }

            // bbbb
            for (auto aiter = occ_beta_orbitals.begin();
                 aiter != occ_beta_orbitals.end(); ++aiter)
                for (auto biter = aiter + 1;
                     biter != occ_beta_orbitals.end(); ++biter)
                {
                    int  aa    = *aiter / 2;
                    int  bb    = *biter / 2;
                    auto value = lval * lval;
                    beta_data2RDM[idx2RDM(aa, aa, bb, bb)] += value;
                }

            // aabb
            for (auto a : occ_alpha_orbitals)
                for (auto b : occ_beta_orbitals)
                {
                    int  aa    = a / 2;
                    int  bb    = b / 2;
                    auto value = lval * lval;
                    alphabeta_data2RDM[idx2RDM(aa, aa, bb, bb)] += value;
                }

            for (auto rit = sub_xz.begin(); rit != lit; ++rit)
            {
                auto rdet     = (*rit).first;
                auto diff_det = ldet ^ rdet;
                int  cnt      = diff_det.popcount();

                // No 3 or more body terms in 2RDM
                if (cnt > 4)
                    continue;

                auto rval = (*rit).second[0];
                auto diff_occ_alpha_orbitals =
                    diff_det.get_occupied_alpha_orbitals();
                auto diff_occ_beta_orbitals =
                    diff_det.get_occupied_beta_orbitals();

                if (cnt == 2)
                {
                    int i, a;
                    if (diff_occ_alpha_orbitals.size() > 0)
                    {
                        i = diff_occ_alpha_orbitals[0];
                        a = diff_occ_alpha_orbitals[1];
                    }
                    else
                    {
                        i = diff_occ_beta_orbitals[0];
                        a = diff_occ_beta_orbitals[1];
                    }
                    int      ii    = i / 2;
                    int      aa    = a / 2;
                    auto     value = lval * rval;
                    det_type det(ldet);
                    if (ldet.is_occupied(a))
                    {
                        det.clear_orbital(a);
                        auto sign = det.parity(a, i);
                        if (sign)
                            value = -value;
                    }
                    else
                    {
                        det.clear_orbital(i);
                        auto sign = det.parity(i, a);
                        if (sign)
                            value = -value;
                    }
                    if (diff_occ_alpha_orbitals.size() > 0)
                    {
                        // aaaa
                        for (auto biter = occ_alpha_orbitals.begin();
                             biter != occ_alpha_orbitals.end(); ++biter)
                        {
                            if ((*biter == i) || (*biter == a))
                                continue;
                            int bb = *biter / 2;
                            alpha_data2RDM[idx2RDM(ii, aa, bb, bb)] += value;
                        }

                        // aabb
                        for (auto biter = occ_beta_orbitals.begin();
                             biter != occ_beta_orbitals.end(); ++biter)
                        {
                            int bb = *biter / 2;
                            alphabeta_data2RDM[idx2RDM(ii, aa, bb, bb)] += value;
                        }
                    }
                    else
                    {
                        // bbbb
                        for (auto biter = occ_beta_orbitals.begin();
                             biter != occ_beta_orbitals.end(); ++biter)
                        {
                            if ((*biter == i) || (*biter == a))
                                continue;
                            int bb = *biter / 2;
                            beta_data2RDM[idx2RDM(ii, aa, bb, bb)] += value;
                        }

                        // aabb
                        for (auto biter = occ_alpha_orbitals.begin();
                             biter != occ_alpha_orbitals.end(); ++biter)
                        {
                            int bb = *biter / 2;
                            alphabeta_data2RDM[idx2RDM(bb, bb, ii, aa)] += value;
                        }
                    }
                }

                if (cnt == 4)
                {
                    // aaaa or bbbb
                    if (diff_occ_alpha_orbitals.size() != 2)
                    {
                        std::vector<int> ldiff_occ_orbitals;
                        std::vector<int> ldiff_unocc_orbitals;
                        if (diff_occ_alpha_orbitals.size() == 4)
                        {
                            for (auto pos : diff_occ_alpha_orbitals)
                                if (ldet.is_occupied(pos))
                                    ldiff_occ_orbitals.push_back(pos);
                                else
                                    ldiff_unocc_orbitals.push_back(pos);
                        }
                        else
                        {
                            for (auto pos : diff_occ_beta_orbitals)
                                if (ldet.is_occupied(pos))
                                    ldiff_occ_orbitals.push_back(pos);
                                else
                                    ldiff_unocc_orbitals.push_back(pos);
                        }
                        int      i, j, a, b;
                        auto     value = lval * rval;
                        det_type det(ldet);
                        if (ldiff_occ_orbitals[0] < ldiff_unocc_orbitals[0])
                        {
                            i = ldiff_occ_orbitals[0];
                            j = ldiff_occ_orbitals[1];
                            a = ldiff_unocc_orbitals[0];
                            b = ldiff_unocc_orbitals[1];
                            det.clear_orbital(i);
                            auto sign1 = det.parity(i, a);
                            det.set_orbital(a);
                            det.clear_orbital(j);
                            auto sign2 = det.parity(j, b);
                            if (sign1 ^ sign2)
                                value = -value;
                        }
                        else
                        {
                            i = ldiff_unocc_orbitals[0];
                            j = ldiff_unocc_orbitals[1];
                            a = ldiff_occ_orbitals[0];
                            b = ldiff_occ_orbitals[1];
                            det.clear_orbital(a);
                            auto sign1 = det.parity(a, i);
                            det.set_orbital(i);
                            det.clear_orbital(b);
                            auto sign2 = det.parity(b, j);
                            if (sign1 ^ sign2)
                                value = -value;
                        }
                        int ii = i / 2;
                        int jj = j / 2;
                        int aa = a / 2;
                        int bb = b / 2;
                        if (diff_occ_alpha_orbitals.size() == 4)
                            alpha_data2RDM[idx2RDM(ii, aa, jj, bb)] += value;
                        else
                            beta_data2RDM[idx2RDM(ii, aa, jj, bb)] += value;
                    }
                    else // if (diff_occ_alpha_orbitals.size() == 2)
                    {
                        int      i = diff_occ_alpha_orbitals[0];
                        int      a = diff_occ_alpha_orbitals[1];
                        int      j = diff_occ_beta_orbitals[0];
                        int      b = diff_occ_beta_orbitals[1];
                        det_type det(ldet);
                        auto     value = lval * rval;
                        if (ldet.is_occupied(a))
                        {
                            if (ldet.is_occupied(j))
                            {
                                j = diff_occ_beta_orbitals[1];
                                b = diff_occ_beta_orbitals[0];
                            }
                            det.clear_orbital(a);
                            auto sign1 = det.parity(a, i);
                            det.set_orbital(i);
                            det.clear_orbital(b);
                            auto sign2 = det.parity(b, j);
                            if (sign1 ^ sign2)
                                value = -value;
                        }
                        else
                        {
                            if (ldet.is_occupied(b))
                            {
                                j = diff_occ_beta_orbitals[1];
                                b = diff_occ_beta_orbitals[0];
                            }
                            det.clear_orbital(i);
                            auto sign1 = det.parity(i, a);
                            det.set_orbital(a);
                            det.clear_orbital(j);
                            auto sign2 = det.parity(j, b);
                            if (sign1 ^ sign2)
                                value = -value;
                        }
                        int ii = i / 2;
                        int jj = j / 2;
                        int aa = a / 2;
                        int bb = b / 2;
                        alphabeta_data2RDM[idx2RDM(ii, aa, jj, bb)] += value;
                    }
                }
            }
        }

#ifndef CDFCI_RDM_SERIAL
#pragma omp critical
        {
#endif
            if (ntest == 0)
            {
                for (auto it = 0; it < alpha_data2RDM.size(); ++it)
                    alpha_data2RDM_[it] += alpha_data2RDM[it];
                for (auto it = 0; it < beta_data2RDM.size(); ++it)
                    beta_data2RDM_[it] += beta_data2RDM[it];
                for (auto it = 0; it < alphabeta_data2RDM.size(); ++it)
                    alphabeta_data2RDM_[it] += alphabeta_data2RDM[it];
            }
#ifndef CDFCI_RDM_SERIAL
        } // end omp critical
#endif

#ifndef CDFCI_RDM_SERIAL
    } // end omp parallel
#endif
    std::chrono::duration<NumericalType> time_elpse =
        std::chrono::high_resolution_clock::now() - time_start;
    return time_elpse.count();
}

template <typename W>
void RDM<W>::symmetrize2RDMSqr()
{
    for (auto aa = 0; aa < norb_; ++aa)
        for (auto bb = aa + 1; bb < norb_; ++bb)
        {
            alpha_data2RDM_[idx2RDM(bb, bb, aa, aa)] =
                alpha_data2RDM_[idx2RDM(aa, aa, bb, bb)];
            alpha_data2RDM_[idx2RDM(bb, aa, aa, bb)] =
                -alpha_data2RDM_[idx2RDM(aa, aa, bb, bb)];
            alpha_data2RDM_[idx2RDM(aa, bb, bb, aa)] =
                -alpha_data2RDM_[idx2RDM(aa, aa, bb, bb)];

            beta_data2RDM_[idx2RDM(bb, bb, aa, aa)] =
                beta_data2RDM_[idx2RDM(aa, aa, bb, bb)];
            beta_data2RDM_[idx2RDM(bb, aa, aa, bb)] =
                -beta_data2RDM_[idx2RDM(aa, aa, bb, bb)];
            beta_data2RDM_[idx2RDM(aa, bb, bb, aa)] =
                -beta_data2RDM_[idx2RDM(aa, aa, bb, bb)];
        }

    for (auto ii = 0; ii < norb_; ++ii)
        for (auto aa = ii + 1; aa < norb_; ++aa)
            for (auto bb = 0; bb < norb_; ++bb)
            {
                if ((bb == ii) || (bb == aa))
                    continue;

                alpha_data2RDM_[idx2RDM(bb, aa, ii, bb)] =
                    -alpha_data2RDM_[idx2RDM(ii, aa, bb, bb)];
                alpha_data2RDM_[idx2RDM(ii, bb, bb, aa)] =
                    -alpha_data2RDM_[idx2RDM(ii, aa, bb, bb)];
                alpha_data2RDM_[idx2RDM(bb, bb, ii, aa)] =
                    alpha_data2RDM_[idx2RDM(ii, aa, bb, bb)];
                alpha_data2RDM_[idx2RDM(aa, ii, bb, bb)] =
                    alpha_data2RDM_[idx2RDM(ii, aa, bb, bb)];
                alpha_data2RDM_[idx2RDM(bb, ii, aa, bb)] =
                    -alpha_data2RDM_[idx2RDM(ii, aa, bb, bb)];
                alpha_data2RDM_[idx2RDM(aa, bb, bb, ii)] =
                    -alpha_data2RDM_[idx2RDM(ii, aa, bb, bb)];
                alpha_data2RDM_[idx2RDM(bb, bb, aa, ii)] =
                    alpha_data2RDM_[idx2RDM(ii, aa, bb, bb)];

                beta_data2RDM_[idx2RDM(bb, aa, ii, bb)] =
                    -beta_data2RDM_[idx2RDM(ii, aa, bb, bb)];
                beta_data2RDM_[idx2RDM(ii, bb, bb, aa)] =
                    -beta_data2RDM_[idx2RDM(ii, aa, bb, bb)];
                beta_data2RDM_[idx2RDM(bb, bb, ii, aa)] =
                    beta_data2RDM_[idx2RDM(ii, aa, bb, bb)];
                beta_data2RDM_[idx2RDM(aa, ii, bb, bb)] =
                    beta_data2RDM_[idx2RDM(ii, aa, bb, bb)];
                beta_data2RDM_[idx2RDM(bb, ii, aa, bb)] =
                    -beta_data2RDM_[idx2RDM(ii, aa, bb, bb)];
                beta_data2RDM_[idx2RDM(aa, bb, bb, ii)] =
                    -beta_data2RDM_[idx2RDM(ii, aa, bb, bb)];
                beta_data2RDM_[idx2RDM(bb, bb, aa, ii)] =
                    beta_data2RDM_[idx2RDM(ii, aa, bb, bb)];
            }

    for (auto ii = 0; ii < norb_; ++ii)
        for (auto aa = ii + 1; aa < norb_; ++aa)
            for (auto jj = ii + 1; jj < norb_; ++jj)
                for (auto bb = aa + 1; bb < norb_; ++bb)
                {
                    if (bb == jj)
                        continue;
                    alpha_data2RDM_[idx2RDM(jj, aa, ii, bb)] =
                        -alpha_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    alpha_data2RDM_[idx2RDM(ii, bb, jj, aa)] =
                        -alpha_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    alpha_data2RDM_[idx2RDM(jj, bb, ii, aa)] =
                        alpha_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    alpha_data2RDM_[idx2RDM(aa, ii, bb, jj)] =
                        alpha_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    alpha_data2RDM_[idx2RDM(aa, jj, bb, ii)] =
                        -alpha_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    alpha_data2RDM_[idx2RDM(bb, ii, aa, jj)] =
                        -alpha_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    alpha_data2RDM_[idx2RDM(bb, jj, aa, ii)] =
                        alpha_data2RDM_[idx2RDM(ii, aa, jj, bb)];

                    beta_data2RDM_[idx2RDM(jj, aa, ii, bb)] =
                        -beta_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    beta_data2RDM_[idx2RDM(ii, bb, jj, aa)] =
                        -beta_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    beta_data2RDM_[idx2RDM(jj, bb, ii, aa)] =
                        beta_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    beta_data2RDM_[idx2RDM(aa, ii, bb, jj)] =
                        beta_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    beta_data2RDM_[idx2RDM(aa, jj, bb, ii)] =
                        -beta_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    beta_data2RDM_[idx2RDM(bb, ii, aa, jj)] =
                        -beta_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    beta_data2RDM_[idx2RDM(bb, jj, aa, ii)] =
                        beta_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                }

    for (auto ii = 0; ii < norb_; ++ii)
        for (auto aa = ii + 1; aa < norb_; ++aa)
            for (auto bb = 0; bb < norb_; ++bb)
                alphabeta_data2RDM_[idx2RDM(bb, bb, aa, ii)] =
                    alphabeta_data2RDM_[idx2RDM(bb, bb, ii, aa)];

    for (auto ii = 0; ii < norb_; ++ii)
        for (auto aa = ii + 1; aa < norb_; ++aa)
            for (auto jj = 0; jj < norb_; ++jj)
                for (auto bb = 0; bb < norb_; ++bb)
                    alphabeta_data2RDM_[idx2RDM(aa, ii, bb, jj)] =
                        alphabeta_data2RDM_[idx2RDM(ii, aa, jj, bb)];

    if (spinsep_)
        return;

    for (auto ii = 0; ii < norb_; ++ii)
        for (auto aa = 0; aa < norb_; ++aa)
            for (auto jj = 0; jj < norb_; ++jj)
                for (auto bb = 0; bb < norb_; ++bb)
                {
                    data2RDM_[idx2RDM(ii, aa, jj, bb)] +=
                        alpha_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    data2RDM_[idx2RDM(ii, aa, jj, bb)] +=
                        beta_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    data2RDM_[idx2RDM(ii, aa, jj, bb)] +=
                        alphabeta_data2RDM_[idx2RDM(ii, aa, jj, bb)];
                    data2RDM_[idx2RDM(ii, aa, jj, bb)] +=
                        alphabeta_data2RDM_[idx2RDM(jj, bb, ii, aa)];
                }
}

/************************************
  Print RDMs
*************************************/
template <typename W>
void RDM<W>::print1RDM(std::ostream &stream) const
{
    if (spinsep_)
    {
        for (auto val : alpha_data1RDM_)
            stream << std::setprecision(16) << val << " ";
        stream << "\n";
        for (auto val : beta_data1RDM_)
            stream << std::setprecision(16) << val << " ";
        stream << "\n";
    }
    else
    {
        for (auto val : data1RDM_)
            stream << std::setprecision(16) << val << " ";
        stream << "\n";
    }
    return;
}

template <typename W>
void RDM<W>::print2RDM(std::ostream &stream) const
{
    if (spinsep_)
    {
        for (auto val : alpha_data2RDM_)
            stream << std::setprecision(16) << val << " ";
        stream << "\n";

        for (auto val : beta_data2RDM_)
            stream << std::setprecision(16) << val << " ";
        stream << "\n";

        for (auto val : alphabeta_data2RDM_)
            stream << std::setprecision(16) << val << " ";
        stream << "\n";
    }
    else
    {
        for (auto val : data2RDM_)
            stream << std::setprecision(16) << val << " ";
        stream << "\n";
    }
    return;
}

template <typename W>
void RDM<W>::output_rdm(NumericalType &zero_rdm, Matrix &one_rdm, Tensor4 &two_rdm)
{
    zero_rdm = data0RDM_;
    one_rdm  = Eigen::Map<Matrix>(data1RDM_.data(), norb_, norb_);
    two_rdm  = Eigen::TensorMap<Tensor4>(data2RDM_.data(), norb_, norb_, norb_, norb_);
}

#endif
