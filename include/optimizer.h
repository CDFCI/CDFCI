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

#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <iostream>
#include <random>
#include "option.h"
#include "config.h"

namespace contracter
{
    Matrix rotate_one_int(const Matrix &one_int, const Matrix &U)
    {
        /*
            ┌───────┐
            │one_int│
            └┬─────┬┘
            ┌┴┐   ┌┴┐
            │U│   │U│
            └┬┘   └┬┘
        */
        return U.transpose() * one_int * U;
    }

    Tensor4 rotate_two_int(const Tensor4 &two_int, const Matrix &U)
    {
        /*
            ┌──────────┐
            │  two_int │
            └┬──┬──┬──┬┘
            ┌┴┐┌┴┐┌┴┐┌┴┐
            │U││U││U││U│
            └┬┘└┬┘└┬┘└┬┘
        */
        Eigen::Tensor<NumericalType, 2> U_tensor(U.rows(), U.cols());
        for (int i = 0; i < U.rows(); ++i)
        {
            for (int j = 0; j < U.cols(); ++j)
            {
                U_tensor(i, j) = U(i, j);
            }
        }
        Eigen::array<Eigen::IndexPair<int>, 1> contract_dims = {Eigen::IndexPair<int>(0, 0)};
        Tensor4                                contract_1    = two_int.contract(U_tensor, contract_dims);
        Tensor4                                contract_2    = contract_1.contract(U_tensor, contract_dims);
        Tensor4                                contract_3    = contract_2.contract(U_tensor, contract_dims);
        Tensor4                                contract_4    = contract_3.contract(U_tensor, contract_dims);
        return contract_4;
    }

    NumericalType compute_one_body_energy(const Matrix &one_int, const Matrix &U, const Matrix &one_rdm)
    {
        /*
            ┌───────┐
            │one_int│
            └┬─────┬┘
            ┌┴┐   ┌┴┐
            │U│   │U│
            └┬┘   └┬┘
            ┌┴─────┴┐
            │one_rdm│
            └───────┘
        */
        Matrix contract_1 = rotate_one_int(one_int, U);
        Matrix contract_2 = contract_1 * one_rdm;
        return contract_2.trace();
    }

    NumericalType compute_two_body_energy(const Tensor4 &two_int, const Matrix &U, const Tensor4 &two_rdm)
    {
        /*
            ┌──────────┐
            │  two_int │
            └┬──┬──┬──┬┘
            ┌┴┐┌┴┐┌┴┐┌┴┐
            │U││U││U││U│
            └┬┘└┬┘└┬┘└┬┘
            ┌┴──┴──┴──┴┐
            │  two_rdm │
            └──────────┘
        */
        Tensor4                                contract_1    = rotate_two_int(two_int, U);
        Eigen::array<Eigen::IndexPair<int>, 4> contract_dims = {
            Eigen::IndexPair<int>(0, 0), Eigen::IndexPair<int>(1, 1),
            Eigen::IndexPair<int>(2, 2), Eigen::IndexPair<int>(3, 3)};
        Eigen::Tensor<NumericalType, 0> contract_2 = contract_1.contract(two_rdm, contract_dims);
        return contract_2(0);
    }

    Matrix compute_one_body_gradient(const Matrix &one_int, const Matrix &U, const Matrix &one_rdm)
    {
        /*
            ┌───────┐
            │one_int│
            └┬─────┬┘
            ┌┴┐
            │U│
            └┬┘
            ┌┴─────┴┐
            │one_rdm│
            └───────┘
        */
        return one_int * U * one_rdm;
    }

    Matrix compute_two_body_gradient(const Tensor4 &two_int, const Matrix &U, const Tensor4 &two_rdm)
    {
        /*
            ┌──────────┐
            │  two_int │
            └┬──┬──┬──┬┘
            ┌┴┐┌┴┐┌┴┐
            │U││U││U│
            └┬┘└┬┘└┬┘
            ┌┴──┴──┴──┴┐
            │  two_rdm │
            └──────────┘
        */
        Eigen::Tensor<NumericalType, 2> U_tensor(U.rows(), U.cols());
        for (int i = 0; i < U.rows(); ++i)
        {
            for (int j = 0; j < U.cols(); ++j)
            {
                U_tensor(i, j) = U(i, j);
            }
        }
        Eigen::array<Eigen::IndexPair<int>, 1> contract_dims   = {Eigen::IndexPair<int>(0, 0)};
        Tensor4                                contract_1      = two_int.contract(U_tensor, contract_dims);
        Tensor4                                contract_2      = contract_1.contract(U_tensor, contract_dims);
        Tensor4                                contract_3      = contract_2.contract(U_tensor, contract_dims);
        Eigen::array<Eigen::IndexPair<int>, 3> contract_dims_3 = {Eigen::IndexPair<int>(1, 0), Eigen::IndexPair<int>(2, 1), Eigen::IndexPair<int>(3, 2)};
        Eigen::Tensor<NumericalType, 2>        contract_4      = contract_3.contract(two_rdm, contract_dims_3);
        Matrix                                 V(U.rows(), U.cols());
        for (int i = 0; i < U.rows(); ++i)
        {
            for (int j = 0; j < U.cols(); ++j)
            {
                V(i, j) = contract_4(i, j);
            }
        }
        return V;
    }
}

class Optimizer
{
public:
    Option        opt;
    NumericalType zero_int; // core_energy
    Matrix        one_int;  // morb × morb, index: (i, a)
    Tensor4       two_int;  // morb⁴      , index: (i, a, j, b)
    NumericalType zero_rdm; // <Φ|Φ>
    Matrix        one_rdm;  // norb × norb, index: (i, a)
    Tensor4       two_rdm;  // norb⁴      , index: (i, a, j, b)

    NumericalType final_value;

public:
    Optimizer() {}

    Optimizer(Option opt, NumericalType zero_int, Matrix one_int, Tensor4 two_int)
    {
        this->opt      = opt;
        this->zero_int = zero_int;
        this->one_int  = one_int;
        this->two_int  = two_int;
    }

    Optimizer(Option opt, NumericalType zero_int, Matrix one_int, Tensor4 two_int, NumericalType zero_rdm, Matrix one_rdm, Tensor4 two_rdm)
    {
        this->opt      = opt;
        this->zero_int = zero_int;
        this->one_int  = one_int;
        this->two_int  = two_int;
        this->zero_rdm = zero_rdm;
        this->one_rdm  = one_rdm;
        this->two_rdm  = two_rdm;
    }

    ~Optimizer() {}

    void update_rdm(NumericalType zero_rdm, Matrix one_rdm, Tensor4 two_rdm)
    {
        this->zero_rdm = zero_rdm;
        this->one_rdm  = one_rdm;
        this->two_rdm  = two_rdm;
    }

    /*
                                        ┌───────┐         ┌──────────┐
                                        │one_int│         │  two_int │
                                        └┬─────┬┘         └┬──┬──┬──┬┘
                                        ┌┴┐   ┌┴┐         ┌┴┐┌┴┐┌┴┐┌┴┐
       <Φ|H|Φ> = core_energy * <Φ|Φ>  + │U│   │U│ + 1/2 * │U││U││U││U│
                                        └┬┘   └┬┘         └┬┘└┬┘└┬┘└┬┘
                                        ┌┴─────┴┐         ┌┴──┴──┴──┴┐
                                        │one_rdm│         │  two_rdm │
                                        └───────┘         └──────────┘
     */

    inline NumericalType f(const Matrix &U)
    {
        return zero_int +
               +(contracter::compute_one_body_energy(one_int, U, one_rdm) + contracter::compute_two_body_energy(two_int, U, two_rdm) / 2) / zero_rdm;
    }

    inline Matrix gradient(const Matrix &U)
    {
        return 2 * contracter::compute_one_body_gradient(one_int, U, one_rdm) + 2 * contracter::compute_two_body_gradient(two_int, U, two_rdm);
    }

    inline Matrix gradient_step(const Matrix &U, const NumericalType &tau, const Matrix &grad)
    {
        return U - tau * grad;
    }

    inline NumericalType matrix_inner_product(const Matrix &A, const Matrix &B)
    {
        return (A.array() * B.array()).sum();
    }

    Matrix project_to_manifold(const Matrix &A)
    {
        Eigen::JacobiSVD<Matrix> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);

        Matrix          U = svd.matrixU();
        Matrix          V = svd.matrixV();
        Eigen::VectorXd S = svd.singularValues();
        return U * V.transpose();
    }

    Matrix add_Gaussian_noise(const Matrix &A, const NumericalType &stddev)
    {
        std::random_device         rd{};
        std::mt19937               gen{rd()};
        std::normal_distribution<> d{0, stddev};

        Matrix noise(A.rows(), A.cols());
        for (int i = 0; i < noise.rows(); ++i)
        {
            for (int j = 0; j < noise.cols(); ++j)
            {
                noise(i, j) = d(gen);
            }
        }

        return A + noise;
    }

    NumericalType compute_bb_stepsize(const Matrix &U1, const Matrix &U0, const Matrix &G1, const Matrix &G0, bool is_odd)
    {
        Matrix        S  = U1 - U0;
        Matrix        V  = G1 - G0;
        NumericalType SV = matrix_inner_product(S, V);
        NumericalType tau;
        if (is_odd)
            tau = fabs(SV) / matrix_inner_product(V, V);
        else
            tau = matrix_inner_product(S, S) / fabs(SV);
        return tau;
    }

    Matrix optimize_U(Matrix U0)
    {
        int           max_iterations = opt["max_iterations"];
        int           report_freq    = opt["report_freq"];
        NumericalType tol            = opt["tolerance"];
        NumericalType tau            = opt["initial_stepsize"];

        Matrix        U           = U0;
        Matrix        grad0       = gradient(U0);
        Matrix        grad        = grad0;
        NumericalType E0          = f(U0);
        NumericalType E           = E0;
        NumericalType accum_error = 0.0;

        printf("%6d: %20.16f, %16e\n", 0, E, accum_error);
        for (int it = 0; it < max_iterations / report_freq; ++it)
        {
            for (int jt = 0; jt < report_freq; ++jt)
            {
                int iter = it * report_freq + jt + 1;
                U        = gradient_step(U, tau, grad);
                U        = project_to_manifold(U);
                grad     = gradient(U);
                tau      = compute_bb_stepsize(U, U0, grad, grad0, iter % 2 == 1);

                U0    = U;
                grad0 = grad;
            }
            E0 = E;
            E  = f(U);
            if (it == 0)
                accum_error = fabs(E - E0);
            else
                accum_error = 0.2 * fabs(E - E0) + 0.8 * accum_error;

            printf("%6d: %20.16f, %16e\n", (it + 1) * report_freq, E, accum_error);

            if (accum_error < tol)
            {
                final_value = E;
                return U;
            }
        }
        final_value = E;
        return U;
    }
};

#endif