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

#ifndef CDFCI_COORDINATE_STRATEGY_H
#define CDFCI_COORDINATE_STRATEGY_H 1

#include "config.h"
#include "vecmath.h"
#include <memory>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_set>
#include <string>
#include <stdexcept>
#include <cassert>

/* Coordinate Pick functors*/
/* Mapping from WaveFunctionFragment to WaveFunctionFragment*/
template <typename W>
class CoordinatePick
{
    using wff_type = typename W::wff_type;

public:
    static std::unique_ptr<CoordinatePick<W>>
    init(const std::string strategy);

    virtual wff_type
    operator()(const W &vec_xz, wff_type &det_list, size_t num_coordinates = 1) const = 0;

    virtual ~CoordinatePick() {}
};

template<typename W>
class CoordinatePickGcdGrad : public CoordinatePick<W>
{
    using wff_type   = typename W::wff_type;
    using key_type   = typename wff_type::key_type;
    using value_type = typename wff_type::value_type;
    using data_type  = typename W::data_type;
    using data_matrix_type  = typename W::data_matrix_type;
    using hasher     = typename W::hasher;
    using key_equal  = typename W::key_equal;

public:
    ~CoordinatePickGcdGrad(){};

    /** Straightforward implementation of N choose 1. */
    wff_type
    operator()(const W &vec_xz, wff_type &det_list, size_t num_coordinates = 1) const
    {
        wff_type   result;
        value_type max_det;
        // Warning: assume det_list is non-empty.
        NumericalType abs_grad, max_abs_grad = 0;
        data_matrix_type xx = vec_xz.get_xx_double();
        for (const auto &det : det_list)
        {
            data_type x = vecmath::slice_first_half(det.second);
            data_type z = vecmath::slice_second_half(det.second);
            // data_type g = vec_xz.get_grad(x, z);
            data_type g = {};
            if constexpr (W::NSTATES_val == 1) {
                g = (x * xx + z);
            } else {
                for (int i = 0; i < W::NSTATES_val; ++i) {
                    g[i] = z[i];
                    for (int k = 0; k < W::NSTATES_val; ++k)
                        g[i] += x[k] * xx[k + i * W::NSTATES_val];
                }
            }

            abs_grad = vecmath::two_norm(g);
            if (abs_grad >= max_abs_grad)
            {
                max_abs_grad = abs_grad;
                max_det      = det;
            }
        }
        result.push_back(max_det.first, max_det.second);
        return result;
    }
};

template<typename W, typename = std::enable_if_t<W::NSTATES_val == 1>>
class CoordinatePickBlockGcdGrad : public CoordinatePick<W>
{
    using wff_type   = typename W::wff_type;
    using key_type   = typename wff_type::key_type;
    using value_type = typename wff_type::value_type;
    using hasher     = typename W::hasher;
    using key_equal  = typename W::key_equal;

public:
    ~CoordinatePickBlockGcdGrad(){};

    /** Implementation of N choose M using Min Heap + remove duplicates, single-threaded */
    wff_type
    operator()(const W &vec_xz, wff_type &det_list, size_t num_coordinates = 1) const
    {
        if (num_coordinates > det_list.size())
            num_coordinates = det_list.size();

        // storing the smallest value on the top
        using HeapElement             = std::pair<NumericalType, value_type>;
        const auto compareHeapElement = [](const HeapElement &lhs, const HeapElement &rhs)
        {
            return lhs.first > rhs.first;
        };
        std::priority_queue<HeapElement, std::vector<HeapElement>, decltype(compareHeapElement)> min_heap(compareHeapElement);

        // maintaining a set of existing determinants in the heap
        std::unordered_set<key_type, hasher, key_equal> existing_det;

        NumericalType w = vec_xz.get_xx() * vec_xz.get_scale() * vec_xz.get_scale();

        for (auto &det : det_list)
        {
            auto          x        = det.second[0];
            auto          z        = det.second[1];
            NumericalType abs_grad = fabs(x * w + z);
            // early stop
            if ((min_heap.size() >= num_coordinates) && (min_heap.top().first > abs_grad))
                continue;

            if (existing_det.find(det.first) != existing_det.end())
                continue;

            // fetch the correct z
            z             = vec_xz.get_z(det.first);
            det.second[1] = z;
            abs_grad      = fabs(x * w + z);

            if (min_heap.size() < num_coordinates)
            {
                existing_det.insert(det.first);
                min_heap.emplace(abs_grad, det);
            }
            else if (min_heap.top().first < abs_grad)
            {
                existing_det.erase(min_heap.top().second.first);
                min_heap.pop();
                existing_det.insert(det.first);
                min_heap.emplace(abs_grad, det);
            }
        }

        assert(min_heap.size() == existing_det.size());

        // return result
        wff_type result;
        while (!min_heap.empty())
        {
            result.push_back(min_heap.top().second);
            min_heap.pop();
        }

        return result;
    }

    /** Implementation of N choose M using Min Heap + remove duplicates, multi-threaded
     * with bug, not in use
    wff_type
    operator()(const W &vec_xz, wff_type &det_list, size_t num_coordinates = 1) const
    {
        size_t list_size = det_list.size();
        if (num_coordinates > list_size)
            num_coordinates = list_size;

        using HeapElement = std::pair<NumericalType, value_type>;
        const auto compareHeapElement = [](const HeapElement &lhs, const HeapElement &rhs) {
            return lhs.first > rhs.first;
        };
        NumericalType w = vec_xz.get_xx() * vec_xz.get_scale() * vec_xz.get_scale(); // truncate to low precision

#ifndef CDFCI_SOLVER_SERIAL
        std::vector<HeapElement> global_vector;
        size_t p = omp_get_max_threads();
#else
        size_t p = 1;
#endif

        // assigning tasks
        size_t chunk_size = (list_size + p - 1) / p;
        std::vector<typename wff_type::iterator> chunk_starts(p + 1, det_list.end());
        chunk_starts[0] = det_list.begin();
        for (size_t i = 1; i <= p; ++i) {
            auto it = chunk_starts[i - 1];
            std::advance(it, chunk_size);
            if (it == det_list.end())
                break;
            else
                chunk_starts[i] = it;
        }

#ifndef CDFCI_SOLVER_SERIAL
#pragma omp parallel firstprivate(w)
{
#endif
        // storing the smallest value on the top
        std::priority_queue<HeapElement, std::vector<HeapElement>, decltype(compareHeapElement)> min_heap(compareHeapElement);

        // maintaining a set of existing determinants in the heap
        std::unordered_set<key_type, hasher, key_equal> existing_det;

#ifndef CDFCI_SOLVER_SERIAL
        size_t tid = omp_get_thread_num();
#else
        size_t tid = 0;
#endif
        auto start_it = chunk_starts[tid];
        auto end_it = chunk_starts[tid + 1];

        for (auto it = start_it; it != end_it; ++it)
        {
            auto &det = *it;
            auto x = det.second[0];
            auto z = det.second[1];
            NumericalType abs_grad = fabs(x * w + z);
            // early stop
            if ((min_heap.size() >= num_coordinates) && (min_heap.top().first > abs_grad))
                continue;

            if (existing_det.find(det.first) != existing_det.end())
                continue;

            // fetch the correct z
            z = vec_xz.get_z(det.first);
            det.second[1] = z;
            abs_grad = fabs(x * w + z);

            if (min_heap.size() < num_coordinates) {
                existing_det.insert(det.first);
                min_heap.emplace(abs_grad, det);
            }
            else if (min_heap.top().first < abs_grad) {
                existing_det.erase(min_heap.top().second.first);
                min_heap.pop();
                existing_det.insert(det.first);
                min_heap.emplace(abs_grad, det);
            }
        }

        assert(min_heap.size() == existing_det.size());

#ifndef CDFCI_SOLVER_SERIAL
        #pragma omp critical
        {
            while (!min_heap.empty()) {
                global_vector.push_back(min_heap.top());
                min_heap.pop();
            }
        }
}
        std::sort(global_vector.begin(), global_vector.end(), [](const auto &a, const auto &b) {
            return a.first > b.first;
        });

        // return result (remove duplicates)
        wff_type result;
        key_equal KeyEqual;
        key_type last_seen;
        for (size_t idx = 0; result.size() != num_coordinates && idx < global_vector.size(); ++idx) {
            if (KeyEqual(global_vector[idx].second.first, last_seen))
                continue;
            last_seen = global_vector[idx].second.first;
            result.push_back(global_vector[idx].second);
        }
        return result;
#else
        // return result
        wff_type result;
        while (!min_heap.empty()) {
            result.push_back(min_heap.top().second);
            min_heap.pop();
        }

        return result;
#endif
    } */
};

template <typename W>
std::unique_ptr<CoordinatePick<W>>
CoordinatePick<W>::init(const std::string strategy)
{
    if (strategy == "gcd_grad")
        return std::make_unique<CoordinatePickGcdGrad<W>>();

    if constexpr (W::NSTATES_val == 1)
        if (strategy == "block_gcd_grad")
            return std::make_unique<CoordinatePickBlockGcdGrad<W>>();

    throw std::invalid_argument("coordinate_pick " + strategy +
                                " is invalid in the input file.");
}

/* Coordinate Update functors*/
/* Update WaveFunctionFragment in-place and return ScaleFactor*/
template <typename H, typename W>
class CoordinateUpdate
{
    using wff_type   = typename W::wff_type;
    using value_type = typename wff_type::value_type;

public:
    static std::unique_ptr<CoordinateUpdate<H, W>>
    init(const std::string strategy);

    virtual ScaleFactorType
    operator()(const H &h, const W &vec_xz,
               wff_type &det_picked) const = 0;

    virtual ~CoordinateUpdate(){};
};

/* Use Line Search at each coordinate to determine stepsize. */
template<typename H, typename W>
class CoordinateUpdateLS : public CoordinateUpdate<H,W>
{
    using wff_type   = typename W::wff_type;
    using value_type = typename wff_type::value_type;
    using data_type  = typename W::data_type;
    using data_matrix_type  = typename W::data_matrix_type;

public:
    ~CoordinateUpdateLS(){};

    /**
     * Solve the cubic equation: x^3 + b*x^2 + c*x + d = 0
     * with Newton's Iteration
     *
     * Returns one real root based on the following rule:
     *   1. If there is exactly one real root (and two complex),
     *      return that root.
     *   2. If there are two real roots (one simple, one double),
     *      return the simple root.
     *   3. If there are three distinct real roots, return the one that is
     *      furthest from the middle root (i.e., one of the two side roots).
     */
    double solve_cubic(double b, double c, double d,
                       const double depsilon,
                       const int max_iter) const {
        // Transform to depressed cubic form: y^3 + p*y + q = 0
        // using the substitution x = y - b / 3 to eliminate the
        // quadratic term
        double p = c - b * b / 3.0;
        double q = d + b * (2.0 / 27.0 * b * b - c / 3.0);
        double shift = - b / 3.0; // final x = y + shift

        // Compute discriminant to determine root types
        double p3 = p / 3.0;
        double q2 = q / 2.0;
        double delta = q2 * q2 + p3 * p3 * p3;
        double root = 0;
        const double pi = M_PI;

        if (delta >= 0) {
            // One real root (and possibly a double complex root)
            double sqrt_delta = sqrt(delta);
            root = cbrt(-q2 + sqrt_delta) + cbrt(-q2 - sqrt_delta);
        } else {
            // Three distinct real roots — select one further from the center
            double r = sqrt(-p3);
            double qrtd = sqrt(-delta);

            double theta;
            if (q2 >= 0) {
                // Adjust angle to pick side root instead of the middle one
                theta = (atan2(-qrtd, -q2) - 2.0 * pi) / 3.0;
            } else {
                theta = atan2(qrtd, -q2) / 3.0;
            }

            root = 2.0 * r * cos(theta);
        }

        double x = root + shift;

        // Use Newton iteration to improve accuracy since root formula
        // may be numerically unstable
        // auto f = [](double z){
        //     return z * (z * (z + b) + c) + d;
        // };
        // auto df = [](double z){
        //     return z * (3.0 * z + 2.0 * b) + c;
        // };
        auto newton_one_step = [b, c, d](double z){
            double f = z * (z * (z + b) + c) + d;
            double df = z * (3.0 * z + 2.0 * b) + c;
            return z - f / df;
        };

        int iter = 0;

        double xn = newton_one_step(x);
        while (fabs((xn - x)/x) > depsilon && iter < max_iter)
        {
            x = xn;
            xn = newton_one_step(x);
            ++iter;
        }
        return xn;
    }

    double solve_cubic(double a, double b, double c, double d,
                       const double depsilon,
                       const int max_iter) const
    {
        return solve_cubic(b / a, c / a, d / a, depsilon, max_iter);
    }

    ScaleFactorType
    operator()(const H &h, const W &vec_xz,
               wff_type &det_list) const
    {
        data_matrix_type xx = vec_xz.get_xx_double();
        const double depsilon = 1e-10; // convergence criterion
        const int max_iter = 10; // maximum iterations for Newton's method

        for (auto &det_picked : det_list)
        {
            auto det = det_picked.first;
            auto x   = vecmath::slice_first_half(det_picked.second);
            auto z   = vecmath::slice_second_half(det_picked.second);
            auto g   = vec_xz.get_grad(x, z);
            auto dA  = h.get_diagonal(det);

            data_type result;

            // Solve: x^3 + b*x^2 + c*x + d = 0
            // with the following b, c, d
            if constexpr (W::NSTATES_val == 1) {
                auto b = 3 * x;
                auto c = dA + 2 * x * x + xx;
                auto d = g;

                result = solve_cubic(b, c, d, depsilon, max_iter);
            } else {
                auto x2 = vecmath::inner_product(x, x);
                auto g2 = vecmath::inner_product(g, g);
                auto gTx = vecmath::inner_product(x, g);

                auto gTxxg = 0.0;
                for (int i = 0; i < W::NSTATES_val; i++)
                    for (int j = 0; j < W::NSTATES_val; j++)
                        gTxxg += g[i] * xx[i + j * W::NSTATES_val] * g[j];

                auto a = g2;
                auto b = 3 * gTx;
                auto c = dA + gTx * gTx / g2 + x2 + gTxxg / g2;
                auto d = 1;

                auto tau = solve_cubic(a, b, c, d, depsilon, max_iter);
                result = vecmath::multiply(g, tau);
            }

            vecmath::assign_first_half(det_picked.second, result);
        }

        return 1.0;
    }
};

/* Solve a local eigenvalue problem to determine stepsize. */
template<typename H, typename W, typename = std::enable_if_t<W::NSTATES_val == 1>>
class CoordinateUpdateEig : public CoordinateUpdate<H,W>
{
    using wff_type   = typename W::wff_type;
    using value_type = typename wff_type::value_type;

public:
    ~CoordinateUpdateEig(){};

    ScaleFactorType
    operator()(const H &h, const W &vec_xz,
               wff_type &det_list) const
    {
        size_t          n = det_list.size(); // problem size
        Eigen::VectorXd x_vec(n);            // x[I]
        Eigen::VectorXd z_vec(n);            // z[I]
        Eigen::MatrixXd H_mat(n, n);         // H[I, I]
        size_t          i = 0, j;
        for (auto it = det_list.begin(); it != det_list.end(); ++it, ++i)
        {
            auto det    = it->first;
            x_vec(i)    = it->second[0];
            z_vec(i)    = it->second[1];
            H_mat(i, i) = h.get_diagonal(det);
            j           = i + 1;
            for (auto it_j = it + 1; it_j != det_list.end(); ++it_j, ++j)
            {
                H_mat(j, i) = h.get_entry(det, it_j->first);
            }
        }
        Eigen::SelfAdjointView<Eigen::MatrixXd, Eigen::Lower> H_sym = H_mat.selfadjointView<Eigen::Lower>();

        // All the computations below don't involve scale, which is cancelled out.
        auto xx     = vec_xz.get_xx();     // x'*x
        auto xz     = vec_xz.get_xz();     // x'*z
        auto xTx    = x_vec.squaredNorm(); // x[I]'*x[I]
        auto xTz    = x_vec.dot(z_vec);    // x[I]'*z[I]
        auto Hx_vec = H_sym * x_vec;       // H[I, I]*x[I]
        auto xTHx   = x_vec.dot(Hx_vec);   // x[I]'*H[I, I]*x[I]

        auto yy = xx - xTx;
        auto y1 = sqrt(yy);

        QUAD_PRECISION lambda, gamma, scale;
        if (fabs(y1) > std::numeric_limits<double>::epsilon())
        {
            Eigen::MatrixXd A_mat(n + 1, n + 1);
            A_mat(0, 0)                                                 = (xz - 2 * xTz + xTHx) / yy;
            A_mat.block(1, 0, n, 1)                                     = (z_vec - Hx_vec) / y1;
            A_mat.block(1, 1, n, n)                                     = H_mat;
            Eigen::SelfAdjointView<Eigen::MatrixXd, Eigen::Lower> A_sym = A_mat.selfadjointView<Eigen::Lower>();
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>        solver(A_sym);

            if (solver.info() != Eigen::Success)
            {
                throw std::runtime_error("Eigen decomposition failed");
            }

            NumericalType   eval = solver.eigenvalues()[0];      // the minimum eigenvalue
            Eigen::VectorXd evec = solver.eigenvectors().col(0); // corresponding eigenvector

            lambda = -eval;
            gamma  = evec(0) / y1;
            scale  = gamma * sqrt(lambda);

            // dx = x_new - x_old = a / scale - x;
            i = 0;
            for (auto it = det_list.begin(); it != det_list.end(); ++it, ++i)
                it->second[0] = evec(i + 1) / gamma - x_vec(i);
        }
        else
        {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H_sym);

            if (solver.info() != Eigen::Success)
            {
                throw std::runtime_error("Eigen decomposition failed");
            }

            NumericalType   eval = solver.eigenvalues()[0];      // the minimum eigenvalue
            Eigen::VectorXd evec = solver.eigenvectors().col(0); // corresponding eigenvector

            scale = 1.0;
            // dx = x_new - x_old = a - x;
            i = 0;
            for (auto it = det_list.begin(); it != det_list.end(); ++it, ++i)
                it->second[0] = evec(i) * sqrt(-eval) - x_vec(i);
        }
        return scale;
    }
};

template <typename H, typename W>
std::unique_ptr<CoordinateUpdate<H, W>>
CoordinateUpdate<H, W>::init(const std::string strategy)
{
    if constexpr (W::NSTATES_val == 1) {
        if (strategy == "ls")
        {
            return std::make_unique<CoordinateUpdateLS<H, W>>();
        }
        if (strategy == "eig")
        {
            return std::make_unique<CoordinateUpdateEig<H, W>>();
        }
        throw std::invalid_argument("coordinate_ls " + strategy +
                                    " is invalid in the input file.");
        return nullptr;
    } else {
        if (strategy == "ls")
        {
            return std::make_unique<CoordinateUpdateLS<H, W>>();
        }
        throw std::invalid_argument("coordinate_ls " + strategy +
                                    " is invalid in the input file.");
    }
}

#endif
