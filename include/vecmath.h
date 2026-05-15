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

#ifndef CDFCI_VECMATH_H
#define CDFCI_VECMATH_H

#include <cmath>
#include "config.h"

/*
 * vecmath.h
 *
 * Namespace: vecmath
 * Description: A lightweight linear algebra utility library for scalar and fixed-size array operations.
 *
 * Naming Conventions:
 *   - All functions use `snake_case`.
 *   - Operations that modify the input in-place are named with `_assign`, e.g. `add_assign`.
 *   - Operations restricted to array halves are suffixed with `_first_half` or `_second_half`.
 *   - Read-only utility functions follow mathematical naming, e.g. `max_norm`, `multiply`.
 *   - Scalar and array overloads are supported consistently.
 *
 * Example Usage:
 *   using namespace vecmath;
 *   std::array<double, 4> x = {1, 2, 0, 0};
 *   if (is_zero_on_second_half(x)) ...
 */

namespace vecmath {

    // Extend std::is_floating_point to include QUAD_PRECISION
    template <typename T>
    struct is_floating_point_ex : std::is_floating_point<T> {};

    template <>
    struct is_floating_point_ex<QUAD_PRECISION> : std::true_type {};

    template <typename T>
    constexpr bool is_floating_point_ex_v = is_floating_point_ex<T>::value;

    // Scalar addition
    template <typename T>
    std::enable_if_t<is_floating_point_ex_v<T>, void>
    add_assign(T& lhs, const T& rhs) {
        lhs += rhs;
    }

    // Array addition
    template <typename T, size_t N>
    void add_assign(std::array<T, N>& lhs, const std::array<T, N>& rhs) {
        for (size_t i = 0; i < N; ++i)
            lhs[i] += rhs[i];
    }

    // Array partial addition (first half)
    template <typename T, size_t N, typename RhsType>
    void add_assign_first_half(std::array<T, N>& lhs, const RhsType& rhs) {
        if constexpr (std::is_same_v<RhsType, T>) {
            for (size_t i = 0; i < N / 2; ++i)
                lhs[i] += rhs;
        } else {
            for (size_t i = 0; i < N / 2; ++i)
                lhs[i] += rhs[i];
        }
    }

    // Array partial addition (second half)
    template <typename T, size_t N, typename RhsType>
    void add_assign_second_half(std::array<T, N>& lhs, const RhsType& rhs) {
        if constexpr (std::is_same_v<RhsType, T>) {
            for (size_t i = N / 2; i < N; ++i)
                lhs[i] += rhs;
        } else {
            for (size_t i = N / 2; i < N; ++i)
                lhs[i] += rhs[i - N / 2];
        }
    }

    // Set values to first half
    template <typename T, size_t N, typename RhsType>
    void assign_first_half(std::array<T, N>& lhs, const RhsType& rhs) {
        if constexpr (std::is_same_v<RhsType, T>) {
            for (size_t i = 0; i < N / 2; ++i)
                lhs[i] = rhs;
        } else {
            for (size_t i = 0; i < N / 2; ++i)
                lhs[i] = rhs[i];
        }
    }

    // Set values to second half
    template <typename T, size_t N, typename RhsType>
    void assign_second_half(std::array<T, N>& lhs, const RhsType& rhs) {
        if constexpr (std::is_same_v<RhsType, T>) {
            for (size_t i = N / 2; i < N; ++i)
                lhs[i] = rhs;
        } else {
            for (size_t i = N / 2; i < N; ++i)
                lhs[i] = rhs[i - N / 2];
        }
    }

    // Scalar multiplication
    template <typename T>
    T multiply(const T& lhs, const T& rhs) {
        return lhs * rhs;
    }

    // Array * scalar multiplication
    template <typename T, size_t N>
    std::array<T, N> multiply(const std::array<T, N>& lhs, const T& rhs) {
        std::array<T, N> result{};
        for (size_t i = 0; i < N; ++i)
            result[i] = lhs[i] * rhs;
        return result;
    }

    // Slice first half
    template <typename T, size_t N>
    typename std::enable_if<(N > 2), std::array<T, N / 2>>::type
    slice_first_half(const std::array<T, N>& arr) {
        std::array<T, N / 2> result{};
        for (size_t i = 0; i < N / 2; ++i)
            result[i] = arr[i];
        return result;
    }

    template <typename T>
    T slice_first_half(const std::array<T, 2>& arr) {
        return arr[0];
    }

    // Slice second half
    template <typename T, size_t N>
    typename std::enable_if<(N > 2), std::array<T, N / 2>>::type
    slice_second_half(const std::array<T, N>& arr) {
        std::array<T, N / 2> result{};
        for (size_t i = N / 2; i < N; ++i)
            result[i - N / 2] = arr[i];
        return result;
    }

    template <typename T>
    T slice_second_half(const std::array<T, 2>& arr) {
        return arr[1];
    }

    // Check if first half is zero
    template <typename T, size_t N>
    bool is_zero_on_first_half(const std::array<T, N>& arr) {
        for (size_t i = 0; i < N / 2; ++i)
            if (arr[i] != 0) return false;
        return true;
    }

    template <typename T>
    bool is_zero_on_first_half(const std::array<T, 2>& arr) {
        return (arr[0] == 0);
    }

    // Check if second half is zero
    template <typename T, size_t N>
    bool is_zero_on_second_half(const std::array<T, N>& arr) {
        for (size_t i = N / 2; i < N; ++i)
            if (arr[i] != 0) return false;
        return true;
    }

    template <typename T>
    bool is_zero_on_second_half(const std::array<T, 2>& arr) {
        return (arr[1] == 0);
    }

    // index
    template <typename T>
    T index(const T& data, std::size_t /*unused*/) {
        return data;
    }

    template <typename T, std::size_t N>
    T index(const std::array<T, N>& data, std::size_t idx) {
        return data[idx];
    }

    // Two-norm of array or scalar
    template <typename T>
    std::enable_if_t<is_floating_point_ex_v<T>, NumericalType>
    two_norm(const T& val) {
        return fabs(val);
    }

    template <typename T, size_t N>
    NumericalType two_norm(const std::array<T, N>& arr) {
        NumericalType result = 0.0;
        for (size_t i = 0; i < N; ++i)
            result += arr[i] * arr[i];
        return sqrt(result);
    }

    // Max-norm of array or scalar
    template <typename T>
    std::enable_if_t<is_floating_point_ex_v<T>, NumericalType>
    max_norm(const T& val) {
        return fabs(val);
    }

    template <typename T, size_t N>
    NumericalType max_norm(const std::array<T, N>& arr) {
        NumericalType result = 0.0;
        for (size_t i = 0; i < N; ++i)
            result = std::max(result, std::fabs(arr[i]));
        return result;
    }

    template <typename T, size_t N>
    NumericalType inner_product(const std::array<T, N>& x,
                                const std::array<T, N>& y) {
        NumericalType result = 0.0;
        for (size_t i = 0; i < N; ++i)
            result += x[i] * y[i];
        return result;
    }
}

#endif