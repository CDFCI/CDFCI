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

#ifndef CDFCI_CONFIG_H
#define CDFCI_CONFIG_H

#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <array>
#include <istream>
#include <ostream>
#include <type_traits>

/* version */
#define CDFCI_VERSION_MAJOR 1 // for incompatible API changes
#define CDFCI_VERSION_MINOR 0 // for adding functionality in a backwards-compatible manner
#define CDFCI_VERSION_PATCH 0 // for backwards-compatible bug fixes

/* macros */

// Count Trailing Zero
// GCC builtin function
#ifdef __x86_64__
#define CTZ __builtin_ctzll
#define PARITY __builtin_parityll
#define POPCOUNT __builtin_popcountll
#else
#define CTZ __builtin_ctz
#define PARITY __builtin_parity
#define POPCOUNT __builtin_popcount
#endif

// OpenMP
#ifdef _OPENMP
#include <omp.h>
#else
#define CDFCI_SOLVER_SERIAL
#define CDFCI_RDM_SERIAL
#endif

/* Self defined timer */
class Timer {
public:
    void measure(const std::string& label, const std::function<void()>& func) {
        auto start = std::chrono::high_resolution_clock::now();
        func();
        auto end = std::chrono::high_resolution_clock::now();
        times[label] += std::chrono::duration<double>(end - start).count();
    }

    void print_results() const {
        for (const auto& [label, time] : times) {
            std::cout << std::setw(20) << std::fixed << std::setprecision(2) << label << ": " << time << " seconds" << std::endl;
        }
    }

private:
    std::unordered_map<std::string, double> times;
};

// float128 support
// temporary solution for clang: use long double instead
#if defined(__clang__)
#define CDFCI_USE_LONG_DOUBLE
#define QUAD_PRECISION long double
#else
#define QUAD_PRECISION __float128
#include <quadmath.h>

inline QUAD_PRECISION fabs(QUAD_PRECISION x) { return fabsq(x); }
inline QUAD_PRECISION sqrt(QUAD_PRECISION x) { return sqrtq(x); }
inline QUAD_PRECISION cbrt(QUAD_PRECISION x) { return cbrtq(x); }
inline QUAD_PRECISION cos(QUAD_PRECISION x) { return cosq(x); }
inline QUAD_PRECISION atan2(QUAD_PRECISION y, QUAD_PRECISION x) { return atan2q(y, x); }
#endif


/* Global type definition */
using Orbital         = int; // spin-orbitals
using OrbitalPair     = std::pair<Orbital, Orbital>;
using NumericalType   = double;
using ScaleFactorType = QUAD_PRECISION;

/* Eigen Type Aliases */
#include "lib/eigen/Eigen/Dense"
#include "lib/eigen/unsupported/Eigen/CXX11/Tensor"

using Matrix  = Eigen::Matrix<NumericalType, Eigen::Dynamic, Eigen::Dynamic>;
using Tensor4 = Eigen::Tensor<NumericalType, 4>;


// --------- Read/Write Tool ---------
// ------------------ detection utilities (SFINAE) --------------------
template <typename, typename = void>
struct has_member_dump : std::false_type {};
template <typename T>
struct has_member_dump<T, std::void_t<
    decltype(std::declval<const T&>().dump(std::declval<std::ostream&>()))
>> : std::true_type {};

template <typename, typename = void>
struct has_member_load : std::false_type {};
template <typename T>
struct has_member_load<T, std::void_t<
    decltype(std::declval<T&>().load(std::declval<std::istream&>()))
>> : std::true_type {};

template <typename T>
constexpr bool has_member_dump_v = has_member_dump<T>::value;
template <typename T>
constexpr bool has_member_load_v = has_member_load<T>::value;

// is_std_array<T> trait
template <typename T> struct is_std_array : std::false_type {};
template <typename U, std::size_t N>
struct is_std_array<std::array<U,N>> : std::true_type {};
template <typename T>
constexpr bool is_std_array_v = is_std_array<std::decay_t<T>>::value;

// -------------------- unified dump/load --------------------
template <typename T>
inline bool dump(std::ostream& os, const T& obj);

// load must be declared before used recursively
template <typename T>
inline bool load(std::istream& is, T& obj);

// dump: raw -> array -> member -> static_assert
template <typename T>
inline bool dump(std::ostream& os, const T& obj) {
    using U = std::decay_t<T>;

    if constexpr (std::is_trivially_copyable<U>::value && !is_std_array_v<U>) {
        // trivially copyable non-array: write raw bytes
        os.write(reinterpret_cast<const char*>(&obj), static_cast<std::streamsize>(sizeof(U)));
        return static_cast<bool>(os);
    } else if constexpr (is_std_array_v<U>) {
        // std::array path
        using Elem = typename U::value_type;
        if constexpr (std::is_trivially_copyable<Elem>::value) {
            const auto bytes = static_cast<std::streamsize>(sizeof(Elem) * obj.size());
            os.write(reinterpret_cast<const char*>(obj.data()), bytes);
            return static_cast<bool>(os);
        } else {
            for (const auto& e : obj)
                if (!dump(os, e)) return false; // recurse
            return true;
        }
    } else if constexpr (has_member_dump_v<U>) {
        // custom type defines member dump(std::ostream&)
        obj.dump(os);
        return static_cast<bool>(os);
    } else {
        static_assert(std::is_trivially_copyable<U>::value
                      || is_std_array_v<U>
                      || has_member_dump_v<U>,
                      "dump(os, obj): unsupported type. "
                      "Make it trivially copyable, or std::array thereof, "
                      "or provide obj.dump(std::ostream&).");
        return false;
    }
}

// load: raw -> array -> member -> static_assert
template <typename T>
inline bool load(std::istream& is, T& obj) {
    using U = std::decay_t<T>;

    if constexpr (std::is_trivially_copyable<U>::value && !is_std_array_v<U>) {
        is.read(reinterpret_cast<char*>(&obj), static_cast<std::streamsize>(sizeof(U)));
        return static_cast<bool>(is);
    } else if constexpr (is_std_array_v<U>) {
        using Elem = typename U::value_type;
        if constexpr (std::is_trivially_copyable<Elem>::value) {
            const auto bytes = static_cast<std::streamsize>(sizeof(Elem) * obj.size());
            is.read(reinterpret_cast<char*>(obj.data()), bytes);
            return static_cast<bool>(is);
        } else {
            for (auto& e : obj)
                if (!load(is, e)) return false; // recurse
            return true;
        }
    } else if constexpr (has_member_load_v<U>) {
        obj.load(is);
        return static_cast<bool>(is);
    } else {
        static_assert(std::is_trivially_copyable<U>::value
                      || is_std_array_v<U>
                      || has_member_load_v<U>,
                      "load(is, obj): unsupported type. "
                      "Make it trivially copyable, or std::array thereof, "
                      "or provide obj.load(std::istream&).");
        return false;
    }
}

// --------- Result struct ---------
struct Result {
    NumericalType energy = 0.0;
    std::vector<NumericalType> state_energies;
    size_t iterations = 0;
    size_t report_interval = 0;
    std::vector<NumericalType> energy_history;
    std::vector<size_t> x_size_history;
    std::vector<size_t> z_size_history;
    std::vector<double> time_history;
};

#endif
