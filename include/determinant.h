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

#ifndef CDFCI_DETERMINANT_H
#define CDFCI_DETERMINANT_H 1

#include "config.h"
#include "lib/robin_hood/robin_hood.h" // only use the hash function
#include <climits>
#include <vector>
#include <memory>

template <int N = 1>
class Determinant
{
public:
    using det_type =
        size_t; // Use of size_t is consistent with the macro definition above.

    constexpr static int size     = N;
    constexpr static int det_size = sizeof(det_type) * CHAR_BIT;

    // Masks for parity calculation.
    static std::array<std::array<det_type, size>, det_size * size> masks;

    // Q: I would like to declare them as protected. But then how do I construct
    // DeterminantDecoded?
    // Bit representation of the Slater determinant, using N numbers
    // of det_type integers.
    std::array<det_type, size> repr;

    Determinant()
    {
        // Initialize repr to 0.
        std::fill(repr.begin(), repr.end(), 0);
    }
    Determinant(std::array<det_type, N> input_repr) : repr{input_repr} {}
    Determinant(const Determinant<N> &det) { repr = det.repr; }

    // Need to be called explicitly to construct the masks.
    static void constuct_masks()
    {
        for (auto i = 0; i < masks.size(); ++i)
        {
            masks[i][i / det_size] = (static_cast<det_type>(1) << (i % det_size)) - 1;
            for (auto j = 0; j < i / det_size; ++j)
            {
                masks[i][j] = ~static_cast<det_type>(0); // all bits are 1.
            }
            for (auto j = i / det_size + 1; j < size; ++j)
            {
                masks[i][j] = 0;
            }
        }
        return;
    }

    // FLAG: efficiency needed.
    // Return the occupied orbitals in a vector, from the smallest to the largest one.
    std::vector<Orbital> get_occupied_orbitals() const
    {
        std::vector<Orbital> orbitals;
        for (auto i = 0; i < size; ++i)
        {
            det_type det = repr[i];
            while (det)
            {
                int orb = CTZ(det); // the index of the least significant bit
                orbitals.push_back(orb + det_size * i);
                det &= (det - 1); // unset the least significant bit
            }
        }
        return orbitals;
    }

    std::vector<Orbital> get_occupied_alpha_orbitals() const
    {
        std::vector<Orbital> orbitals;
        for (auto i = 0; i < size; ++i)
        {
            det_type det = repr[i];
            while (det)
            {
                int orb = CTZ(det); // the index of the least significant bit
                if ((orb + det_size * i) % 2 == 0)
                    orbitals.push_back(orb + det_size * i);
                det &= (det - 1); // unset the least significant bit
            }
        }
        return orbitals;
    }

    std::vector<Orbital> get_occupied_beta_orbitals() const
    {
        std::vector<Orbital> orbitals;
        for (auto i = 0; i < size; ++i)
        {
            det_type det = repr[i];
            while (det)
            {
                int orb = CTZ(det); // the index of the least significant bit
                if ((orb + det_size * i) % 2 == 1)
                    orbitals.push_back(orb + det_size * i);
                det &= (det - 1); // unset the least significant bit
            }
        }
        return orbitals;
    }

    std::vector<Orbital> get_unoccupied_orbitals(int norb) const
    {
        std::vector<Orbital> orbitals;
        for (auto i = 0; i < size; ++i)
        {
            // Left here
            det_type det = ~repr[i];
            while (det)
            {
                int orb = CTZ(det); // the index of the least significant bit
                if (orb + det_size * i >= norb)
                    return orbitals;
                orbitals.push_back(orb + det_size * i);
                det &= (det - 1); // unset the least significant bit
            }
        }
        return orbitals;
    }

    std::vector<Orbital> get_unoccupied_alpha_orbitals(int norb) const
    {
        std::vector<Orbital> orbitals;
        for (auto i = 0; i < size; ++i)
        {
            // Left here
            det_type det = ~repr[i];
            while (det)
            {
                int orb = CTZ(det); // the index of the least significant bit
                if (orb + det_size * i >= norb)
                    return orbitals;
                if ((orb + det_size * i) % 2 == 0)
                    orbitals.push_back(orb + det_size * i);
                det &= (det - 1); // unset the least significant bit
            }
        }
        return orbitals;
    }

    std::vector<Orbital> get_unoccupied_beta_orbitals(int norb) const
    {
        std::vector<Orbital> orbitals;
        for (auto i = 0; i < size; ++i)
        {
            // Left here
            det_type det = ~repr[i];
            while (det)
            {
                int orb = CTZ(det); // the index of the least significant bit
                if (orb + det_size * i >= norb)
                    return orbitals;
                if ((orb + det_size * i) % 2 == 1)
                    orbitals.push_back(orb + det_size * i);
                det &= (det - 1); // unset the least significant bit
            }
        }
        return orbitals;
    }

    // FLAG: efficiency needed.
    bool is_occupied(Orbital i) const
    {
        return repr[i / det_size] & (static_cast<det_type>(1) << (i % det_size));
    }

    void set_orbital(Orbital i)
    {
        repr[i / det_size] |= (static_cast<det_type>(1) << (i % det_size));
    }

    void clear_orbital(Orbital i)
    {
        repr[i / det_size] &= ~(static_cast<det_type>(1) << (i % det_size));
    }

    // Return the number of 1's mod 2 in bits min(i, a) <= k < max(i, a).
    int parity(Orbital i, Orbital a) const
    {
        int result = 0;
        for (auto t = 0; t < size; ++t)
        {
            result ^= PARITY(repr[t] & (masks[i][t] ^ masks[a][t]));
        }
        return result;
    }

    int popcount() const
    {
        int result = 0;
        for (auto t = 0; t < size; ++t)
        {
            result += POPCOUNT(repr[t]);
        }
        return result;
    }

    Determinant<N> operator&(Determinant<N> const &det) const
    {
        Determinant<N> diff_det;
        for (auto t = 0; t < N; ++t)
        {
            diff_det.repr[t] = det.repr[t] & repr[t];
        }
        return diff_det;
    }

    Determinant<N> operator^(Determinant<N> const &det) const
    {
        Determinant<N> diff_det;
        for (auto t = 0; t < N; ++t)
        {
            diff_det.repr[t] = det.repr[t] ^ repr[t];
        }
        return diff_det;
    }

    void dump(std::ostream& os) const {
        os.write(reinterpret_cast<const char*>(repr.data()), sizeof(repr));
    }

    void load(std::istream& is) {
        is.read(reinterpret_cast<char*>(repr.data()), sizeof(repr));
    }
};

template <>
int Determinant<1>::parity(Orbital i, Orbital a) const
{
    return PARITY(repr[0] & (((static_cast<det_type>(1) << i) - 1) ^
                             ((static_cast<det_type>(1) << a) - 1)));
}

template <int N>
std::array<std::array<typename Determinant<N>::det_type, Determinant<N>::size>,
           Determinant<N>::det_size * Determinant<N>::size>
    Determinant<N>::masks;

template <int N = 1>
class DeterminantDecoded : public Determinant<N>
{
public:
    std::vector<Orbital> occupied_orbitals;

    DeterminantDecoded(const Determinant<N> &det)
    {
        this->repr        = det.repr;
        occupied_orbitals = det.get_occupied_orbitals();
    }
};

// /* Hash function mapping determinants to size_t */
// template <typename T>
// inline void hash_combine(std::size_t &seed, const T &value)
// {
//     // Mixes the hash of 'value' into 'seed' using bitwise operations
//     // to reduce collisions and improve hash distribution.

//     // 0x9e3779b9 is a constant derived from the golden ratio:
//     //   0x9e3779b9 ≈ 2654435769 ≈ 2^32 / φ, where φ ≈ 1.6180339887...
//     // This constant is commonly used in hashing algorithms and number mixing techniques.

//     std::hash<T> hasher;
//     seed ^= hasher(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
// }

// template <int N>
// struct DeterminantHash
// {
//     size_t operator()(const Determinant<N> &det) const
//     {
//         size_t seed = 0;
//         for (size_t i = 0; i < N; ++i)
//         {
//             hash_combine(seed, det.repr[i]);
//         }
//         return seed;
//     }
// };

static inline uint64_t rotl64(uint64_t x, int r) {
    return (x << r) | (x >> (64 - r));
}

static inline uint64_t fmix64(uint64_t k) {
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33;
    return k;
}

template <int N>
struct DeterminantHash {
    size_t operator()(const Determinant<N>& det) const {
        uint64_t h = 0x243F6A8885A308D3ULL ^ (uint64_t)N;
        constexpr uint64_t C = 0x9ddfea08eb382d69ULL;
        for (size_t i = 0; i < N; ++i) {
            uint64_t v = (uint64_t)det.repr[i];
            uint64_t x = fmix64(v + 0x9e3779b97f4a7c15ULL * (i + 1));
            h ^= x;
            h *= C;
            h = rotl64(h, 27);
        }
        h = fmix64(h);
        if constexpr (sizeof(size_t) == 8) return (size_t)h;
        else return (size_t)(h ^ (h >> 32));
    }
};

template <int N>
struct DeterminantHashRobinhood
{
    size_t operator()(const Determinant<N> &det) const
    {
        DeterminantHash<N> hash1;
        return (robin_hood::hash_int(hash1(det)));
    }
};

/* Determinant Equal functor */
template <int N>
struct DeterminantEqual
{
    bool operator()(const Determinant<N> &det1, const Determinant<N> &det2) const
    {
        bool result = true;
        for (auto i = 0; i < N; ++i)
        {
            result = result && (det1.repr[i] == det2.repr[i]);
        }
        return result;
    }
};

#endif
