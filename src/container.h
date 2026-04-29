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

#ifndef CDFCI_CONTAINER_H
#define CDFCI_CONTAINER_H

#include "lib/libcuckoo/cuckoohash_map.hh"
#include "lib/robin_hood/robin_hood.h"
#include "config.h"
#include <fstream>
#include <type_traits>
#include <unordered_map>
#include <optional>

/* Single thread hash table */

// wrapper of robin_hood::unordered_flat_map, single thread, fast
template <typename Key, typename T, typename Hash = std::hash<Key>,
          typename Pred = std::equal_to<Key>>
class ContainerRobinhood
{
public:
    using hash_map       = robin_hood::unordered_flat_map<Key, T, Hash, Pred>;
    using key_type       = typename hash_map::key_type;
    using mapped_type    = typename hash_map::mapped_type;
    using value_type     = typename hash_map::value_type;
    using size_type      = typename hash_map::size_type;
    using hasher         = Hash;
    using key_equal      = Pred;
    using iterator       = typename hash_map::iterator;
    using const_iterator = typename hash_map::const_iterator;

protected:
    hash_map data_;

public:
    /* Destructor */
    virtual ~ContainerRobinhood() {}

    /* Iterators */
    const_iterator begin() const { return data_.begin(); }
    const_iterator end() const { return data_.end(); }
    iterator       begin() { return data_.begin(); }
    iterator       end() { return data_.end(); }
    const_iterator find(const key_type &key) const { return data_.find(key); }
    iterator       find(const key_type &key) { return data_.find(key); }

    /* find_val */
    std::optional<mapped_type> find_val(const key_type &key) const
    {
        auto it = data_.find(key);
        if (it != data_.end())
        {
            return it->second;
        }
        return std::nullopt;
    }

    /* other member functions */
    size_type                 size() const { return data_.size(); }
    bool                      empty() const { return data_.empty(); }
    size_type                 capacity() const { return data_.capacity(); }
    void                      clear() { data_.clear(); }
    void                      reserve(size_type n) { data_.reserve(n); }
    std::pair<iterator, bool> insert(const value_type &keyval) { return data_.insert(keyval); }
    std::pair<iterator, bool> insert(value_type &&keyval) { return data_.insert(std::move(keyval)); }
    void                      emplace(const key_type &key, const mapped_type &val) { data_.emplace(key, val); }

    /* loop function*/
    template <typename Func>
    void loop(Func&& f) {
        for (auto& elem : data_)
            f(elem);
    }

    template <typename Func>
    void loop(Func&& f) const {
        for (const auto& elem : data_)
            f(elem);
    }

    int load_from_file(std::ifstream& file) {
        // 1) read count
        uint64_t n = 0;
        file.read(reinterpret_cast<char*>(&n), sizeof(n));
        if (!file) {
            std::cerr << "Read error: failed to read size.\n";
            return -2;
        }

        std::cout << "Loading " << n << " elements into the hash map.\n";

        // 2) read into a temporary map for strong exception safety
        hash_map tmp;
        // Reserve to avoid rehash while loading
        tmp.reserve(static_cast<size_type>(n));

        for (uint64_t i = 0; i < n; ++i) {
            key_type    k{};
            mapped_type v{};

            // Unified load interface for arbitrary Key / Value
            if (!load(file, k)) {
                std::cerr << "Read error while loading key at index " << i << ".\n";
                return -3;
            }
            if (!load(file, v)) {
                std::cerr << "Read error while loading value at index " << i << ".\n";
                return -4;
            }
            // Insert into temporary map
            auto it_ins = tmp.emplace(std::move(k), std::move(v));
        }

        if (!file.good() && !file.eof()) {
            std::cerr << "Read error: stream went bad.\n";
            return -5;
        }

        data_ = std::move(tmp);
        return 0;
    }

    int dump_to_file(std::ofstream& file) {
        const uint64_t n = static_cast<uint64_t>(data_.size());
        dump(file, n);

        for (const auto& kv : data_) {
            dump(file, kv.first);
            dump(file, kv.second);
            if (!file) {
                std::cerr << "Write error while dumping data.\n";
                return -2;
            }
        }
        return file.good() ? 0 : -3;
    }

};

// std::unordered_map, single thread, slow
template <typename Key, typename T, typename Hash = std::hash<Key>,
          typename Pred = std::equal_to<Key>>
using ContainerStd = std::unordered_map<Key, T, Hash, Pred>;

/* Concurrent hash table */

// wrapper of cuckoohash_map, multiple threads, medium speed, sensitive to the hash
// function wrapper is used to provide interfaces as std::unordered_map
// it does not have native const_iterator find.
template <typename Key, typename T, typename Hash = std::hash<Key>,
          typename Pred = std::equal_to<Key>>
class ContainerCuckoo
{
public:
    using hash_map        = libcuckoo::cuckoohash_map<Key, T, Hash, Pred>;
    using key_type        = typename hash_map::key_type;
    using mapped_type     = typename hash_map::mapped_type;
    using value_type      = typename hash_map::value_type;
    using size_type       = typename hash_map::size_type;
    using difference_type = typename hash_map::difference_type;
    using hasher          = Hash;
    using key_equal       = Pred;
    using reference       = typename hash_map::reference;
    using const_reference = typename hash_map::const_reference;
    using pointer         = typename hash_map::pointer;
    using const_pointer   = typename hash_map::const_pointer;
    using iterator        = typename hash_map::locked_table::iterator;
    using const_iterator  = typename hash_map::locked_table::const_iterator;

protected:
    hash_map data_;

public:
    /* Destructor */
    virtual ~ContainerCuckoo() {}

    /* Iterators */
    // WARNING:
    //   Iterators may be invalid as long as one thread modifies the table.
    //   Users should take care of the validity of iterators.
    const_iterator begin() const
    {
        auto           locked_table = data_.lock_table();
        const_iterator iter         = locked_table.begin();
        locked_table.unlock();
        return iter;
    }

    const_iterator end() const
    {
        auto           locked_table = data_.lock_table();
        const_iterator iter         = locked_table.end();
        locked_table.unlock();
        return iter;
    }

    iterator begin()
    {
        auto     locked_table = data_.lock_table();
        iterator iter         = locked_table.begin();
        locked_table.unlock();
        return iter;
    }

    iterator end()
    {
        auto     locked_table = data_.lock_table();
        iterator iter         = locked_table.end();
        locked_table.unlock();
        return iter;
    }

    const_iterator find(const key_type &key) const
    {
        auto           locked_table = data_.lock_table();
        const_iterator iter         = locked_table.find(key);
        locked_table.unlock();
        return iter;
    }

    iterator find(const key_type &key)
    {
        auto     locked_table = data_.lock_table();
        iterator iter         = locked_table.find(key);
        locked_table.unlock();
        return iter;
    }

    /* find_val */
    std::optional<mapped_type> find_val(const key_type &key) const
    {
        mapped_type val;
        bool        found = data_.find(key, val);
        if (found)
            return val;
        else
            return std::nullopt;
    }

    /* Other interfaces to the member functions */
    size_type size() const { return data_.size(); }
    bool      empty() const { return data_.empty(); }
    size_type capacity() const { return data_.capacity(); }
    void      clear() { data_.clear(); }
    void      reserve(size_type n) { data_.reserve(n); }

    template <typename K, typename F, typename V>
    bool upsert(K &&key, F fn, V &&val)
    {
        return data_.upsert(key, fn, val);
    }

    template <typename K, typename F>
    bool update_fn(K &&key, F fn)
    {
        return data_.update_fn(key, fn);
    }

    /* loop function */
    template <typename Func>
    void loop(Func&& f) {
        auto locked_table = data_.lock_table();
        for (auto& elem : locked_table)
            f(elem);
        locked_table.unlock();
    }

    template <typename Func>
    void loop(Func&& f) const {
        auto locked_table = data_.lock_table();
        for (const auto& elem : data_)
            f(elem);
        locked_table.unlock();
    }

    int load_from_file(std::ifstream& file) {
        // 1) read count
        uint64_t n = 0;
        file.read(reinterpret_cast<char*>(&n), sizeof(n));
        if (!file) {
            std::cerr << "Read error: failed to read size.\n";
            return -2;
        }
        std::cout << "Loading " << n << " elements into the hash map.\n";

        // 2) read into a temporary map for strong exception safety
        hash_map tmp;
        // Reserve to avoid rehash while loading
        tmp.reserve(static_cast<size_type>(n));

        for (uint64_t i = 0; i < n; ++i) {
            key_type    k{};
            mapped_type v{};

            // Unified load interface for arbitrary Key / Value
            if (!load(file, k)) {
                std::cerr << "Read error while loading key at index " << i << ".\n";
                return -3;
            }
            if (!load(file, v)) {
                std::cerr << "Read error while loading value at index " << i << ".\n";
                return -4;
            }
            // Insert into temporary map
            // use upsert: if key exists, overwrite with new value
            tmp.upsert(std::move(k),
                    [&](mapped_type& existing) { existing = v; },
                    std::move(v));
        }

        if (!file.good() && !file.eof()) {
            std::cerr << "Read error: stream went bad.\n";
            return -5;
        }

        data_ = std::move(tmp);
        return 0;
    }

    int dump_to_file(std::ofstream& file) {
        const uint64_t n = static_cast<uint64_t>(data_.size());
        dump(file, n);

        auto lt = data_.lock_table();
        for (const auto& kv : lt) {
            dump(file, kv.first);
            dump(file, kv.second);
            if (!file) {
                std::cerr << "Write error while dumping data.\n";
                return -2;
            }
        }
        lt.unlock();
        return file.good() ? 0 : -3;
    }

};

/* container using std::vector */
template <typename Key, typename T, typename Hash = std::hash<Key>,
          typename Pred = std::equal_to<Key>> // Hash is not used
class ContainerVector
{
public:
    using hash_map = std::vector<std::pair<Key, T>>;

    using key_type        = Key;
    using mapped_type     = T;
    using value_type      = std::pair<Key, T>;
    using size_type       = typename hash_map::size_type;
    using difference_type = typename hash_map::difference_type;
    using hasher          = Hash;
    using key_equal       = Pred;
    using reference       = typename hash_map::reference;
    using const_reference = typename hash_map::const_reference;
    using pointer         = typename hash_map::pointer;
    using const_pointer   = typename hash_map::const_pointer;
    using iterator        = typename hash_map::iterator;
    using const_iterator  = typename hash_map::const_iterator;

protected:
    hash_map data_;

public:
    /* Destructor */
    virtual ~ContainerVector() {}

    const_iterator begin() const { return data_.begin(); }
    const_iterator end() const { return data_.end(); }
    iterator       begin() { return data_.begin(); }
    iterator       end() { return data_.end(); }

    // interfaces to the member functions
    size_type size() const { return data_.size(); } // why lock table
    bool      empty() const { return data_.empty(); }
    size_type capacity() const { return data_.capacity(); }
    void      clear() { data_.clear(); }
    void      reserve(size_type n) { data_.reserve(n); }
    void      push_back(value_type &&val) { data_.push_back(std::move(val)); }
    void      push_back(const value_type &val) { data_.push_back(val); }
    void      emplace_back(key_type &&det, mapped_type &&val) { data_.emplace_back(std::move(det), std::move(val)); }
    void      emplace_back(const key_type &det, const mapped_type &val) { data_.emplace_back(det, val); }

    iterator find(key_type &&key)
    {
        for (auto iter = data_.begin(); iter != data_.end(); ++iter)
            if (Pred::operator()((*iter).first, key))
                return iter;
        return data_.end();
    }

    template <typename C>
    void append(C &sub_container)
    {
        data_.insert(data_.end(), sub_container.begin(), sub_container.end());
        return;
    }

    /* loop function */
    template <typename Func>
    void loop(Func&& f) {
        for (auto& elem : data_)
            f(elem);
    }

    template <typename Func>
    void loop(Func&& f) const {
        for (const auto& elem : data_)
            f(elem);
    }
};

/* Define upsert functions for all containers */
// check if the container support upsert()
template <typename T>
struct exist_upsert
{
private:
    template <typename X>
    static auto check(X x)
        -> decltype(x.upsert(typename X::key_type(), 1, typename X::mapped_type()),
                    std::true_type());

    static std::false_type check(...);

public:
    using type = decltype(check(std::declval<T>()));
    bool value = type::value;
};

// the container supports upsert
template <typename T, typename K, typename F, typename V,
          typename std::enable_if<exist_upsert<T>::type::value, int>::type = 0>
bool upsert(T &table, K &&key, F fn, V &&val)
{
    return table.upsert(std::forward<K>(key), fn, std::forward<V>(val));
}

// the container does not support upsert
// K = T::key_type, V = T::mapped_type
template <typename T, typename K, typename F, typename V,
          typename std::enable_if<!exist_upsert<T>::type::value, int>::type = 0>
bool upsert(T &table, K &&key, F fn, V &&val)
{
    auto iter = table.find(std::forward<K>(key));
    if (iter != table.end())
    {
        fn(iter->second);
        return false;
    }
    table.insert({std::forward<K>(key), std::forward<V>(val)});
    return true;
}

/* define update_fn for all containers */
// check if the container support update_fn()
template <typename T>
struct exist_update_fn
{
private:
    template <typename X>
    static auto check(X x)
        -> decltype(x.update_fn(typename X::key_type(), 1), std::true_type());

    static std::false_type check(...);

public:
    using type = decltype(check(std::declval<T>()));
    bool value = type::value;
};

// the container supports update_fn
template <typename T, typename K, typename F,
          typename std::enable_if<exist_update_fn<T>::type::value, int>::type = 0>
bool update_fn(T &table, K &&key, F fn)
{
    return table.update_fn(std::forward<K>(key), fn);
}

// the container does not support update_fn
template <typename T, typename K, typename F,
          typename std::enable_if<!exist_update_fn<T>::type::value, int>::type = 0>
bool update_fn(T &table, K &&key, F fn)
{
    auto iter = table.find(std::forward<K>(key));
    if (iter != table.end())
    {
        fn(iter->second);
        return true;
    }
    return false;
}

#endif
