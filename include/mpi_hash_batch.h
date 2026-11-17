#ifndef CDFCI_MPI_HASH_BATCH_H
#define CDFCI_MPI_HASH_BATCH_H

#include "mpi_utils.h"
#include "lib/robin_hood/robin_hood.h"
#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

namespace cdfci_mpi
{

// Reuse transport storage across epochs. MPI's int-count interface is checked
// before narrowing; an unsupported batch aborts collectively rather than letting
// one rank throw while its peers are inside a collective.
class HashByteExchange
{
    std::vector<int> send_counts_, recv_counts_, send_displs_, recv_displs_;
    std::vector<unsigned char> packed_;
    std::vector<size_t> copy_offsets_;

    static int checked_count(size_t n, MPI_Comm comm)
    {
        if (n > static_cast<size_t>(std::numeric_limits<int>::max()))
        {
            std::cerr << "MPI hash exchange exceeds the int-count limit; "
                         "reduce the coordinate batch size." << std::endl;
            MPI_Abort(comm, 1);
            throw std::overflow_error("MPI hash exchange int-count overflow.");
        }
        return static_cast<int>(n);
    }

public:
    void exchange(const std::vector<std::vector<unsigned char>> &send,
                  std::vector<unsigned char> &recv, MPI_Comm comm)
    {
        const size_t ranks = send.size();
        if (ranks == 1)
        {
            recv.assign(send.front().begin(), send.front().end());
            return;
        }
        send_counts_.resize(ranks);
        recv_counts_.resize(ranks);
        send_displs_.resize(ranks);
        recv_displs_.resize(ranks);
        size_t total_send = 0, total_recv = 0;
        for (size_t r = 0; r < ranks; ++r)
        {
            send_displs_[r] = checked_count(total_send, comm);
            send_counts_[r] = checked_count(send[r].size(), comm);
            total_send += send[r].size();
        }
        checked_count(total_send, comm);
        check(MPI_Alltoall(send_counts_.data(), 1, MPI_INT,
                          recv_counts_.data(), 1, MPI_INT, comm),
              "MPI_Alltoall hash batch counts");
        for (size_t r = 0; r < ranks; ++r)
        {
            recv_displs_[r] = checked_count(total_recv, comm);
            total_recv += static_cast<size_t>(recv_counts_[r]);
        }
        checked_count(total_recv, comm);
        packed_.resize(total_send);
        for (size_t r = 0; r < ranks; ++r)
            if (!send[r].empty())
                std::memcpy(packed_.data() + send_displs_[r], send[r].data(), send[r].size());
        recv.resize(total_recv);
        check(MPI_Alltoallv(packed_.empty() ? nullptr : packed_.data(),
                           send_counts_.data(), send_displs_.data(), MPI_BYTE,
                           recv.empty() ? nullptr : recv.data(),
                           recv_counts_.data(), recv_displs_.data(), MPI_BYTE, comm),
              "MPI_Alltoallv hash batch records");
    }

    // Pack thread-local buffers directly into the final rank-contiguous
    // transport buffer.  The rank -> thread order is the same as the old
    // merge path, so record order and MPI displacements remain unchanged.
    void exchange_threaded(
        const std::vector<std::vector<std::vector<unsigned char>>> &thread_send,
        std::vector<unsigned char> &recv, MPI_Comm comm)
    {
        const size_t threads = thread_send.size();
        const size_t ranks = static_cast<size_t>(cdfci_mpi::size(comm));
        for (const auto &buffers : thread_send)
        {
            if (buffers.size() != ranks)
            {
                MPI_Abort(comm, 1);
                throw std::invalid_argument("Inconsistent MPI threaded buffer dimensions.");
            }
        }
        send_counts_.assign(ranks, 0);
        recv_counts_.assign(ranks, 0);
        send_displs_.assign(ranks, 0);
        recv_displs_.assign(ranks, 0);

        // The rank-major/thread-minor prefix is the same as the serial copy
        // order. Every block gets one disjoint final destination range.
        copy_offsets_.resize(ranks * threads);
        size_t total_send = 0;
        for (size_t r = 0; r < ranks; ++r)
        {
            send_displs_[r] = checked_count(total_send, comm);
            for (size_t t = 0; t < threads; ++t)
            {
                copy_offsets_[r * threads + t] = total_send;
                total_send += thread_send[t][r].size();
            }
            send_counts_[r] = checked_count(
                total_send - static_cast<size_t>(send_displs_[r]), comm);
        }
        checked_count(total_send, comm);

        packed_.resize(total_send);
        // Small messages retain a serial path so thread-team startup does not
        // dominate. Flatten rank x thread blocks: parallelism is not limited
        // to the number of destination ranks on small-node layouts.
        constexpr size_t parallel_copy_min_bytes = 1024 * 1024;
        const bool parallel_copy = threads > 1 && total_send >= parallel_copy_min_bytes;
        const size_t blocks = ranks * threads;
        if (parallel_copy)
        {
#ifdef _OPENMP
#pragma omp parallel for num_threads(static_cast<int>(threads)) schedule(static)
#endif
            for (size_t block = 0; block < blocks; ++block)
            {
                const size_t r = block / threads;
                const size_t t = block % threads;
                const auto &src = thread_send[t][r];
                if (!src.empty())
                    std::memcpy(packed_.data() + copy_offsets_[block], src.data(), src.size());
            }
        }
        else
        {
            for (size_t block = 0; block < blocks; ++block)
            {
                const auto &src = thread_send[block % threads][block / threads];
                if (!src.empty())
                    std::memcpy(packed_.data() + copy_offsets_[block], src.data(), src.size());
            }
        }

        if (ranks == 1)
        {
            recv.assign(packed_.begin(), packed_.end());
            return;
        }

        check(MPI_Alltoall(send_counts_.data(), 1, MPI_INT,
                           recv_counts_.data(), 1, MPI_INT, comm),
              "MPI_Alltoall hash batch counts");
        size_t total_recv = 0;
        for (size_t r = 0; r < ranks; ++r)
        {
            recv_displs_[r] = checked_count(total_recv, comm);
            total_recv += static_cast<size_t>(recv_counts_[r]);
        }
        checked_count(total_recv, comm);
        recv.resize(total_recv);
        check(MPI_Alltoallv(
                  packed_.empty() ? nullptr : packed_.data(),
                  send_counts_.data(), send_displs_.data(), MPI_BYTE,
                  recv.empty() ? nullptr : recv.data(), recv_counts_.data(),
                  recv_displs_.data(), MPI_BYTE, comm),
              "MPI_Alltoallv threaded hash batch records");
    }
};

// Owner-local assembly for an epoch of the distributed wavefunction hash table.
// Keep (target, selected-column, H_ji) records until BOTH dz and exact selected z
// have been accumulated. x is immutable throughout this epoch. Only flush()
// changes z; no table iterators or references survive an insertion.
template <typename W>
class HashUpdateBatch
{
    using Det = typename W::key_type;
    struct Entry
    {
        double x = 0.0;
        double old_z = 0.0;
        double dz = 0.0;
        bool nonzero_contribution = false;
    };
    using Map = robin_hood::unordered_flat_map<Det, Entry,
                                              typename W::hasher, typename W::key_equal>;
    struct Shard
    {
        Map pending;
        std::vector<size_t> offsets;
        std::vector<double> exact_z;
        std::vector<Det> kept;
        long long z_count_delta = 0;
    };
    std::vector<Shard> shards_;

public:
    static constexpr size_t record_bytes =
        determinant_bytes<Det>() + sizeof(int) + sizeof(double);

    size_t prepare(const std::vector<unsigned char> &records,
                   const std::vector<double> &dx, double scale,
                   const W &wf, int workers, std::vector<double> &exact_z)
    {
        if (records.size() % record_bytes != 0)
            throw std::runtime_error("Corrupted MPI hash batch record size.");
        workers = std::max(1, workers);
        shards_.resize(static_cast<size_t>(workers));
        for (auto &shard : shards_)
        {
            shard.pending.clear();
            shard.offsets.clear();
            shard.exact_z.assign(dx.size(), 0.0);
            shard.kept.clear();
        }

        // Stable parallel partition: count each input block by shard, assign
        // disjoint prefix-sum ranges, then fill those ranges in input order.
        // This preserves the accumulation order within every shard.
        typename W::hasher hash;
        const size_t record_count = records.size() / record_bytes;
        const int block_count = workers;
        std::vector<size_t> block_counts(
            static_cast<size_t>(block_count) * workers, 0);
        int invalid_index = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(block_count) schedule(static) reduction(|:invalid_index)
#endif
        for (int b = 0; b < block_count; ++b)
        {
            const size_t begin = record_count * static_cast<size_t>(b) / block_count;
            const size_t end = record_count * static_cast<size_t>(b + 1) / block_count;
            for (size_t i = begin; i < end; ++i)
            {
                const auto *src = records.data() + i * record_bytes;
                int p = 0;
                std::memcpy(&p, src + determinant_bytes<Det>(), sizeof(p));
                if (p < 0 || static_cast<size_t>(p) >= dx.size())
                    invalid_index = 1;
                if (workers == 1)
                    continue;
                Det det = unpack_determinant<Det>(src);
                size_t shard = robin_hood::hash_int(hash(det)) % shards_.size();
                ++block_counts[static_cast<size_t>(b) * workers + shard];
            }
        }
        // Throw on the calling thread, before any selected-index access.
        if (invalid_index)
            throw std::runtime_error("Corrupted MPI hash batch selected index.");
        std::vector<size_t> block_offsets(
            static_cast<size_t>(block_count) * workers, 0);
        for (int shard = 0; shard < workers; ++shard)
        {
            size_t prefix = 0;
            for (int b = 0; b < block_count; ++b)
            {
                block_offsets[static_cast<size_t>(b) * workers + shard] = prefix;
                prefix += block_counts[static_cast<size_t>(b) * workers + shard];
            }
            shards_[static_cast<size_t>(shard)].offsets.resize(prefix);
        }
#ifdef _OPENMP
#pragma omp parallel for num_threads(block_count) schedule(static)
#endif
        for (int b = 0; b < block_count; ++b)
        {
            if (workers == 1)
                continue;
            std::vector<size_t> local_counts(static_cast<size_t>(workers), 0);
            const size_t begin = record_count * static_cast<size_t>(b) / block_count;
            const size_t end = record_count * static_cast<size_t>(b + 1) / block_count;
            for (size_t i = begin; i < end; ++i)
            {
                const auto *src = records.data() + i * record_bytes;
                Det det = unpack_determinant<Det>(src);
                size_t shard = robin_hood::hash_int(hash(det)) % shards_.size();
                size_t pos = block_offsets[static_cast<size_t>(b) * workers + shard]
                           + local_counts[shard]++;
                shards_[shard].offsets[pos] = i * record_bytes;
            }
        }

        // Work-sharing over logical shards is correct even if OpenMP supplies
        // fewer physical threads than requested.
#ifdef _OPENMP
#pragma omp parallel for num_threads(workers) schedule(static)
#endif
        for (int s = 0; s < workers; ++s)
        {
            auto &shard = shards_[static_cast<size_t>(s)];
            const size_t count = workers == 1 ? records.size() / record_bytes
                                               : shard.offsets.size();
            for (size_t i = 0; i < count; ++i)
            {
                const size_t offset = workers == 1 ? i * record_bytes : shard.offsets[i];
                const auto *src = records.data() + offset;
                Det det = unpack_determinant<Det>(src);
                int p = 0;
                double hij = 0.0;
                std::memcpy(&p, src + determinant_bytes<Det>(), sizeof(p));
                std::memcpy(&hij, src + determinant_bytes<Det>() + sizeof(p), sizeof(hij));
                auto insertion = shard.pending.try_emplace(det);
                auto &entry = insertion.first->second;
                if (insertion.second)
                {
                    auto value = wf.local_get_value_or_zero(det);
                    entry.x = value[0];
                    entry.old_z = value[1];
                }
                double dz = hij * dx[static_cast<size_t>(p)];
                entry.dz += dz;
                // At zero threshold, cancellation must not remove a candidate
                // which the direct per-contribution path would have retained.
                entry.nonzero_contribution |= std::fabs(scale * dz) > 0.0;
                shard.exact_z[static_cast<size_t>(p)] += hij * entry.x;
            }
        }
        size_t unique = 0;
        exact_z.assign(dx.size(), 0.0);
        for (const auto &shard : shards_)
        {
            unique += shard.pending.size();
            for (size_t p = 0; p < dx.size(); ++p)
                exact_z[p] += shard.exact_z[p];
        }
        return unique;
    }

    template <typename CandidateSet>
    size_t flush(W &wf, double scale, double threshold, CandidateSet &candidates)
    {
        int workers = static_cast<int>(shards_.size());
        // Serial containers may be used with OpenMP column generation too.
#if defined(_OPENMP) && !defined(CDFCI_SOLVER_SERIAL)
#pragma omp parallel for num_threads(workers) schedule(static)
#endif
        for (int s = 0; s < workers; ++s)
        {
            auto &shard = shards_[static_cast<size_t>(s)];
            long long count_delta = 0;
            for (const auto &item : shard.pending)
            {
                const auto &entry = item.second;
                bool keep = entry.x != 0.0 || entry.old_z != 0.0
                         || (threshold == 0.0 ? entry.nonzero_contribution
                                              : std::fabs(scale * entry.dz) > threshold);
                if (keep)
                {
                    Det det = item.first;
                    count_delta += wf.local_add_z_deferred_count(det, entry.dz);
                    shard.kept.push_back(det);
                }
            }
            shard.z_count_delta = count_delta;
        }
        long long z_count_delta = 0;
        for (const auto &shard : shards_)
            z_count_delta += shard.z_count_delta;
        wf.local_commit_z_count_delta(z_count_delta);
        size_t writes = 0;
        for (const auto &shard : shards_)
        {
            writes += shard.kept.size();
            for (const auto &det : shard.kept)
                candidates.insert(det);
        }
        return writes;
    }
};

} // namespace cdfci_mpi
#endif
