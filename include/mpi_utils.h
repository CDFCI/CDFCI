#ifndef CDFCI_MPI_UTILS_H
#define CDFCI_MPI_UTILS_H

#ifndef CDFCI_ENABLE_MPI
#error "mpi_utils.h requires CDFCI_ENABLE_MPI"
#endif

#include "determinant.h"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <mpi.h>

namespace cdfci_mpi
{

inline int rank(MPI_Comm comm = MPI_COMM_WORLD)
{
    int initialized = 0;
    int finalized = 0;
    MPI_Initialized(&initialized);
    if (!initialized)
        return 0;
    MPI_Finalized(&finalized);
    if (finalized)
        return 0;
    int value = 0;
    MPI_Comm_rank(comm, &value);
    return value;
}

inline int size(MPI_Comm comm = MPI_COMM_WORLD)
{
    int initialized = 0;
    int finalized = 0;
    MPI_Initialized(&initialized);
    if (!initialized)
        return 1;
    MPI_Finalized(&finalized);
    if (finalized)
        return 1;
    int value = 1;
    MPI_Comm_size(comm, &value);
    return value;
}

inline bool is_root(MPI_Comm comm = MPI_COMM_WORLD, int root = 0)
{
    int initialized = 0;
    int finalized = 0;
    MPI_Initialized(&initialized);
    if (!initialized)
        return true;
    MPI_Finalized(&finalized);
    if (finalized)
        return true;
    return rank(comm) == root;
}

inline int env_int(const char *name)
{
    const char *text = std::getenv(name);
    if (text == nullptr || *text == '\0')
        return -1;

    char *end = nullptr;
    long value = std::strtol(text, &end, 10);
    if (end == text)
        return -1;
    return static_cast<int>(value);
}

inline int launcher_rank_hint()
{
    const char *names[] = {
        "SLURM_PROCID",
        "PMI_RANK",
        "PMIX_RANK",
        "OMPI_COMM_WORLD_RANK",
        "MV2_COMM_WORLD_RANK"};
    for (const char *name : names)
    {
        int value = env_int(name);
        if (value >= 0)
            return value;
    }
    return -1;
}

inline int launcher_size_hint()
{
    const char *names[] = {
        "SLURM_NTASKS",
        "PMI_SIZE",
        "PMIX_SIZE",
        "OMPI_COMM_WORLD_SIZE",
        "MV2_COMM_WORLD_SIZE"};
    for (const char *name : names)
    {
        int value = env_int(name);
        if (value > 0)
            return value;
    }
    return -1;
}

inline bool is_io_root(MPI_Comm comm = MPI_COMM_WORLD, int root = 0)
{
    // For stdout/stderr only.  Trust the MPI rank whenever MPI actually formed
    // a multi-rank communicator.  The launcher-rank fallback is only for the
    // mislaunch case where every process is a singleton MPI rank 0.
    if (!is_root(comm, root))
        return false;
    if (size(comm) > 1)
        return true;

    int launcher_rank = launcher_rank_hint();
    return launcher_rank < 0 || launcher_rank == root;
}

inline void check(int ierr, const char *where);

inline int openmp_max_threads()
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

inline void print_rank_thread_info(MPI_Comm comm = MPI_COMM_WORLD, int root = 0)
{
    int my_rank = rank(comm);
    int nranks = size(comm);
    int threads = openmp_max_threads();

    std::vector<int> all_threads;
    if (my_rank == root)
        all_threads.resize(static_cast<size_t>(nranks), 1);

    check(MPI_Gather(&threads,
                     1,
                     MPI_INT,
                     all_threads.empty() ? nullptr : all_threads.data(),
                     1,
                     MPI_INT,
                     root,
                     comm),
          "MPI_Gather OpenMP threads");

    char processor_name[MPI_MAX_PROCESSOR_NAME] = {};
    int processor_name_len = 0;
    check(MPI_Get_processor_name(processor_name, &processor_name_len),
          "MPI_Get_processor_name");
    if (processor_name_len < MPI_MAX_PROCESSOR_NAME)
        processor_name[processor_name_len] = '\0';
    else
        processor_name[MPI_MAX_PROCESSOR_NAME - 1] = '\0';

    std::vector<char> all_names;
    if (my_rank == root)
        all_names.resize(static_cast<size_t>(nranks) * MPI_MAX_PROCESSOR_NAME, '\0');

    check(MPI_Gather(processor_name,
                     MPI_MAX_PROCESSOR_NAME,
                     MPI_CHAR,
                     all_names.empty() ? nullptr : all_names.data(),
                     MPI_MAX_PROCESSOR_NAME,
                     MPI_CHAR,
                     root,
                     comm),
          "MPI_Gather processor names");

    if (!is_io_root(comm, root))
        return;

    std::cout << "MPI: enabled with " << nranks << " ranks" << std::endl;
    for (int r = 0; r < nranks; ++r)
    {
        const char *name = all_names.data()
                         + static_cast<size_t>(r) * MPI_MAX_PROCESSOR_NAME;
        std::cout << "MPI rank " << std::setw(4) << r
                  << ": " << all_threads[static_cast<size_t>(r)]
                  << " OpenMP thread(s), host " << name << std::endl;
    }

    int launcher_size = launcher_size_hint();
    if (nranks == 1 && launcher_size > 1)
    {
        std::cout << "Warning: launcher requested " << launcher_size
                  << " tasks, but MPI_COMM_WORLD has size 1. "
                  << "Check the MPI launcher and runtime libraries on every node."
                  << std::endl;
    }
}

inline void check(int ierr, const char *where)
{
    if (ierr != MPI_SUCCESS)
        throw std::runtime_error(std::string("MPI call failed at ") + where);
}

template <typename Det>
constexpr int determinant_bytes()
{
    return static_cast<int>(sizeof(typename Det::det_type) * Det::size);
}

template <typename Det>
void pack_determinant(const Det &det, unsigned char *dst)
{
    std::memcpy(dst, det.repr.data(), determinant_bytes<Det>());
}

template <typename Det>
Det unpack_determinant(const unsigned char *src)
{
    Det det;
    std::memcpy(det.repr.data(), src, determinant_bytes<Det>());
    return det;
}

template <int N>
int owner(const Determinant<N> &det, int nranks)
{
    DeterminantHash<N> hash;
    return static_cast<int>(hash(det) % static_cast<size_t>(nranks));
}

template <typename Det>
bool determinant_less(const Det &lhs, const Det &rhs)
{
    for (int i = 0; i < Det::size; ++i)
    {
        if (lhs.repr[i] < rhs.repr[i])
            return true;
        if (rhs.repr[i] < lhs.repr[i])
            return false;
    }
    return false;
}

template <typename Det>
struct Candidate
{
    Det det;
    double score = 0.0;
    bool valid = false;
};

template <typename Det>
int candidate_record_bytes()
{
    return static_cast<int>(sizeof(int) + sizeof(double) + determinant_bytes<Det>());
}

template <typename Det>
void pack_candidate(const Candidate<Det> &candidate, unsigned char *dst)
{
    int valid = candidate.valid ? 1 : 0;
    std::memcpy(dst, &valid, sizeof(int));
    std::memcpy(dst + sizeof(int), &candidate.score, sizeof(double));
    pack_determinant(candidate.det, dst + sizeof(int) + sizeof(double));
}

template <typename Det>
Candidate<Det> unpack_candidate(const unsigned char *src)
{
    Candidate<Det> candidate;
    int valid = 0;
    std::memcpy(&valid, src, sizeof(int));
    std::memcpy(&candidate.score, src + sizeof(int), sizeof(double));
    candidate.det = unpack_determinant<Det>(src + sizeof(int) + sizeof(double));
    candidate.valid = (valid != 0);
    return candidate;
}

template <typename Det>
std::vector<Candidate<Det>>
allgather_candidates(const std::vector<Candidate<Det>> &local_top,
                     size_t block_size,
                     MPI_Comm comm = MPI_COMM_WORLD)
{
    const int record_bytes = candidate_record_bytes<Det>();
    const int send_bytes = static_cast<int>(block_size) * record_bytes;
    std::vector<unsigned char> send(send_bytes, 0);

    for (size_t i = 0; i < local_top.size() && i < block_size; ++i)
        pack_candidate(local_top[i], send.data() + static_cast<int>(i) * record_bytes);

    int nranks = size(comm);
    std::vector<unsigned char> recv(static_cast<size_t>(send_bytes) * nranks, 0);
    check(MPI_Allgather(send.data(), send_bytes, MPI_BYTE,
                        recv.data(), send_bytes, MPI_BYTE, comm),
          "MPI_Allgather candidates");

    std::vector<Candidate<Det>> result;
    result.reserve(static_cast<size_t>(nranks) * block_size);
    for (int r = 0; r < nranks; ++r)
    {
        for (size_t i = 0; i < block_size; ++i)
        {
            const auto *ptr = recv.data()
                            + static_cast<size_t>(r) * send_bytes
                            + i * record_bytes;
            auto candidate = unpack_candidate<Det>(ptr);
            if (candidate.valid)
                result.push_back(candidate);
        }
    }
    return result;
}

template <typename Det>
std::vector<Det>
select_global_top(std::vector<Candidate<Det>> candidates, size_t block_size)
{
    std::sort(candidates.begin(), candidates.end(),
              [](const auto &lhs, const auto &rhs) {
                  if (lhs.score != rhs.score)
                      return lhs.score > rhs.score;
                  return determinant_less(rhs.det, lhs.det);
              });

    std::vector<Det> selected;
    selected.reserve(block_size);
    for (const auto &candidate : candidates)
    {
        bool seen = false;
        for (const auto &det : selected)
        {
            if (!determinant_less(det, candidate.det)
                && !determinant_less(candidate.det, det))
            {
                seen = true;
                break;
            }
        }
        if (!seen)
            selected.push_back(candidate.det);
        if (selected.size() == block_size)
            break;
    }
    return selected;
}

template <typename Det>
void bcast_determinants(std::vector<Det> &dets,
                        int root = 0,
                        MPI_Comm comm = MPI_COMM_WORLD)
{
    int count = static_cast<int>(dets.size());
    check(MPI_Bcast(&count, 1, MPI_INT, root, comm), "MPI_Bcast determinant count");
    dets.resize(static_cast<size_t>(count));
    if (count == 0)
        return;

    const int bytes_per_det = determinant_bytes<Det>();
    std::vector<unsigned char> buffer(static_cast<size_t>(count) * bytes_per_det, 0);
    if (rank(comm) == root)
    {
        for (int i = 0; i < count; ++i)
            pack_determinant(dets[static_cast<size_t>(i)],
                             buffer.data() + static_cast<size_t>(i) * bytes_per_det);
    }
    check(MPI_Bcast(buffer.data(),
                    static_cast<int>(buffer.size()),
                    MPI_BYTE,
                    root,
                    comm),
          "MPI_Bcast determinants");
    if (rank(comm) != root)
    {
        for (int i = 0; i < count; ++i)
            dets[static_cast<size_t>(i)] =
                unpack_determinant<Det>(buffer.data() + static_cast<size_t>(i) * bytes_per_det);
    }
}

inline void allreduce_sum(std::vector<double> &local,
                          std::vector<double> &global,
                          MPI_Comm comm = MPI_COMM_WORLD)
{
    global.assign(local.size(), 0.0);
    if (local.empty())
        return;
    check(MPI_Allreduce(local.data(),
                        global.data(),
                        static_cast<int>(local.size()),
                        MPI_DOUBLE,
                        MPI_SUM,
                        comm),
          "MPI_Allreduce double vector");
}

inline unsigned long long allreduce_sum_ull(unsigned long long local,
                                            MPI_Comm comm = MPI_COMM_WORLD)
{
    unsigned long long global = 0;
    check(MPI_Allreduce(&local,
                        &global,
                        1,
                        MPI_UNSIGNED_LONG_LONG,
                        MPI_SUM,
                        comm),
          "MPI_Allreduce unsigned long long");
    return global;
}

inline void bcast_scalars(double &xx,
                          double &xz,
                          double &energy,
                          double &scale,
                          int root = 0,
                          MPI_Comm comm = MPI_COMM_WORLD)
{
    double values[4] = {xx, xz, energy, scale};
    check(MPI_Bcast(values, 4, MPI_DOUBLE, root, comm), "MPI_Bcast scalars");
    xx = values[0];
    xz = values[1];
    energy = values[2];
    scale = values[3];
}

} // namespace cdfci_mpi

#endif
