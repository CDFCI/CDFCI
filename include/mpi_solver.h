#ifndef CDFCI_MPI_SOLVER_H
#define CDFCI_MPI_SOLVER_H

#ifndef CDFCI_ENABLE_MPI
#error "mpi_solver.h requires CDFCI_ENABLE_MPI"
#endif

#include "config.h"
#include "hamiltonian.h"
#include "lib/robin_hood/robin_hood.h"
#include "mpi_utils.h"
#include "mpi_hash_batch.h"
#include "option.h"
#include "wavefunction.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

// MPI timer labels accumulate wall-clock time over all iterations.  The final
// report uses MPI_MAX, so each printed value is the slowest rank for that stage.
enum class MPITimerLabel : size_t
{
    SelectLocal = 0,   // Score local candidates and keep the rank-local top block.
    SelectComm,        // Gather local top lists, select global top block, broadcast it.
    GetSelected,       // Read selected x/z values from owner ranks and allreduce them.
    BlockSolve,        // Build H_II on root, solve the coordinate update, broadcast dx.
    Enum,              // Enumerate connected determinants, evaluate H_ij, and pack records.
    RouteComm,         // Route column contributions to determinant-owner ranks.
    PendingUpdate,     // Accumulate local pending dz updates and exact selected-z pieces.
    HashUpdate,        // Apply pending dz values to the local wavefunction hash table.
    ExactZAllreduce,   // Allreduce exact selected z values and overwrite selected entries.
    Candidate,         // Build the next local candidate set K.
    ScalarUpdate,      // Update/broadcast xx, xz, energy, and the global scale factor.
    Count
};

class MPITimer
{
public:
    MPITimer()
    {
        times_.fill(0.0);
    }

    template <typename Func>
    void measure(MPITimerLabel label, Func &&func)
    {
        auto start = std::chrono::high_resolution_clock::now();
        func();
        auto end = std::chrono::high_resolution_clock::now();
        times_[static_cast<size_t>(label)] +=
            std::chrono::duration<double>(end - start).count();
    }

    void print_results(MPI_Comm comm, int root = 0) const
    {
        constexpr size_t nlabels = static_cast<size_t>(MPITimerLabel::Count);
        static const std::array<const char *, nlabels> labels = {
            "T_select_local",
            "T_select_comm",
            "T_get_selected",
            "T_block_solve",
            "T_enum_hij_pack",
            "T_route_comm",
            "T_pending_update",
            "T_hash_update",
            "T_exact_z_allreduce",
            "T_candidate",
            "T_scalar_update"};

        const double local_total =
            std::accumulate(times_.begin(), times_.end(), 0.0);
        double max_total = 0.0;
        cdfci_mpi::check(MPI_Reduce(&local_total,
                                     &max_total,
                                     1,
                                     MPI_DOUBLE,
                                     MPI_MAX,
                                     root,
                                     comm),
                          "MPI_Reduce timer total");

        std::array<double, nlabels> min_times = {};
        std::array<double, nlabels> sum_times = {};
        std::array<double, nlabels> max_times = {};
        cdfci_mpi::check(MPI_Reduce(times_.data(),
                                     min_times.data(),
                                     static_cast<int>(nlabels),
                                     MPI_DOUBLE,
                                     MPI_MIN,
                                     root,
                                     comm),
                          "MPI_Reduce timer mins");
        cdfci_mpi::check(MPI_Reduce(times_.data(),
                                     sum_times.data(),
                                     static_cast<int>(nlabels),
                                     MPI_DOUBLE,
                                     MPI_SUM,
                                     root,
                                     comm),
                          "MPI_Reduce timer sums");
        cdfci_mpi::check(MPI_Reduce(times_.data(),
                                     max_times.data(),
                                     static_cast<int>(nlabels),
                                     MPI_DOUBLE,
                                     MPI_MAX,
                                     root,
                                     comm),
                          "MPI_Reduce timer maxes");

        if (cdfci_mpi::is_io_root(comm, root))
        {
            std::cout << "MPI timing summary (seconds over all iterations)"
                      << std::endl;
            std::cout << std::setw(24) << "Stage"
                      << std::setw(14) << "Min"
                      << std::setw(14) << "Avg"
                      << std::setw(14) << "Max"
                      << std::setw(10) << "Max %"
                      << std::endl;
        }

        for (size_t i = 0; i < labels.size(); ++i)
        {
            if (cdfci_mpi::is_io_root(comm, root))
            {
                double avg_time = sum_times[i] / cdfci_mpi::size(comm);
                double pct = (max_total > 0.0) ? 100.0 * max_times[i] / max_total : 0.0;
                std::cout << std::setw(24) << labels[i]
                          << std::setw(14) << std::fixed << std::setprecision(2)
                          << min_times[i]
                          << std::setw(14) << avg_time
                          << std::setw(14) << max_times[i]
                          << std::setw(9) << std::setprecision(1) << pct
                          << "%"
                          << std::endl;
            }
        }
    }

private:
    std::array<double, static_cast<size_t>(MPITimerLabel::Count)> times_;
};

enum class MPICounterLabel : size_t
{
    RoutedSendBytes = 0,
    RoutedRecvBytes,
    RoutedSendRecords,
    RoutedRecvRecords,
    PendingUniqueDets,
    LocalCandidateDets,
    Count
};

class MPICounters
{
public:
    MPICounters()
    {
        values_.fill(0);
    }

    void add(MPICounterLabel label, unsigned long long value)
    {
        values_[static_cast<size_t>(label)] += value;
    }

    void print_results(MPI_Comm comm, int root = 0) const
    {
        constexpr size_t nlabels = static_cast<size_t>(MPICounterLabel::Count);
        static const std::array<const char *, nlabels> labels = {
            "routed_send_bytes",
            "routed_recv_bytes",
            "routed_send_records",
            "routed_recv_records",
            "pending_unique_dets",
            "local_candidate_dets"};

        std::array<unsigned long long, nlabels> global_values = {};
        std::array<unsigned long long, nlabels> max_values = {};
        cdfci_mpi::check(MPI_Reduce(values_.data(),
                                     global_values.data(),
                                     static_cast<int>(nlabels),
                                     MPI_UNSIGNED_LONG_LONG,
                                     MPI_SUM,
                                     root,
                                     comm),
                          "MPI_Reduce counter sums");
        cdfci_mpi::check(MPI_Reduce(values_.data(),
                                     max_values.data(),
                                     static_cast<int>(nlabels),
                                     MPI_UNSIGNED_LONG_LONG,
                                     MPI_MAX,
                                     root,
                                     comm),
                          "MPI_Reduce counter maxes");

        if (cdfci_mpi::is_io_root(comm, root))
        {
            std::cout << "MPI data summary (cumulative over all iterations)"
                      << std::endl;
            std::cout << std::setw(24) << "Counter"
                      << std::setw(18) << "Global sum"
                      << std::setw(18) << "Avg/rank"
                      << std::setw(18) << "Max rank"
                      << std::endl;
        }

        for (size_t i = 0; i < labels.size(); ++i)
        {
            if (cdfci_mpi::is_io_root(comm, root))
            {
                double avg_rank = static_cast<double>(global_values[i]) / cdfci_mpi::size(comm);
                if (i == static_cast<size_t>(MPICounterLabel::RoutedSendBytes)
                    || i == static_cast<size_t>(MPICounterLabel::RoutedRecvBytes))
                {
                    const double mib = 1024.0 * 1024.0;
                    std::cout << std::setw(24) << labels[i]
                              << std::setw(18) << std::fixed << std::setprecision(2)
                              << static_cast<double>(global_values[i]) / mib
                              << std::setw(18) << avg_rank / mib
                              << std::setw(18) << static_cast<double>(max_values[i]) / mib
                              << " MiB"
                              << std::endl;
                }
                else
                {
                    std::cout << std::setw(24) << labels[i]
                              << std::setw(18) << global_values[i]
                              << std::setw(18) << std::fixed << std::setprecision(1)
                              << avg_rank
                              << std::setw(18) << max_values[i]
                              << std::endl;
                }
            }
        }
    }

private:
    std::array<unsigned long long, static_cast<size_t>(MPICounterLabel::Count)> values_;
};

template <typename H, typename W>
class MPICDFCISolver
{
    static_assert(W::NSTATES_val == 1,
                  "The first MPI implementation supports only ground-state CDFCI.");

public:
    using det_type = typename W::key_type;
    using mapped_type = typename W::mapped_type;
    using data_type = typename W::data_type;
    using hasher = typename W::hasher;
    using key_equal = typename W::key_equal;

    Option opt = {{"num_iterations", 10},
                  {"report_interval", 1},
                  {"num_coordinates", 1},
                  {"coordinate_update", "auto"},
                  {"ref_det_occ", Option::array()},
                  {"z_threshold", 0.0},
                  {"z_threshold_search", false},
                  {"hash_updates", "buffered"},
                  {"max_wavefunction_size", 1000000},
                  {"max_load_factor", 0.79},
                  {"stopping_dx_threshold", 1e-10},
                  {"stopping_dx_damping_factor", 0.0},
                  {"load_configuration", false},
                  {"load_configuration_location", "./wavefunction.in"},
                  {"recalculate_z_on_load", false},
                  {"save_configuration", false},
                  {"save_configuration_location", "./wavefunction.out"},
                  {"verbose", 2}};

    size_t num_iter = 10;
    size_t report_interval = 1;
    size_t num_coordinates = 1;
    std::string coord_update_str = "ls";
    std::vector<Orbital> ref_det_occ;
    NumericalType z_threshold = 0.0;
    bool z_threshold_search = false;
    bool buffered_hash_updates = true;
    size_t max_wavefunction_size = 1000000;
    NumericalType max_load_factor = 0.79;
    NumericalType stopping_dx_damping_factor = 0.0;
    NumericalType stopping_dx_threshold = 1e-10;
    int verbose = 2;

    NumericalType ground_state_energy = 0.0;
    size_t iterations = 0;
    std::vector<size_t> x_size_history;
    std::vector<size_t> z_size_history;
    std::vector<double> time_history;

    MPI_Comm comm = MPI_COMM_WORLD;
    int mpi_rank = 0;
    int mpi_size = 1;
    int root = 0;
    int omp_threads = 1;
    int total_omp_workers = 1;
    std::vector<int> rank_omp_threads;
    std::vector<int> rank_thread_offsets;

    MPICDFCISolver() = default;

    MPICDFCISolver(Option option, MPI_Comm comm_ = MPI_COMM_WORLD)
        : comm(comm_)
    {
        opt.merge_patch(option);
        init();
        validate_option();
    }

    void init()
    {
        mpi_rank = cdfci_mpi::rank(comm);
        mpi_size = cdfci_mpi::size(comm);

        num_iter = opt["num_iterations"];
        report_interval = opt["report_interval"];
        num_coordinates = opt["num_coordinates"];
        ref_det_occ = opt.at("ref_det_occ").get<std::vector<int>>();
        z_threshold = opt["z_threshold"];
        z_threshold_search = opt["z_threshold_search"];
        const auto hash_updates = opt.at("hash_updates").get<std::string>();
        if (hash_updates != "buffered" && hash_updates != "direct")
            throw std::invalid_argument("mpi_cdfci.hash_updates must be buffered or direct.");
        buffered_hash_updates = hash_updates == "buffered";
        max_wavefunction_size = opt["max_wavefunction_size"];
        max_load_factor = opt["max_load_factor"];
        stopping_dx_damping_factor = opt["stopping_dx_damping_factor"];
        stopping_dx_threshold = opt["stopping_dx_threshold"];
        verbose = opt["verbose"];

        coord_update_str = opt.value("coordinate_update", std::string("auto"));
        if (coord_update_str == "auto")
            coord_update_str = (num_coordinates == 1) ? "ls" : "eig";

        configure_openmp_runtime();
        refresh_parallel_layout();
    }

    void validate_option() const
    {
        if (report_interval < 1)
            throw std::invalid_argument("mpi_cdfci.report_interval should be at least 1.");
        if (num_iter < report_interval)
            throw std::invalid_argument("mpi_cdfci.num_iterations should be no smaller than report_interval.");
        if (num_coordinates < 1)
            throw std::invalid_argument("mpi_cdfci.num_coordinates should be at least 1.");
        if (z_threshold < 0)
            throw std::invalid_argument("mpi_cdfci.z_threshold should be nonnegative.");
        if (stopping_dx_threshold < 0)
            throw std::invalid_argument("mpi_cdfci.stopping_dx_threshold should be nonnegative.");
        if (stopping_dx_damping_factor < 0 || stopping_dx_damping_factor > 1)
            throw std::invalid_argument("mpi_cdfci.stopping_dx_damping_factor should be between 0 and 1.");
        if (max_load_factor <= 0 || max_load_factor >= 1)
            throw std::invalid_argument("mpi_cdfci.max_load_factor should be between 0 and 1.");
        if (coord_update_str == "ls" && num_coordinates != 1)
            throw std::invalid_argument("mpi_cdfci.coordinate_update = ls is only valid for num_coordinates = 1.");
        if (coord_update_str != "ls" && coord_update_str != "eig")
            throw std::invalid_argument("mpi_cdfci.coordinate_update must be auto, ls, or eig.");
        if (z_threshold_search)
            throw std::invalid_argument("mpi_cdfci.z_threshold_search is not supported in the pure MPI baseline.");
        if (opt.value("load_configuration", false) || opt.value("save_configuration", false))
            throw std::invalid_argument("mpi_cdfci load/save wavefunction is not implemented in the first MPI baseline.");
    }

    void print_solver_header() const
    {
        if (verbose < 2 || !cdfci_mpi::is_io_root(comm, root))
            return;

        std::cout << "Solver information" << std::endl;
        std::cout << "-------------------" << std::endl;
        Option opt_ = {{"mpi_cdfci", ""}};
        opt_["mpi_cdfci"] = opt;
        opt_["mpi_cdfci"]["coordinate_update"] = coord_update_str;
        std::cout << opt_.dump(4) << std::endl << std::endl;
    }

    W solve(H &h)
    {
        W local_wf(max_wavefunction_size, max_load_factor);
        local_wf.num_states = 1;

        det_type hf = get_ref_det(h);
        std::vector<det_type> local_candidates;
        if (cdfci_mpi::owner(hf, mpi_size) == mpi_rank)
            local_candidates.push_back(hf);

        if (verbose > 0 && cdfci_mpi::is_io_root(comm, root))
        {
            std::cout << "MPI CDFCI calculation" << std::endl;
            std::cout << "---------------------" << std::endl;
            int slurm_nodes = cdfci_mpi::env_int("SLURM_NNODES");
            int slurm_cpus_per_task = cdfci_mpi::env_int("SLURM_CPUS_PER_TASK");
            std::cout << "MPI ranks: " << mpi_size << std::endl;
            std::cout << "Slurm nodes: ";
            if (slurm_nodes > 0)
                std::cout << slurm_nodes;
            else
                std::cout << "not provided";
            std::cout << std::endl;
            std::cout << "Slurm CPUs per MPI rank: ";
            if (slurm_cpus_per_task > 0)
                std::cout << slurm_cpus_per_task;
            else
                std::cout << "not provided";
            std::cout << std::endl;
            std::cout << "OpenMP threads on this rank: " << omp_threads << std::endl;
            std::cout << "Total MPI/OpenMP work units: " << total_omp_workers << std::endl;
            std::cout << "Resource control: nodes/ranks/threads come from sbatch, srun, and OpenMP environment." << std::endl;
            if (!rank_omp_threads.empty())
            {
                std::cout << "Rank/thread layout:";
                for (int r = 0; r < mpi_size; ++r)
                    std::cout << " rank" << r << "="
                              << rank_omp_threads[static_cast<size_t>(r)];
                std::cout << std::endl;
            }
            int launcher_size = cdfci_mpi::launcher_size_hint();
            if (mpi_size == 1 && launcher_size > 1)
            {
                std::cout << "Warning: launcher requested " << launcher_size
                          << " tasks, but MPI_COMM_WORLD has size 1. "
                          << "Use mpirun -np N or the Slurm MPI mode matching "
                          << "the loaded MPI module." << std::endl;
            }
            std::cout << "Block size: " << num_coordinates << std::endl;
            std::cout << "Reference determinant occupied spin-orbitals:" << std::endl;
            for (auto &orb : hf.get_occupied_orbitals())
                std::cout << orb << "  ";
            std::cout << std::endl;
            std::cout << "Reference energy: " << std::fixed << std::setprecision(10)
                      << h.get_diagonal(hf) << std::endl
                      << std::endl;
        }

        double xx = 0.0;
        double xz = 0.0;
        double energy = 0.0;
        double scale = 1.0;
        double dx_accumulated = 1.0;
        bool early_stop = false;

        MPITimer timer;
        MPICounters counters;
        auto time_start = std::chrono::high_resolution_clock::now();

        using PendingMap =
            robin_hood::unordered_flat_map<det_type, double, hasher, key_equal>;
        using CandidateSet =
            robin_hood::unordered_flat_set<det_type, hasher, key_equal>;
        PendingMap pending_dz;
        CandidateSet next_seen;
        std::vector<double> local_exact_z_I;
        std::vector<std::vector<unsigned char>> routed_send_buffers;
        std::vector<unsigned char> routed_recv_buffer;
        cdfci_mpi::HashUpdateBatch<W> hash_batch;
        cdfci_mpi::HashByteExchange hash_exchange;
        std::vector<std::vector<std::vector<unsigned char>>> thread_send_buffers;

        output_header();

        for (size_t outer = 0; outer < num_iter / report_interval; ++outer)
        {
            for (size_t inner = 0; inner < report_interval; ++inner)
            {
                iterations = outer * report_interval + inner + 1;

                std::vector<cdfci_mpi::Candidate<det_type>> local_top;
                timer.measure(MPITimerLabel::SelectLocal, [&]() {
                    local_top = select_local_top(local_wf, local_candidates, xx, scale);
                });

                std::vector<det_type> selected;
                timer.measure(MPITimerLabel::SelectComm, [&]() {
                    auto all_candidates =
                        cdfci_mpi::allgather_candidates(local_top, num_coordinates, comm);
                    if (mpi_rank == root)
                        selected = cdfci_mpi::select_global_top(all_candidates, num_coordinates);
                    cdfci_mpi::bcast_determinants(selected, root, comm);
                });

                if (selected.empty())
                {
                    early_stop = true;
                    break;
                }

                std::vector<double> x_I;
                std::vector<double> z_I;
                timer.measure(MPITimerLabel::GetSelected, [&]() {
                    get_selected_values(local_wf, selected, x_I, z_I);
                });

                std::vector<double> delta_x_I(selected.size(), 0.0);
                std::vector<double> H_II;
                double new_scale = scale;
                timer.measure(MPITimerLabel::BlockSolve, [&]() {
                    if (mpi_rank == root)
                    {
                        H_II = build_hii(h, selected);
                        new_scale = solve_block_update(H_II, x_I, z_I, xx, xz, delta_x_I);
                    }
                    cdfci_mpi::check(MPI_Bcast(&new_scale, 1, MPI_DOUBLE, root, comm),
                                      "MPI_Bcast new scale");
                    cdfci_mpi::check(MPI_Bcast(delta_x_I.data(),
                                               static_cast<int>(delta_x_I.size()),
                                               MPI_DOUBLE,
                                               root,
                                               comm),
                                      "MPI_Bcast delta_x_I");
                });

                for (size_t p = 0; p < selected.size(); ++p)
                {
                    if (cdfci_mpi::owner(selected[p], mpi_size) == mpi_rank)
                    {
                        det_type det = selected[p];
                        local_wf.local_add_x(det, delta_x_I[p]);
                    }
                }

                pending_dz.clear();
                next_seen.clear();
                local_exact_z_I.assign(selected.size(), 0.0);
                // Keep existing byte elements alive across epochs.  The
                // exchange resizes the vector to the exact receive length and
                // MPI_Alltoallv overwrites every byte in [0, total_recv).
                // Clearing here would make resize() value-initialize the full
                // receive payload on every iteration.

                timer.measure(MPITimerLabel::Enum, [&]() {
                    build_routed_column_buffers(h, selected, routed_send_buffers, thread_send_buffers);
                });
                unsigned long long local_send_bytes = 0;
                if (omp_threads > 1)
                {
                    for (const auto &thread_buffers : thread_send_buffers)
                        for (const auto &buffer : thread_buffers)
                            local_send_bytes += static_cast<unsigned long long>(buffer.size());
                }
                else
                {
                    for (const auto &buffer : routed_send_buffers)
                        local_send_bytes += static_cast<unsigned long long>(buffer.size());
                }
                unsigned long long record_bytes =
                    static_cast<unsigned long long>(routed_column_record_bytes());
                counters.add(MPICounterLabel::RoutedSendBytes, local_send_bytes);
                counters.add(MPICounterLabel::RoutedSendRecords,
                             record_bytes == 0 ? 0 : local_send_bytes / record_bytes);

                timer.measure(MPITimerLabel::RouteComm, [&]() {
#ifdef _OPENMP
                    if (omp_threads > 1)
                        hash_exchange.exchange_threaded(thread_send_buffers, routed_recv_buffer, comm);
                    else
                        hash_exchange.exchange(routed_send_buffers, routed_recv_buffer, comm);
#else
                    hash_exchange.exchange(routed_send_buffers, routed_recv_buffer, comm);
#endif
                });
                unsigned long long local_recv_bytes =
                    static_cast<unsigned long long>(routed_recv_buffer.size());
                counters.add(MPICounterLabel::RoutedRecvBytes, local_recv_bytes);
                counters.add(MPICounterLabel::RoutedRecvRecords,
                             record_bytes == 0 ? 0 : local_recv_bytes / record_bytes);

                size_t pending_unique_count = 0;
                timer.measure(MPITimerLabel::PendingUpdate, [&]() {
                    if (buffered_hash_updates)
                    {
                        pending_unique_count = hash_batch.prepare(routed_recv_buffer,
                            delta_x_I, new_scale, local_wf, omp_threads, local_exact_z_I);
                    }
                    else if (z_threshold == 0.0)
                    {
                        size_t recv_records =
                            record_bytes == 0
                                ? 0
                                : routed_recv_buffer.size()
                                      / static_cast<size_t>(record_bytes);
                        next_seen.reserve(recv_records + selected.size());
                        consume_routed_column_buffer_direct(routed_recv_buffer,
                                                            delta_x_I,
                                                            new_scale,
                                                            local_wf,
                                                            local_exact_z_I,
                                                            next_seen);
                        pending_unique_count = next_seen.size();
                    }
                    else
                    {
                        consume_routed_column_buffer(routed_recv_buffer,
                                                     delta_x_I,
                                                     local_wf,
                                                     pending_dz,
                                                     local_exact_z_I);
                        pending_unique_count = pending_dz.size();
                    }
                });
                counters.add(MPICounterLabel::PendingUniqueDets,
                             static_cast<unsigned long long>(pending_unique_count));

                timer.measure(MPITimerLabel::HashUpdate, [&]() {
                    if (buffered_hash_updates)
                    {
                        hash_batch.flush(local_wf, new_scale, z_threshold, next_seen);
                        return;
                    }
                    if (z_threshold == 0.0)
                        return;
                    next_seen.reserve(pending_dz.size() + selected.size());
                    for (auto &entry : pending_dz)
                    {
                        det_type target = entry.first;
                        double dz = entry.second;
                        mapped_type old_value = local_wf.local_get_value_or_zero(target);
                        bool keep = old_value[0] != 0.0
                                 || old_value[1] != 0.0
                                 || std::fabs(new_scale * dz) > z_threshold;
                        if (keep)
                        {
                            local_wf.local_add_z(target, dz);
                            next_seen.insert(target);
                        }
                    }
                });

                std::vector<double> exact_z_I;
                timer.measure(MPITimerLabel::ExactZAllreduce, [&]() {
                    cdfci_mpi::allreduce_sum(local_exact_z_I, exact_z_I, comm);
                    for (size_t p = 0; p < selected.size(); ++p)
                    {
                        if (cdfci_mpi::owner(selected[p], mpi_size) == mpi_rank)
                        {
                            det_type det = selected[p];
                            local_wf.local_set_z(det, exact_z_I[p]);
                            next_seen.insert(det);
                        }
                    }
                });

                timer.measure(MPITimerLabel::ScalarUpdate, [&]() {
                    if (mpi_rank == root)
                    {
                        double x_dot_d = dot(x_I, delta_x_I);
                        double d_dot_d = dot(delta_x_I, delta_x_I);
                        double d_dot_z = dot(delta_x_I, exact_z_I);
                        double d_h_d = quadratic_form(delta_x_I, H_II);
                        xx = xx + 2.0 * x_dot_d + d_dot_d;
                        xz = xz + 2.0 * d_dot_z - d_h_d;
                        scale = new_scale;
                        energy = xz / xx;
                    }
                    cdfci_mpi::bcast_scalars(xx, xz, energy, scale, root, comm);
                    local_wf.set_global_scalars(xx, xz, static_cast<ScaleFactorType>(scale));
                });

                timer.measure(MPITimerLabel::Candidate, [&]() {
                    local_candidates.clear();
                    local_candidates.reserve(next_seen.size());
                    for (const auto &det : next_seen)
                        local_candidates.push_back(det);
                });
                counters.add(MPICounterLabel::LocalCandidateDets,
                             static_cast<unsigned long long>(local_candidates.size()));

                if (check_convergence(delta_x_I, scale, dx_accumulated))
                {
                    early_stop = true;
                    break;
                }

                local_wf.check_overflow();
            }

            unsigned long long local_x_size =
                static_cast<unsigned long long>(local_wf.size_x());
            unsigned long long local_z_size =
                static_cast<unsigned long long>(local_wf.size_z());
            auto global_x_size = cdfci_mpi::allreduce_sum_ull(local_x_size, comm);
            auto global_z_size = cdfci_mpi::allreduce_sum_ull(local_z_size, comm);
            auto global_candidate_size = cdfci_mpi::allreduce_sum_ull(
                static_cast<unsigned long long>(local_candidates.size()), comm);

            auto time_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_elapsed = time_end - time_start;

            output_line(iterations,
                        energy,
                        dx_accumulated,
                        static_cast<size_t>(global_x_size),
                        static_cast<size_t>(global_z_size),
                        static_cast<size_t>(global_candidate_size),
                        time_elapsed.count());

            if (early_stop)
                break;
        }

        output_final(energy, iterations);
        timer.print_results(comm, root);
        counters.print_results(comm, root);
        if (xx != 0.0)
            local_wf.set_global_scalars(xx, xz, static_cast<ScaleFactorType>(scale));
        return local_wf;
    }

    Result get_result() const
    {
        Result result;
        result.energy = ground_state_energy;
        result.iterations = iterations;
        result.x_size_history = x_size_history;
        result.z_size_history = z_size_history;
        result.time_history = time_history;
        return result;
    }

private:
    void configure_openmp_runtime()
    {
#ifdef _OPENMP
        omp_set_dynamic(0);
        int slurm_cpus_per_task = cdfci_mpi::env_int("SLURM_CPUS_PER_TASK");
        if (slurm_cpus_per_task > 0)
            omp_set_num_threads(slurm_cpus_per_task);
        omp_threads = std::max(1, omp_get_max_threads());
#else
        omp_threads = 1;
#endif
    }

    void refresh_parallel_layout()
    {
        rank_omp_threads.assign(static_cast<size_t>(mpi_size), 1);
        int local_threads = std::max(1, omp_threads);
        cdfci_mpi::check(MPI_Allgather(&local_threads,
                                       1,
                                       MPI_INT,
                                       rank_omp_threads.data(),
                                       1,
                                       MPI_INT,
                                       comm),
                          "MPI_Allgather OpenMP thread layout");

        rank_thread_offsets.assign(static_cast<size_t>(mpi_size), 0);
        total_omp_workers = 0;
        for (int r = 0; r < mpi_size; ++r)
        {
            rank_thread_offsets[static_cast<size_t>(r)] = total_omp_workers;
            total_omp_workers += std::max(1, rank_omp_threads[static_cast<size_t>(r)]);
        }
        if (total_omp_workers < 1)
            total_omp_workers = mpi_size;
    }

    bool can_parallel_wavefunction_updates() const
    {
#if defined(_OPENMP) && !defined(CDFCI_SOLVER_SERIAL)
        return omp_threads > 1;
#else
        return false;
#endif
    }

    int routed_column_record_bytes() const
    {
        return cdfci_mpi::determinant_bytes<det_type>() + static_cast<int>(sizeof(int))
             + static_cast<int>(sizeof(double));
    }

    void append_routed_column_record(std::vector<unsigned char> &buffer,
                                     const det_type &target,
                                     int selected_index,
                                     double hij) const
    {
        const int bytes_per_det = cdfci_mpi::determinant_bytes<det_type>();
        const int record_bytes = routed_column_record_bytes();
        const size_t offset = buffer.size();
        buffer.resize(offset + static_cast<size_t>(record_bytes));

        unsigned char *dst = buffer.data() + offset;
        cdfci_mpi::pack_determinant(target, dst);
        std::memcpy(dst + bytes_per_det, &selected_index, sizeof(int));
        std::memcpy(dst + bytes_per_det + sizeof(int), &hij, sizeof(double));
    }

    void append_column_records(const typename H::Column &column,
                               int selected_index,
                               std::vector<std::vector<unsigned char>> &send_buffers) const
    {
        for (const auto &entry : column)
        {
            const det_type &target = entry.first;
            double hij = entry.second;
            int target_owner = cdfci_mpi::owner(target, mpi_size);
            append_routed_column_record(send_buffers[static_cast<size_t>(target_owner)],
                                        target,
                                        selected_index,
                                        hij);
        }
    }

    void build_routed_column_buffers(H &h,
                                     const std::vector<det_type> &selected,
                                     std::vector<std::vector<unsigned char>> &send_buffers,
                                     std::vector<std::vector<std::vector<unsigned char>>> &thread_send_buffers) const
    {
        if (send_buffers.size() != static_cast<size_t>(mpi_size))
            send_buffers.resize(static_cast<size_t>(mpi_size));
        for (auto &buffer : send_buffers)
            buffer.clear();

        // Rotate ownership across coordinates while keeping each pair unique.
        const size_t worker_count = static_cast<size_t>(total_omp_workers);
        const size_t pair_count = static_cast<size_t>(h.nelec) *
                                  static_cast<size_t>(h.nelec - 1) / 2;
        auto pair_worker = [&](int global_worker, size_t coordinate) {
            const size_t offset = ((coordinate % worker_count) *
                                   (pair_count % worker_count)) % worker_count;
            return static_cast<int>((static_cast<size_t>(global_worker) + worker_count - offset) % worker_count);
        };

        auto build_serial = [&]() {
            for (size_t p = 0; p < selected.size(); ++p)
            {
                det_type source = selected[p];
                int selected_index = static_cast<int>(p);

                // Diagonal and single-excitation terms are much smaller than the
                // double-excitation part, so assign them to one deterministic rank.
                if (mpi_rank == static_cast<int>(p % static_cast<size_t>(mpi_size)))
                {
                    typename H::Column serial_column = h.get_column_serial_part(source);
                    append_column_records(serial_column, selected_index, send_buffers);
                }

                const int pair_id = pair_worker(rank_thread_offsets[static_cast<size_t>(mpi_rank)], p);
                if (static_cast<size_t>(pair_id) < pair_count) {
                    typename H::Column parallel_column = h.get_column_parallel_part(source, pair_id, total_omp_workers);
                    append_column_records(parallel_column, selected_index, send_buffers);
                }
            }
        };

#ifdef _OPENMP
        if (omp_threads <= 1)
        {
            build_serial();
            return;
        }

        thread_send_buffers.resize(static_cast<size_t>(omp_threads));
        for (auto &thread_buffers : thread_send_buffers)
        {
            thread_buffers.resize(static_cast<size_t>(mpi_size));
            for (auto &buffer : thread_buffers)
                buffer.clear();
        }

#pragma omp parallel for num_threads(omp_threads) schedule(static)
        for (int tid = 0; tid < omp_threads; ++tid)
        {
            const int global_tid =
                rank_thread_offsets[static_cast<size_t>(mpi_rank)] + tid;
            auto &local_send_buffers = thread_send_buffers[static_cast<size_t>(tid)];

            for (size_t p = 0; p < selected.size(); ++p)
            {
                det_type source = selected[p];
                int selected_index = static_cast<int>(p);

                // Keep the small diagonal/single part deterministic and unique.
                if (tid == 0
                    && mpi_rank == static_cast<int>(p % static_cast<size_t>(mpi_size)))
                {
                    typename H::Column serial_column = h.get_column_serial_part(source);
                    append_column_records(serial_column,
                                          selected_index,
                                          local_send_buffers);
                }

                const int pair_id = pair_worker(global_tid, p);
                if (static_cast<size_t>(pair_id) < pair_count) {
                    typename H::Column parallel_column = h.get_column_parallel_part(source, pair_id, total_omp_workers);
                    append_column_records(parallel_column, selected_index, local_send_buffers);
                }
            }
        }

        // For threaded execution, keep thread-local buffers as the transport
        // source. exchange_threaded() performs the only required copy into
        // the final rank-contiguous packed buffer.
#else
        for (size_t p = 0; p < selected.size(); ++p)
        {
            det_type source = selected[p];
            int selected_index = static_cast<int>(p);

            // Diagonal and single-excitation terms are much smaller than the
            // double-excitation part, so assign them to one deterministic rank.
            if (mpi_rank == static_cast<int>(p % static_cast<size_t>(mpi_size)))
            {
                typename H::Column serial_column = h.get_column_serial_part(source);
                append_column_records(serial_column, selected_index, send_buffers);
            }

            const int pair_id = pair_worker(mpi_rank, p);
            if (static_cast<size_t>(pair_id) < pair_count) {
                typename H::Column parallel_column = h.get_column_parallel_part(source, pair_id, mpi_size);
                append_column_records(parallel_column, selected_index, send_buffers);
            }
        }
#endif
    }

    template <typename PendingMap>
    void consume_routed_column_buffer(const std::vector<unsigned char> &recv_buffer,
                                      const std::vector<double> &delta_x_I,
                                      const W &local_wf,
                                      PendingMap &pending_dz,
                                      std::vector<double> &local_exact_z_I) const
    {
        const int bytes_per_det = cdfci_mpi::determinant_bytes<det_type>();
        const int record_bytes = routed_column_record_bytes();
        if (recv_buffer.empty())
            return;
        if (recv_buffer.size() % static_cast<size_t>(record_bytes) != 0)
            throw std::runtime_error("Corrupted MPI routed column buffer size.");

        const size_t num_records = recv_buffer.size() / static_cast<size_t>(record_bytes);

#ifdef _OPENMP
        if (omp_threads > 1)
        {
            std::vector<PendingMap> thread_pending(static_cast<size_t>(omp_threads));
            std::vector<std::vector<double>> thread_exact_z(
                static_cast<size_t>(omp_threads),
                std::vector<double>(local_exact_z_I.size(), 0.0));

#pragma omp parallel num_threads(omp_threads)
            {
                const int tid = omp_get_thread_num();
                const int team_size = omp_get_num_threads();
                const size_t begin_record =
                    num_records * static_cast<size_t>(tid) / static_cast<size_t>(team_size);
                const size_t end_record =
                    num_records * static_cast<size_t>(tid + 1) / static_cast<size_t>(team_size);
                auto &local_pending = thread_pending[static_cast<size_t>(tid)];
                auto &local_exact_z = thread_exact_z[static_cast<size_t>(tid)];
                local_pending.reserve(end_record - begin_record);

                for (size_t record = begin_record; record < end_record; ++record)
                {
                    const unsigned char *src =
                        recv_buffer.data() + record * static_cast<size_t>(record_bytes);
                    det_type target = cdfci_mpi::unpack_determinant<det_type>(src);
                    int selected_index = 0;
                    double hij = 0.0;
                    std::memcpy(&selected_index, src + bytes_per_det, sizeof(int));
                    std::memcpy(&hij, src + bytes_per_det + sizeof(int), sizeof(double));

                    if (selected_index < 0
                        || static_cast<size_t>(selected_index) >= delta_x_I.size()
                        || static_cast<size_t>(selected_index) >= local_exact_z_I.size())
                    {
                        throw std::runtime_error("Corrupted MPI routed column selected index.");
                    }

                    const size_t p = static_cast<size_t>(selected_index);
                    local_pending[target] += hij * delta_x_I[p];
                    mapped_type value = local_wf.local_get_value_or_zero(target);
                    local_exact_z[p] += hij * value[0];
                }
            }

            pending_dz.reserve(pending_dz.size() + num_records);
            for (int t = 0; t < omp_threads; ++t)
            {
                for (auto &entry : thread_pending[static_cast<size_t>(t)])
                    pending_dz[entry.first] += entry.second;
                for (size_t p = 0; p < local_exact_z_I.size(); ++p)
                    local_exact_z_I[p] += thread_exact_z[static_cast<size_t>(t)][p];
            }
            return;
        }
#endif

        pending_dz.reserve(pending_dz.size() + num_records);
        for (size_t offset = 0; offset < recv_buffer.size(); offset += static_cast<size_t>(record_bytes))
        {
            const unsigned char *src = recv_buffer.data() + offset;
            det_type target = cdfci_mpi::unpack_determinant<det_type>(src);
            int selected_index = 0;
            double hij = 0.0;
            std::memcpy(&selected_index, src + bytes_per_det, sizeof(int));
            std::memcpy(&hij, src + bytes_per_det + sizeof(int), sizeof(double));

            if (selected_index < 0
                || static_cast<size_t>(selected_index) >= delta_x_I.size()
                || static_cast<size_t>(selected_index) >= local_exact_z_I.size())
            {
                throw std::runtime_error("Corrupted MPI routed column selected index.");
            }

            pending_dz[target] += hij * delta_x_I[static_cast<size_t>(selected_index)];
            mapped_type value = local_wf.local_get_value_or_zero(target);
            local_exact_z_I[static_cast<size_t>(selected_index)] += hij * value[0];
        }
    }

    template <typename CandidateSet>
    void consume_routed_column_buffer_direct(const std::vector<unsigned char> &recv_buffer,
                                             const std::vector<double> &delta_x_I,
                                             double new_scale,
                                             W &local_wf,
                                             std::vector<double> &local_exact_z_I,
                                             CandidateSet &next_seen) const
    {
        const int bytes_per_det = cdfci_mpi::determinant_bytes<det_type>();
        const int record_bytes = routed_column_record_bytes();
        if (recv_buffer.empty())
            return;
        if (recv_buffer.size() % static_cast<size_t>(record_bytes) != 0)
            throw std::runtime_error("Corrupted MPI routed column buffer size.");

#if defined(_OPENMP) && !defined(CDFCI_SOLVER_SERIAL)
        if (can_parallel_wavefunction_updates())
        {
            const size_t num_records =
                recv_buffer.size() / static_cast<size_t>(record_bytes);
            std::vector<std::vector<double>> thread_exact_z(
                static_cast<size_t>(omp_threads),
                std::vector<double>(local_exact_z_I.size(), 0.0));
            std::vector<CandidateSet> thread_seen(static_cast<size_t>(omp_threads));

#pragma omp parallel num_threads(omp_threads)
            {
                const int tid = omp_get_thread_num();
                const int team_size = omp_get_num_threads();
                const size_t begin_record =
                    num_records * static_cast<size_t>(tid) / static_cast<size_t>(team_size);
                const size_t end_record =
                    num_records * static_cast<size_t>(tid + 1) / static_cast<size_t>(team_size);
                auto &local_exact_z = thread_exact_z[static_cast<size_t>(tid)];
                auto &local_seen = thread_seen[static_cast<size_t>(tid)];
                local_seen.reserve(end_record - begin_record);

                for (size_t record = begin_record; record < end_record; ++record)
                {
                    const unsigned char *src =
                        recv_buffer.data() + record * static_cast<size_t>(record_bytes);
                    det_type target = cdfci_mpi::unpack_determinant<det_type>(src);
                    int selected_index = 0;
                    double hij = 0.0;
                    std::memcpy(&selected_index, src + bytes_per_det, sizeof(int));
                    std::memcpy(&hij, src + bytes_per_det + sizeof(int), sizeof(double));

                    if (selected_index < 0
                        || static_cast<size_t>(selected_index) >= delta_x_I.size()
                        || static_cast<size_t>(selected_index) >= local_exact_z_I.size())
                    {
                        throw std::runtime_error("Corrupted MPI routed column selected index.");
                    }

                    const size_t p = static_cast<size_t>(selected_index);
                    double dz = hij * delta_x_I[p];
                    if (std::fabs(new_scale * dz) > 0.0)
                    {
                        double x = local_wf.local_add_z_and_get_x(target, dz);
                        local_exact_z[p] += hij * x;
                        local_seen.insert(target);
                        continue;
                    }
                    mapped_type value = local_wf.local_get_value_or_zero(target);
                    local_exact_z[p] += hij * value[0];

                    bool keep = value[0] != 0.0
                             || value[1] != 0.0
                             || std::fabs(new_scale * dz) > 0.0;
                    if (keep)
                    {
                        local_wf.local_add_z(target, dz);
                        local_seen.insert(target);
                    }
                }
            }

            for (int t = 0; t < omp_threads; ++t)
            {
                for (size_t p = 0; p < local_exact_z_I.size(); ++p)
                    local_exact_z_I[p] += thread_exact_z[static_cast<size_t>(t)][p];
                for (const auto &det : thread_seen[static_cast<size_t>(t)])
                    next_seen.insert(det);
            }
            return;
        }
#endif

        for (size_t offset = 0; offset < recv_buffer.size(); offset += static_cast<size_t>(record_bytes))
        {
            const unsigned char *src = recv_buffer.data() + offset;
            det_type target = cdfci_mpi::unpack_determinant<det_type>(src);
            int selected_index = 0;
            double hij = 0.0;
            std::memcpy(&selected_index, src + bytes_per_det, sizeof(int));
            std::memcpy(&hij, src + bytes_per_det + sizeof(int), sizeof(double));

            if (selected_index < 0
                || static_cast<size_t>(selected_index) >= delta_x_I.size()
                || static_cast<size_t>(selected_index) >= local_exact_z_I.size())
            {
                throw std::runtime_error("Corrupted MPI routed column selected index.");
            }

            double dz = hij * delta_x_I[static_cast<size_t>(selected_index)];
            if (std::fabs(new_scale * dz) > 0.0)
            {
                double x = local_wf.local_add_z_and_get_x(target, dz);
                local_exact_z_I[static_cast<size_t>(selected_index)] += hij * x;
                next_seen.insert(target);
                continue;
            }
            mapped_type value = local_wf.local_get_value_or_zero(target);
            local_exact_z_I[static_cast<size_t>(selected_index)] += hij * value[0];

            bool keep = value[0] != 0.0
                     || value[1] != 0.0
                     || std::fabs(new_scale * dz) > 0.0;
            if (keep)
            {
                local_wf.local_add_z(target, dz);
                next_seen.insert(target);
            }
        }
    }

    det_type get_ref_det(H &h) const
    {
        det_type hf;
        if (ref_det_occ.empty())
            return h.get_hartree_fock();

        for (auto orb : ref_det_occ)
        {
            if ((orb >= 0) && (orb < h.norb))
                hf.set_orbital(orb);
        }
        if (hf.get_occupied_orbitals().size() != static_cast<size_t>(h.nelec))
        {
            throw std::invalid_argument(
                "Reference determinant is invalid. The number of valid orbitals does "
                "not match the number of electrons.");
        }
        return hf;
    }

    std::vector<cdfci_mpi::Candidate<det_type>>
    select_local_top(const W &local_wf,
                     const std::vector<det_type> &local_candidates,
                     double xx,
                     double scale) const
    {
        auto better_candidate = [](const auto &lhs, const auto &rhs) {
            if (lhs.score != rhs.score)
                return lhs.score > rhs.score;
            // Keep the tie rule identical in local pruning and global selection.
            return cdfci_mpi::determinant_less(rhs.det, lhs.det);
        };
        auto trim_scored = [&](auto &scored) {
            if (scored.size() > num_coordinates)
            {
                std::nth_element(scored.begin(),
                                 scored.begin() + num_coordinates,
                                 scored.end(),
                                 better_candidate);
                scored.resize(num_coordinates);
            }
        };

        const double scaled_xx = scale * scale * xx;

#ifdef _OPENMP
        if (omp_threads > 1)
        {
            std::vector<std::vector<cdfci_mpi::Candidate<det_type>>> thread_scored(
                static_cast<size_t>(omp_threads));
            const long long candidate_count =
                static_cast<long long>(local_candidates.size());

#pragma omp parallel num_threads(omp_threads)
            {
                const int tid = omp_get_thread_num();
                auto &scored = thread_scored[static_cast<size_t>(tid)];
                using Candidate = cdfci_mpi::Candidate<det_type>;
                // With "better" as the priority_queue comparator, top is the worst retained item.
                std::vector<Candidate> storage;
                storage.reserve(num_coordinates + 1);
                std::priority_queue<Candidate, std::vector<Candidate>, decltype(better_candidate)>
                    heap(better_candidate, std::move(storage));

#pragma omp for schedule(static)
                for (long long idx = 0; idx < candidate_count; ++idx)
                {
                    const auto &det = local_candidates[static_cast<size_t>(idx)];
                    if (cdfci_mpi::owner(det, mpi_size) != mpi_rank)
                        continue;
                    mapped_type value = local_wf.local_get_value_or_zero(det);
                    double score = std::fabs(value[1] + scaled_xx * value[0]);
                    Candidate item{det, score, true};
                    if (heap.size() < num_coordinates)
                        heap.push(item);
                    else if (num_coordinates && better_candidate(item, heap.top()))
                    {
                        heap.pop();
                        heap.push(item);
                    }
                }
                scored.reserve(heap.size());
                while (!heap.empty()) {
                    scored.push_back(heap.top());
                    heap.pop();
                }
            }

            std::vector<cdfci_mpi::Candidate<det_type>> scored;
            size_t total_scored = 0;
            for (const auto &local_scored : thread_scored)
                total_scored += local_scored.size();
            scored.reserve(total_scored);
            for (auto &local_scored : thread_scored)
                scored.insert(scored.end(), local_scored.begin(), local_scored.end());
            trim_scored(scored);
            std::sort(scored.begin(), scored.end(), better_candidate);
            return scored;
        }
#endif

        std::vector<cdfci_mpi::Candidate<det_type>> scored;
        scored.reserve(local_candidates.size());
        for (const auto &det : local_candidates)
        {
            if (cdfci_mpi::owner(det, mpi_size) != mpi_rank)
                continue;
            mapped_type value = local_wf.local_get_value_or_zero(det);
            double score = std::fabs(value[1] + scaled_xx * value[0]);
            scored.push_back({det, score, true});
        }

        if (scored.size() > num_coordinates)
            trim_scored(scored);
        std::sort(scored.begin(), scored.end(), better_candidate);
        return scored;
    }

    void get_selected_values(const W &local_wf,
                             const std::vector<det_type> &selected,
                             std::vector<double> &x_I,
                             std::vector<double> &z_I) const
    {
        std::vector<double> local(2 * selected.size(), 0.0);
        for (size_t p = 0; p < selected.size(); ++p)
        {
            if (cdfci_mpi::owner(selected[p], mpi_size) == mpi_rank)
            {
                mapped_type value = local_wf.local_get_value_or_zero(selected[p]);
                local[p] = value[0];
                local[p + selected.size()] = value[1];
            }
        }
        std::vector<double> global;
        cdfci_mpi::allreduce_sum(local, global, comm);
        x_I.assign(selected.size(), 0.0);
        z_I.assign(selected.size(), 0.0);
        for (size_t p = 0; p < selected.size(); ++p)
        {
            x_I[p] = global[p];
            z_I[p] = global[p + selected.size()];
        }
    }

    std::vector<double> build_hii(H &h, const std::vector<det_type> &selected) const
    {
        const size_t n = selected.size();
        std::vector<double> H_II(n * n, 0.0);
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j <= i; ++j)
            {
                det_type det_i = selected[i];
                det_type det_j = selected[j];
                double hij = h.get_entry(det_i, det_j);
                H_II[i * n + j] = hij;
                H_II[j * n + i] = hij;
            }
        }
        return H_II;
    }

    double solve_block_update(const std::vector<double> &H_II,
                              const std::vector<double> &x_I,
                              const std::vector<double> &z_I,
                              double xx,
                              double xz,
                              std::vector<double> &delta_x_I) const
    {
        if (coord_update_str == "ls")
        {
            double x = x_I[0];
            double z = z_I[0];
            double dA = H_II[0];
            double b = 3.0 * x;
            double c = dA + 2.0 * x * x + xx;
            double d = x * xx + z;
            delta_x_I[0] = solve_cubic(b, c, d, 1e-10, 10);
            return 1.0;
        }

        const size_t n = x_I.size();
        Eigen::VectorXd x_vec(n);
        Eigen::VectorXd z_vec(n);
        Eigen::MatrixXd H_mat(n, n);
        for (size_t i = 0; i < n; ++i)
        {
            x_vec(static_cast<int>(i)) = x_I[i];
            z_vec(static_cast<int>(i)) = z_I[i];
            for (size_t j = 0; j < n; ++j)
                H_mat(static_cast<int>(i), static_cast<int>(j)) = H_II[i * n + j];
        }

        Eigen::SelfAdjointView<Eigen::MatrixXd, Eigen::Lower> H_sym =
            H_mat.selfadjointView<Eigen::Lower>();

        auto xTx = x_vec.squaredNorm();
        auto xTz = x_vec.dot(z_vec);
        auto Hx_vec = H_sym * x_vec;
        auto xTHx = x_vec.dot(Hx_vec);
        double yy = xx - xTx;
        if (yy < 0.0 && std::fabs(yy) < 1e-12)
            yy = 0.0;

        double scale = 1.0;
        if (std::fabs(yy) > std::numeric_limits<double>::epsilon())
        {
            double y1 = std::sqrt(yy);
            Eigen::MatrixXd A_mat(n + 1, n + 1);
            A_mat.setZero();
            A_mat(0, 0) = (xz - 2.0 * xTz + xTHx) / yy;
            A_mat.block(1, 0, static_cast<int>(n), 1) = (z_vec - Hx_vec) / y1;
            A_mat.block(0, 1, 1, static_cast<int>(n)) =
                A_mat.block(1, 0, static_cast<int>(n), 1).transpose();
            A_mat.block(1, 1, static_cast<int>(n), static_cast<int>(n)) = H_mat;

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(A_mat);
            if (solver.info() != Eigen::Success)
                throw std::runtime_error("MPI block update eigensolver failed.");

            double eval = solver.eigenvalues()[0];
            Eigen::VectorXd evec = solver.eigenvectors().col(0);
            double lambda = -eval;
            double gamma = evec(0) / y1;
            scale = gamma * std::sqrt(lambda);
            for (size_t i = 0; i < n; ++i)
                delta_x_I[i] = evec(static_cast<int>(i + 1)) / gamma - x_vec(static_cast<int>(i));
        }
        else
        {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H_sym);
            if (solver.info() != Eigen::Success)
                throw std::runtime_error("MPI block update eigensolver failed.");

            double eval = solver.eigenvalues()[0];
            Eigen::VectorXd evec = solver.eigenvectors().col(0);
            scale = 1.0;
            for (size_t i = 0; i < n; ++i)
                delta_x_I[i] = evec(static_cast<int>(i)) * std::sqrt(-eval)
                             - x_vec(static_cast<int>(i));
        }
        return scale;
    }

    double solve_cubic(double b,
                       double c,
                       double d,
                       const double depsilon,
                       const int max_iter) const
    {
        double p = c - b * b / 3.0;
        double q = d + b * (2.0 / 27.0 * b * b - c / 3.0);
        double shift = -b / 3.0;
        double p3 = p / 3.0;
        double q2 = q / 2.0;
        double delta = q2 * q2 + p3 * p3 * p3;
        double root = 0.0;
        const double pi = std::acos(-1.0);

        if (delta >= 0.0)
        {
            double sqrt_delta = std::sqrt(delta);
            root = std::cbrt(-q2 + sqrt_delta) + std::cbrt(-q2 - sqrt_delta);
        }
        else
        {
            double r = std::sqrt(-p3);
            double qrtd = std::sqrt(-delta);
            double theta;
            if (q2 >= 0.0)
                theta = (std::atan2(-qrtd, -q2) - 2.0 * pi) / 3.0;
            else
                theta = std::atan2(qrtd, -q2) / 3.0;
            root = 2.0 * r * std::cos(theta);
        }

        double x = root + shift;
        auto newton_one_step = [b, c, d](double z) {
            double f = z * (z * (z + b) + c) + d;
            double df = z * (3.0 * z + 2.0 * b) + c;
            return z - f / df;
        };

        int iter = 0;
        double xn = newton_one_step(x);
        while (x != 0.0 && std::fabs((xn - x) / x) > depsilon && iter < max_iter)
        {
            x = xn;
            xn = newton_one_step(x);
            ++iter;
        }
        return xn;
    }

    double dot(const std::vector<double> &lhs, const std::vector<double> &rhs) const
    {
        double result = 0.0;
        for (size_t i = 0; i < lhs.size(); ++i)
            result += lhs[i] * rhs[i];
        return result;
    }

    double quadratic_form(const std::vector<double> &d,
                          const std::vector<double> &matrix) const
    {
        double result = 0.0;
        const size_t n = d.size();
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                result += d[i] * matrix[i * n + j] * d[j];
        return result;
    }

    bool check_convergence(const std::vector<double> &delta_x_I,
                           double scale,
                           double &dx_accumulated) const
    {
        double dx_norm = 0.0;
        for (double dx : delta_x_I)
            dx_norm += dx * dx;
        dx_norm = std::sqrt(dx_norm) * std::fabs(scale);
        dx_accumulated =
            stopping_dx_damping_factor * dx_accumulated + dx_norm;
        return dx_accumulated < stopping_dx_threshold;
    }

    void output_header() const
    {
        if (verbose == 0 || !cdfci_mpi::is_io_root(comm, root))
            return;

        // |K| is the global number of candidate determinants kept for the next selection.
        std::cout << std::setw(13) << std::left << "Iteration";
        std::cout << std::setw(18) << std::right << "Energy";
        std::cout << std::setw(18) << "dx";
        std::cout << std::setw(15) << "|x|_0";
        std::cout << std::setw(15) << "|z|_0";
        std::cout << std::setw(9) << "|K|";
        std::cout << std::setw(10) << "Time";
        std::cout << std::endl;
    }

    void output_line(size_t iteration,
                     double energy,
                     double dx,
                     size_t x_size,
                     size_t z_size,
                     size_t candidate_size,
                     double time)
    {
        x_size_history.push_back(x_size);
        z_size_history.push_back(z_size);
        time_history.push_back(time);

        if (verbose == 0 || !cdfci_mpi::is_io_root(comm, root))
            return;

        std::cout << std::setw(13) << std::left << iteration;
        std::cout << std::setw(18) << std::right << std::fixed
                  << std::setprecision(10) << energy;
        std::cout << std::setw(18) << std::scientific
                  << std::setprecision(4) << dx;
        std::cout << std::setw(15) << x_size;
        std::cout << std::setw(15) << z_size;
        std::cout << std::setw(9) << candidate_size;
        std::cout << std::setw(10) << std::fixed
                  << std::setprecision(2) << time;
        std::cout << std::endl;
    }

    void output_final(double energy, size_t iterations_)
    {
        ground_state_energy = energy;
        iterations = iterations_;
        if (verbose == 0 || !cdfci_mpi::is_io_root(comm, root))
            return;

        std::cout << "Final FCI Energy: " << std::setw(30)
                  << std::right << std::fixed
                  << std::setprecision(16) << ground_state_energy << std::endl;
        std::cout << "Iterations: " << iterations << std::endl;
    }
};

#endif
