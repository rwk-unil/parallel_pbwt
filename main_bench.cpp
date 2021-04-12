#include<benchmark/benchmark.h>
#include<pbwt_exp.hpp>

constexpr size_t MAX_THREADS = 4;

///////////////////////////
// Benchmark definitions //
///////////////////////////

// Benchmark how long it takes to generate the THREADS-1 required starting a and d arrays
static void BM_generate_a_d_arrays_sequentially(benchmark::State& state) {
    const size_t THREADS = state.range(0); // Number of positions required (to start THREADS threads)
    auto hap_map = read_from_macs_file<bool>("11k.macs"); // TODO CHANGE FILE
    const size_t N = hap_map.size(); // Number of variant sites
    auto positions_to_collect = generate_positions_to_collect(N, THREADS);

    std::vector<a_d_arrays_at_pos> result;

    // BENCHMARKED CODE
    for (auto _ : state) {
        result = generate_a_d_arrays_for_positions_sequentially(hap_map, positions_to_collect);
    }
    // END BENCHMARKED CODE
}

// Benchmark how long it takes to generate the THREADS-1 required starting a and d arrays with THREADS-1 threads
static void BM_generate_a_d_arrays_in_parallel(benchmark::State& state) {
    const size_t THREADS = state.range(0); // Number of positions required (to start THREADS threads)
    auto hap_map = read_from_macs_file<bool>("11k.macs");
    const size_t N = hap_map.size(); // Number of variant sites
    auto positions_to_collect = generate_positions_to_collect(N, THREADS);

    std::vector<a_d_arrays_at_pos> result;

    // BENCHMARKED CODE
    for (auto _ : state) {
        result = generate_a_d_arrays_for_positions_in_parallel(hap_map, positions_to_collect);
    }
    // END BENCHMARKED CODE
}

// Benchmark how long it takes to report set maximal matches
static void BM_report_matches_sequential(benchmark::State& state) {
    auto hap_map = read_from_macs_file<bool>("11k.macs");

    matches_t result;

    // BENCHMARKED CODE
    for (auto _ : state) {
        result = report_matches_sequentially(hap_map);
    }
    // END BENCHMARKED CODE
}

static void BM_report_matches_in_parallel_a_d_sequentially_generated(benchmark::State& state) {
    auto hap_map = read_from_macs_file<bool>("11k.macs");

    std::vector<matches_t> result;

    // BENCHMARKED CODE
    for (auto _ : state) {
        result = report_matches_in_parallel_a_d_sequential(hap_map, state.range(0));
    }
    // END BENCHMARKED CODE
}

static void BM_report_matches_in_parallel(benchmark::State& state) {
    auto hap_map = read_from_macs_file<bool>("11k.macs");

    std::vector<matches_t> result;

    // BENCHMARKED CODE
    for (auto _ : state) {
        result = report_matches_in_parallel(hap_map, state.range(0));
    }
    // END BENCHMARKED CODE
}

////////////////
// Benchmarks //
////////////////

// The ranges are the number of threads that would be used to process the matrix in parallel
// this value-1 is the number of starting positions generated (because first position is natural order)

// Generate a and d sequentially (single thread)
BENCHMARK(BM_generate_a_d_arrays_sequentially)->Range(2, MAX_THREADS)->MeasureProcessCPUTime()->UseRealTime();
// Generate a and d in parallel with sequential fix (range-1 threads + sequential fix)
BENCHMARK(BM_generate_a_d_arrays_in_parallel)->Range(2, MAX_THREADS)->MeasureProcessCPUTime()->UseRealTime();

// Reference
BENCHMARK(BM_report_matches_sequential);
// First strategy
BENCHMARK(BM_report_matches_in_parallel_a_d_sequentially_generated)->Range(2, MAX_THREADS)->MeasureProcessCPUTime()->UseRealTime();
// Second strategy
BENCHMARK(BM_report_matches_in_parallel)->Range(2, MAX_THREADS)->MeasureProcessCPUTime()->UseRealTime();

#if 0
// This will change how much the parallelisation improves the runtime
constexpr size_t SUBSAMPLING_RATE = 800;
constexpr size_t THREADS = 4;

static void BM_alg2(benchmark::State& state) {
    auto hap_map = read_from_macs_file<bool>("11k.macs");

    const size_t ss_rate = SUBSAMPLING_RATE;

    std::vector<alg2_res_t> result;

    for (auto _ : state) {
        result = algorithm_2<false /* Report Set Maximal Matches */>(hap_map, ss_rate);
    }
    //print_vector(result.back().a);
    //print_vector(result.back().d);
}

static void BM_alg2_4(benchmark::State& state) {
    auto hap_map = read_from_macs_file<bool>("11k.macs");

    const size_t ss_rate = SUBSAMPLING_RATE;
    matches_t matches;

    for (auto _ : state) {
        std::vector<alg2_res_t> result;
        result = algorithm_2<true /* Report Set Maximal Matches */>(hap_map, ss_rate, matches);
        //print_vector(result.back().a);
    }

    std::string filename = "normal_matches.txt";
    std::fstream s(filename, s.out);
    if (s.is_open()) {
        for (const auto& m : matches) {
            //std::cout << "Match between " << m.a << " and " << m.b << " from [" <<
            //             m.start << " to " << m.end << "[" << std::endl;

            // Same format as Durbin
            s << "MATCH\t" << m.a << "\t" << m.b << "\t" << m.start << "\t" << m.end << "\t" << m.end-m.start << std::endl;
        }
        s.close();
    } else {
        // Failed to open file for some reason
    }
}

#include <thread>
#include <map>
static void BM_alg2_exp(benchmark::State& state) {
    auto hap_map = read_from_macs_file<bool>("11k.macs");
    const size_t M = hap_map.size();
    const size_t ss_rate = SUBSAMPLING_RATE;
    //const size_t look_back = 0; // DO NOT GO BACK / LOOK BACK, this is crap
    const size_t jumps = M / ss_rate;
    const size_t jumps_per_thread = jumps / THREADS;
    const size_t last_jumps = jumps - jumps_per_thread * THREADS;

    const size_t jump = ss_rate;
    const size_t last_extra_len = M - jump * jumps;

    std::vector<alg2_res_t> results(jumps);

    for (auto _ : state) {
        std::vector<std::thread> workers(THREADS);
        for (size_t i = 0; i < THREADS; ++i) {
            workers[i] = std::thread([=, &hap_map, &results]{
                const size_t jumps_to_do = jumps_per_thread + (i == THREADS-1 ? last_jumps : 0);
                for (size_t j = 0; j < jumps_to_do; ++j) {
                    //const size_t go_back = (i and j) ? look_back : 0;
                    const size_t offset = i*jumps_per_thread+j;
                    const size_t len = jump + /*go_back +*/ ((i == THREADS-1) and (j == jumps_to_do-1) ? last_extra_len : 0);
                    results[offset] = algorithm_2_exp(hap_map, offset*jump/*-go_back*/, len);
                }
                //results[i] = algorithm_2_exp(hap_map, i*jump, jump + ((i == (THREADS-1)) ? last_extra_len : 0));
            });
        }
        for (auto& w : workers) {
            w.join();
        }
        fix_a_d(results);
        //fix_a_d<true>(results); // This shows the number of fixes
        //print_vector(results.back().a);
        //benchmark::DoNotOptimize(encode_all_a(results)); // This is just to measure the cost of this extra operation // Cost is almost nil
    }

    //print_vector(results.back().a);
    //print_vector(results.back().d);

    auto encoded_as = encode_all_a(results);
    // Do some statistics on the encoded a's
    size_t counter = 0;
    for (const auto& e : encoded_as) {
        counter += e.size();
    }

    std::cout << "Normal  a's is " << results.size() * results[0].a.size() << " integer values" << std::endl;
    std::cout << "Encoded a's is " << counter*2 << " integer values" << std::endl;

    // Check what the d vectors look like
    std::map<size_t, size_t> d_map;
    for (const auto& e : results) {
        //print_vector(e.d);
        for (const auto d : e.d) {
            if (d_map.find(d) != d_map.end()) {
                d_map[d]++;
            } else {
                d_map[d] = 1;
            }
        }
    }

    std::cout << "Map size : " << d_map.size() << std::endl;
    // for (const auto& e : d_map) {
    //     std::cout << "{" << e.first << "," << e.second << "}" << std::endl; // This will massively spam
    // }
}

static void BM_alg2_4_exp(benchmark::State& state) {
    auto hap_map = read_from_macs_file<bool>("11k.macs");
    const size_t M = hap_map.size();
    const size_t ss_rate = SUBSAMPLING_RATE;

    const size_t jumps = M / ss_rate;
    const size_t jumps_per_thread = jumps / THREADS;
    const size_t last_jumps = jumps - jumps_per_thread * THREADS;

    const size_t jump = ss_rate;
    const size_t last_extra_len = M - jump * jumps;

    // These two structures are shared between threads but not modified, only contents are touched, therefore no synchronization is needed
    std::vector<alg2_res_t> results(jumps);
    std::vector<matches_t> matches(jumps);

    for (auto _ : state) {
        // These are because the benchmark is run multiple times and the data structures persist
        matches.clear();
        matches.insert(matches.begin(), jumps, {});

        // Compute the estimated a,d's in parallel
        std::vector<std::thread> workers(THREADS);
        for (size_t i = 0; i < THREADS; ++i) {
            workers[i] = std::thread([=, &hap_map, &results]{
                const size_t jumps_to_do = jumps_per_thread + (i == THREADS-1 ? last_jumps : 0);
                for (size_t j = 0; j < jumps_to_do; ++j) {
                    const size_t offset = i*jumps_per_thread+j;
                    const size_t len = jump + ((i == THREADS-1) and (j == jumps_to_do-1) ? last_extra_len : 0);
                    results[offset] = algorithm_2_exp<false /* Rep Matches */>(hap_map, offset*jump, len);
                }
            });
        }
        for (auto& w : workers) {
            w.join();
        }
        // Sequentially fix the a,d's
        fix_a_d(results);

        // Use the fixed a,d's for parallel matching
        for (size_t i = 0; i < THREADS; ++i) {
            workers[i] = std::thread([=, &hap_map, &results, &matches]{
                const size_t jumps_to_do = jumps_per_thread + (i == THREADS-1 ? last_jumps : 0);
                for (size_t j = 0; j < jumps_to_do; ++j) {
                    const size_t offset = i*jumps_per_thread+j;
                    const size_t len = jump + ((i == THREADS-1) and (j == jumps_to_do-1) ? last_extra_len : 0);
                    results[offset] = algorithm_2_exp<true /* Rep Matches */>(hap_map, offset*jump, len, offset ? results[offset-1].a : ppa_t(0), offset ? results[offset-1].d : d_t(0), matches[offset]);
                }
            });
        }
        for (auto& w : workers) {
            w.join();
        }
    }

    std::string filename = "parallel_matches.txt";
    std::fstream s(filename, s.out);
    if (!s.is_open()) return;

    size_t match_counter = 0;
    size_t block_counter = 0;
    for (const auto& m_vect : matches) {
        match_counter += m_vect.size();
        s << "---- Block : " << block_counter++ << std::endl;
        for (const auto& m : m_vect) {
            s << "MATCH\t" << m.a << "\t" << m.b << "\t" << m.start << "\t" << m.end << "\t" << m.end-m.start << std::endl;
        }
    }
    s.close();
    std::cout << "Found " << match_counter << " matches" << std::endl;
    std::cout << "Expected matches should be 373344" << std::endl; // wc -l on file
}

static void BM_alg2_4_exp_thread_optimal(benchmark::State& state) {
    auto hap_map = read_from_macs_file<bool>("11k.macs");
    const size_t M = hap_map.size();

    const size_t ss_rate = M / THREADS;
    const size_t last_extra_len = M - ss_rate * THREADS;

    // These two structures are shared between threads but not modified, only contents are touched, therefore no synchronization is needed
    std::vector<alg2_res_t> results(THREADS);
    std::vector<matches_t> matches(THREADS);

    for (auto _ : state) {
        // These are because the benchmark is run multiple times and the data structures persist
        matches.clear();
        matches.insert(matches.begin(), THREADS, {});

        // Compute the estimated a,d's in parallel
        std::vector<std::thread> workers(THREADS);
        for (size_t i = 0; i < THREADS; ++i) {
            workers[i] = std::thread([=, &hap_map, &results]{
                const size_t len = ss_rate + ((i == THREADS-1) ? last_extra_len : 0);
                results[i] = algorithm_2_exp<false /* Rep Matches */>(hap_map, i*ss_rate, len);
            });
        }
        for (auto& w : workers) {
            w.join();
        }
        // Sequentially fix the a,d's
        fix_a_d(results);

        /// @todo THIS CAN BE IMPROVED BY :
        /// Matching for thread 0 above and only running the code below for the remaining threads, however, the threads above may have to wait on the first one to finish so it could be counter productive

        // Use the fixed a,d's for parallel matching
        for (size_t i = 0; i < THREADS; ++i) {
            workers[i] = std::thread([=, &hap_map, &results, &matches]{
                const size_t len = ss_rate + ((i == THREADS-1) ? last_extra_len : 0);
                results[i] = algorithm_2_exp<true /* Rep Matches */>(hap_map, i*ss_rate, len, i ? results[i-1].a : ppa_t(0), i ? results[i-1].d : d_t(0), matches[i]);
            });
        }
        for (auto& w : workers) {
            w.join();
        }
    }

    std::string filename = "parallel_matches_thread_optimal.txt";
    std::fstream s(filename, s.out);
    if (!s.is_open()) return;

    size_t match_counter = 0;
    size_t block_counter = 0;
    for (const auto& m_vect : matches) {
        match_counter += m_vect.size();
        s << "---- Block : " << block_counter++ << std::endl;
        for (const auto& m : m_vect) {
            s << "MATCH\t" << m.a << "\t" << m.b << "\t" << m.start << "\t" << m.end << "\t" << m.end-m.start << std::endl;
        }
    }
    s.close();
    std::cout << "Found " << match_counter << " matches" << std::endl;
    std::cout << "Expected matches should be 373344" << std::endl; // wc -l on file
}

BENCHMARK(BM_alg2);
BENCHMARK(BM_alg2_exp);
BENCHMARK(BM_alg2_4);
BENCHMARK(BM_alg2_4_exp);
BENCHMARK(BM_alg2_4_exp_thread_optimal);
#endif /* 0 */
BENCHMARK_MAIN();