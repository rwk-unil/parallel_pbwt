#include <cstdint>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <thread>
#include <mutex>
#include <algorithm>

#include "synced_bcf_reader.h"

/// @todo remove
// Dangerous define
//#define size_t uint32_t

#ifndef __PBWT_EXP_HPP__
#define __PBWT_EXP_HPP__

template<typename T>
void print_vector(const std::vector<T>& vector) {
    std::cout << "[";
    for (size_t i = 0; i < vector.size()-1; ++i) {
        std::cout << vector[i] << ", ";
    }
    std::cout << vector.back() << "]" << std::endl;
}

template<typename T>
void print_differences(const std::vector<T>& v1, const std::vector<T>& v2) {
    for (size_t i = 0; i < v1.size(); ++i) {
        if (v1[i] != v2[i]) {
            std::cout << v1[i] << " != " << v2[i] << " at " << i << std::endl;
        }
    }
}

template<typename T = bool, const bool VERBOSE = false>
std::vector<std::vector<T> > read_from_macs_file(const std::string& filename) {
    T conv[256];
    /// @todo conversion for more than 2 alleles
    conv[(size_t)'0'] = 0;
    conv[(size_t)'1'] = 1;

    std::fstream s(filename, s.in);
    if (!s.is_open()) {
        std::cout << "Failed to open " << filename << std::endl;
        return {};
    } else {
        // Getting the header
        /////////////////////
        constexpr size_t line_size_c = 4096;
        char line[line_size_c];
        s.getline(line, line_size_c);
        std::string line_s(line);

        if (line_s.find("COMMAND:") == std::string::npos) {
            std::cout << "MaCS COMMAND line not found" << std::endl;
        }
        std::stringstream ss(line);
        std::vector<std::string> tokens;
        std::string word;
        while(getline(ss, word, ' ')) {
            tokens.push_back(word);
        }

        // Getting the number of haplotypes (samples)
        /////////////////////////////////////////////
        size_t n_samples = stoull(tokens[1]); // Number haplotypes
        if constexpr (VERBOSE) std::cout << "Number of samples is : " << n_samples << std::endl;

        // Getting the number of sites
        //////////////////////////////
        const auto position = s.tellg();
        // Ghetto way of getting number of sites
        /// @todo Make this more robust
        s.seekg(-25, std::ios_base::end);
        size_t last_site = 0;
        s >> last_site;
        size_t n_sites = last_site + 1;
        if constexpr (VERBOSE) std::cout << "sites : " << n_sites << std::endl;

        // Filling the hap map
        //////////////////////
        std::vector<std::vector<T> > hap_map(n_sites, std::vector<T>(n_samples));

        s.seekg(position); // Go back
        s.getline(line, line_size_c); // Get seed line

        const size_t line_length = n_samples + 100;
        char site_line[line_length];
        size_t current_site = 0;
        while(current_site < n_sites) {
            s.getline(site_line, line_length);
            /// @todo optimize this (too many new objects for little to nothing)
            std::stringstream ss((std::string(site_line)));
            /// @todo this is ghetto too
            ss.seekg(5); // Ignore "SITE:"
            {size_t _; ss >> _;} // Throw away int value
            {double _; ss >> _; ss >>_;} // Throw away two float values
            while(ss.peek() != '0' and ss.peek() != '1') {
                ss.get();
            }

            for (size_t i = 0; i < n_samples; ++i) {
                const uint8_t c = ss.get();
                hap_map[current_site][i] = conv[c];
                //std::cout << c;
            }
            //std::cout << std::endl;
            current_site++;
        }
        if constexpr (VERBOSE) std::cout << "Sites parsed : " << current_site << std::endl;
        return hap_map;
    }
}

template <typename T = bool>
std::vector<std::vector<T> > read_from_bcf_file(const std::string& filename, size_t max_m = 0, size_t max_n = 0) {
    std::vector<std::vector<T> > hap_map;

    // https://github.com/samtools/htslib/blob/develop/htslib/synced_bcf_reader.h
    bcf_srs_t * sr = bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;

    bcf_sr_add_reader(sr, filename.c_str());
    size_t n_samples = bcf_hdr_nsamples(sr->readers[0].header);
    size_t count = 0;

    unsigned int nset = 0;
    int *gt_arr = NULL;
    int ngt_arr = 0;
    bcf1_t* line; // https://github.com/samtools/htslib/blob/develop/htslib/vcf.h
    while ((nset = bcf_sr_next_line(sr))) {
        if (max_m and (count >= max_m)) {
            std::cout << "Stopped because of max m : " << max_m << std::endl;
            break;
        }
        line = bcf_sr_get_line(sr, 0 /* file */);
        if (line->n_allele != 2) {
            std::cerr << "Number of alleles is different than 2" << std::endl;
        } else {
            hap_map.push_back(std::vector<T>(0));
            bcf_unpack(line, BCF_UN_STR);
            //std::string chr = bcf_hdr_id2name(sr->readers[0].header, line->rid);
            //int pos = line->pos+1;
            //std::string id = std::string(line->d.id);
            //std::string ref = std::string(line->d.allele[0]);
            //std::string alt = std::string(line->d.allele[1]);
            //unsigned int cref = 0;
            //unsigned int calt = 0;
            int ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt_arr);

            int line_max_ploidy = ngt / n_samples;
            size_t sample_counter = 0;

            for (size_t i = 0; i < n_samples; ++i) {
                if (max_n and (sample_counter >= max_n)) {
                    break;
                }
                int32_t* ptr = gt_arr + i * line_max_ploidy;
                for (int j = 0; j < 2 /* ploidy, @todo check */; ++j) {
                    bool a = (bcf_gt_allele(ptr[j]) == 1);
                    hap_map.back().push_back(a);
                    sample_counter++;
                }
            }
        }
        count++;
    }
    free(gt_arr);
    bcf_sr_destroy(sr);

    return hap_map;
}

typedef std::vector<size_t> ppa_t;
typedef std::vector<size_t> d_t;

typedef struct a_d_arrays_at_pos {
    size_t pos;
    ppa_t a;
    d_t d;
} a_d_arrays_at_pos;

using hap_map_t = std::vector<std::vector<bool> >;
//using hap_map_t = HapMapStdVector<bool>; // Do not use ! Slow

template <const bool VERIFY = false>
inline void fill_rppa(ppa_t& rppa, const ppa_t& ppa) {
    if constexpr(VERIFY) {
        if (ppa.size() != rppa.size()) {
            std::cerr << "vector sizes don't match" << std::endl;
        }
    }

    const size_t N = ppa.size();
    for (size_t i = 0; i < N; ++i) {
        rppa[ppa[i]] = i;
    }
}

typedef struct match_t {
    size_t a;
    size_t b;
    size_t start;
    size_t end; // Note this is k and could be inferred from data structure position (e.g., if all the matches for a given k are in the same data structure)
} match_t;

using matches_t = std::vector<match_t>;

static matches_t __place_holder__; // Trick to still use reference even with no parameter

// Algorithm 2 in Durbin 2014
/// @TODO DOXY
inline void algorithm_2_BuildPrefixAndDivergenceArrays(const std::vector<bool>& x, const size_t& k, ppa_t& a, ppa_t& b, d_t& d, d_t& e) {
    const size_t M = x.size();
    size_t u = 0, v = 0, p = k+1, q = k+1;

    for (size_t i = 0; i < M; ++i) {
        if (d[i] > p) p = d[i];
        if (d[i] > q) q = d[i];
        if (x[a[i]] == 0) { // y[i] is x[a[i]]
            a[u] = a[i]; d[u] = p; u++; p = 0;
        } else {
            b[v] = a[i]; e[v] = q; v++; q = 0;
        }
    }

    // Concatenations
    std::copy(b.begin(), b.begin()+v, a.begin()+u);
    std::copy(e.begin(), e.begin()+v, d.begin()+u);
}

// Algorithm 3 in Durbin 2014
/// @TODO DOXY
inline
void algorithm_3_ReportLongMatches(const std::vector<bool>& x, const size_t N, const size_t& k, const size_t& L, const ppa_t& a, d_t& d, const std::function<void (size_t ai, size_t bi, size_t start, size_t end)> &report) {
    size_t u = 0;
    size_t v = 0;
    size_t ia = 0;
    size_t ib = 0;
    size_t i0 = 0;
    size_t dmin = 0;
    ppa_t a_(N);
    ppa_t b_(N);

    const size_t M = a.size();

    for (size_t i = 0; i < M; ++i) {
        if (d[i] > k-L) {
            if (u && v) { /* Then there is something to report */
                for (ia = i0; ia < i; ++ia) {
                    for (ib = ia+1, dmin = 0; ib < i; ++ib) {
                        if (d[ib] > dmin) dmin = d[ib];
                        // y[i] is x[a[i]]
                        if (x[a[ib]] != x[a[ia]]) {
                            report(a[ia], a[ib], dmin, k);
                        }
                    }
                } /* end for all */
            }
            u = 0;
            v = 0;
            i0 = i;
        }
        if (x[a[i]] == 0) {
            u++;
        } else {
            v++;
        }
    }
}


// Algorithm 4 in Durbin 2014
/// @TODO DOXY
inline
void algorithm_4_ReportSetMaximalMatches(const std::vector<bool>& x, const size_t N, const size_t& k, const ppa_t& a, d_t& d, const std::function<void (size_t ai, size_t bi, size_t start, size_t end)> &report) {
    // Sentinels
    d[0] = k+1;
    d.push_back(k+1);

    const size_t M = a.size();

    for (size_t i = 0; i < M; ++i) {
        size_t m = i-1;
        size_t n = i+1;

        // Scan down the array
        if (d[i] <= d[i+1]) {
            while (d[m+1] <= d[i]) {
                if (k < N and x[a[m]] == x[a[i]]) { // y[i] is x[a[i]]
                    goto next_i;
                }
                m--;
            }
        }

        // Scan up the array
        if (d[i] >= d[i+1]) {
            while (d[n] <= d[i+1]) {
                if (k < N and x[a[n]] == x[a[i]]) { // y[i] is x[a[i]]
                    goto next_i;
                }
                n++;
            }
        }

        // Reporting
        //mutex.lock();
        for (size_t j = m+1; j < i; ++j) {
            // Report (if length > 0)
            if (d[i] != k) {
                report(a[i], a[j], d[i], k);
                //printf("MATCH\t%d\t%d\t%d\t%d\t%d\n", a[i], a[j], d[i], k, k-d[i]);
                //fprintf(pFile, "MATCH\t%d\t%d\t%d\t%d\t%d\n", a[i], a[j], d[i], k, k-d[i]);
                //matches.push_back({
                //    .a = a[i],
                //    .b = a[j],
                //    .start = d[i],
                //    .end = k
                //});
            }
        }
        for (size_t j = i+1; j < n; ++j) {
            // Report (if length > 0)
            if (d[i+1] != k) {
                report(a[i], a[j], d[i+1], k);
                //printf("MATCH\t%d\t%d\t%d\t%d\t%d\n", a[i], a[j], d[i+1], k, k-d[i+1]);
                //fprintf(pFile, "MATCH\t%d\t%d\t%d\t%d\t%d\n", a[i], a[j], d[i+1], k, k-d[i+1]);
                //matches.push_back({
                //    .a = a[i],
                //    .b = a[j],
                //    .start = d[i+1],
                //    .end = k
                //});
            }
        }
        //mutex.unlock();

        next_i:
        ;
    }
    d.pop_back(); // Remove the extra sentinel
}

// This fixes a and d between two indexes [start, stop[ given the previous a and d's
// Algorithm 2 in paper
inline
void fix_a_d_range(const size_t& start, const size_t& stop,
                   const ppa_t& prev_a, const d_t& prev_d, const ppa_t& rppa,
                   ppa_t& a, d_t& d) {
    ppa_t prev_pos_of_a_s_to_fix; // "Arr" in paper

    // Here the a's from "start" to "stop-1" must be ordered given the previous permutations so we get their previous positions
    for (size_t j = start; j < stop; ++j) {
        //std::cout << "fixing id = " << a[j] << " at position " << j << std::endl;
        //std::cout << "position given previous sort is " << rppa[a[j]] << std::endl;
        prev_pos_of_a_s_to_fix.push_back(rppa[a[j]]);
    }

    // Sort the positions, this gives the actual ordering
    std::sort(prev_pos_of_a_s_to_fix.begin(), prev_pos_of_a_s_to_fix.end());

    // Actually fix the a's
    for (size_t j = start; j < stop; ++j) {
        a[j] = prev_a[prev_pos_of_a_s_to_fix[j-start]];
        //std::cout << "fixed id = " << a[j] << std::endl;
        //std::cout << "now at position " << j << std::endl;
    }

    /// @note both operations above could be done with std::sort and a lambda that sorts a given prev_a

    // Fix the d's
    // This requires "scanning" equivalent to the symbol counters of the multibit version
    // Because the fixing "coarse" d's is like having computed the d's on a multibit version
    // Scanning is not the most effective but is required to remove dependencies of previous stage
    // Scans should be short if data has underlying LD structure
    // The more "coarse it is" the smaller the groups that requires scans will be
    // i.e., more groups (more diversity because longer seqs) of smaller size => less scans, and because of LD, scans should not be too long (they may still be enourmous)
    /// @note The scans could maybe be optimized by boolean logic on the scan intervals
    // The first d is always correct because non 0, now fill the 0's with actual values
    for (size_t j = start+1; j < stop; ++j) {
        const size_t scan_start = prev_pos_of_a_s_to_fix[j-start-1] + 1;
        const size_t scan_stop = prev_pos_of_a_s_to_fix[j-start] + 1;

        //std::cout << "scan [start = " << scan_start << " scan stop = [" << scan_stop << std::endl;
        // Scan
        d[j] = *std::max_element(prev_d.begin() + scan_start, prev_d.begin() + scan_stop);
        //std::cout << "fixed d at " << j << " to : " << d[j] << std::endl;
    }
}

// Algorithm 1 in paper
template <const bool DEBUG=false>
void fix_a_d(std::vector<a_d_arrays_at_pos>& a_d_arrays) {
    const size_t M = a_d_arrays[0].a.size(); // Number of haplotypes
    ppa_t rppa(M); // Reverse look-up array
    ppa_t group;

    size_t debug_fix_counter = 0;

    // This first a and d arrays computed are always correct
    for (size_t _ = 1; _ < a_d_arrays.size(); ++_) {
        // To be fixed
        auto& a = a_d_arrays[_].a;
        auto& d = a_d_arrays[_].d;
        // Previous are fixed and can now be used
        const auto& prev_a = a_d_arrays[_-1].a;
        const auto& prev_d = a_d_arrays[_-1].d;

        // Create the reverse ppa (not needed if no fix => optimize this)
        fill_rppa(rppa, prev_a);

        const size_t start = a_d_arrays[_-1].pos; // Started from the previous position
        size_t first_same_seq_index = 0;

        for (size_t i = 0; i < M; ++i) {
            // Check if same sequence as previous
            if (d[i] == start) {
                // ... Nothing to do ...
            } else {
                // Here we have a new sequence, if the last sequence group was bigger than 1 then it may need fixing
                const size_t last_group_size = i - first_same_seq_index;
                if (last_group_size > 1) {
                    fix_a_d_range(first_same_seq_index, i, prev_a, prev_d, rppa, a, d);
                    if constexpr (DEBUG) std::cout << "Fixing between " << first_same_seq_index << " and " << i << std::endl;
                    if constexpr (DEBUG) debug_fix_counter++;
                }

                // Current sequence is different from previous, so this is a new group, update index
                first_same_seq_index = i;
            }
        }

        const size_t last_group_size = M - first_same_seq_index;
        if (last_group_size > 1) { // may need fixing
            fix_a_d_range(first_same_seq_index, M, prev_a, prev_d, rppa, a, d);
            if constexpr (DEBUG) debug_fix_counter++;
        }
    }
    if constexpr (DEBUG) std::cout << "Number of fixes : " << debug_fix_counter << std::endl;
}

// PBWT Process Matrix from position start to position stop (if stop is 0 process to the end)
// If report matches is active this algorithm performs ReportSetMaximalMatches, else it doesn't (resolved at compile time)
template<const bool REPORT_MATCHES = false>
a_d_arrays_at_pos process_matrix_sequentially(const hap_map_t& hap_map, const size_t start, const size_t stop = 0, const ppa_t& given_a = {}, const d_t& given_d = {}, const std::function<void (size_t ai, size_t bi, size_t start, size_t end)> &report = nullptr) {
    const size_t effective_stop = stop ? stop : hap_map.size();
    const size_t N = hap_map.size(); // Number of variant sites (total)
    const size_t M = hap_map[0].size(); // Number of haplotypes

    ppa_t a(M), b(M);
    if (given_a.size()) { // If a is given copy it
        std::copy(given_a.begin(), given_a.end(), a.begin());
    } else { // Else use natural order (0,1,2,...,M-1)
        std::iota(a.begin(), a.end(), 0);
    }

    d_t d(M), e(M);
    if (given_d.size()) { // If d is given copy it
        std::copy(given_d.begin(), given_d.end(), d.begin());
    } else { // Else use a d filled with "start"
        std::fill(d.begin(), d.end(), start);
    }

    for (size_t k = start; k < effective_stop; ++k) {
        if constexpr (REPORT_MATCHES) algorithm_4_ReportSetMaximalMatches(hap_map[k], N, k, a, d, report);
        algorithm_2_BuildPrefixAndDivergenceArrays(hap_map[k], k, a, b, d, e);
    }
    if constexpr (REPORT_MATCHES) {
        if (effective_stop == N) {
            algorithm_4_ReportSetMaximalMatches(hap_map[N-1], N, N, a, d, report); // Special case
        }
    }

    // Return a and d at stop position
    return {effective_stop, a, d};
}

//#define EXPERIMENTAL
#ifdef EXPERIMENTAL
////// EXPERIMENTAL
// Algorithm 2 in Durbin 2014
/// @TODO DOXY
inline void algorithm_2_BuildPrefixAndDivergenceArrays_linked_list(const std::vector<bool>& x, const size_t& k, size_t& start_a, ppa_t& a, ppa_t& b, d_t& d, d_t& e) {
    const size_t M = x.size();
    size_t u = 0, v = 0, p = k+1, q = k+1;

    size_t next = start_a;
    start_a = -1; // No start yet
    size_t last_a = -1;
    size_t start_b = -1;
    size_t last_b = -1;
    if (d[0] > p) p = d[0];
    if (d[0] > q) q = d[0];
    if (x[next] == 0) { // y[i] is x[a[i]]
        // Create list
        start_a = next; // Nop
        last_a = next;
        d[u] = p; u++; p = 0;
    } else {
        // Create list
        start_b = next;
        last_b = next;
        e[v] = q; v++; q = 0;
    }
    next = a[next];
    for (size_t i = 1; i < M; ++i) {
        if (d[i] > p) p = d[i];
        if (d[i] > q) q = d[i];
        if (x[next] == 0) { // y[i] is x[a[i]]
            //a[u] = a[i];
            if (start_a == -1) {
                // Create list
                start_a = next;
            } else {
                // Link list
                a[last_a] = next;
            }
            last_a = next;
            d[u] = p; u++; p = 0;
        } else {
            //b[v] = a[i];
            if (start_b == -1) {
                // Create list
                start_b = next;
            } else {
                // Link list
                a[last_b] = next;
            }
            last_b = next;
            e[v] = q; v++; q = 0;
        }
        next = a[next];
    }

    // Concatenations
    //std::copy(b.begin(), b.begin()+v, a.begin()+u);
    if (start_a == -1) {
        // Only 1's
        start_a = start_b;
    } else {
        if (start_b != -1) {
            // 0's and 1's, link the two lists
            a[last_a] = start_b;
        } else {
            // Only 0's
            // Nothing to do
        }
    }
    std::copy(e.begin(), e.begin()+v, d.begin()+u);
}

// PBWT Process Matrix from position start to position stop (if stop is 0 process to the end)
// If report matches is active this algorithm performs ReportSetMaximalMatches, else it doesn't (resolved at compile time)
template<const bool REPORT_MATCHES = false>
a_d_arrays_at_pos process_matrix_sequentially_linked_list(const hap_map_t& hap_map, const size_t start, const size_t stop = 0, const ppa_t& given_a = {}, const d_t& given_d = {}, const std::function<void (size_t ai, size_t bi, size_t start, size_t end)> &report = nullptr) {
    const size_t effective_stop = stop ? stop : hap_map.size();
    const size_t N = hap_map.size(); // Number of variant sites (total)
    const size_t M = hap_map[0].size(); // Number of haplotypes

    size_t start_a = 0;
    ppa_t a(M), b(M);
    if (given_a.size()) { // If a is given copy it
        start_a = given_a.at(0);
        size_t current_pos = start;
        for (size_t i = 1; i < M; ++i) {
            // Current := next (linked list conversion)
            a[current_pos] = given_a.at(i);
            current_pos = given_a.at(i);
        }
        a[current_pos] = -1; // End of linked list
    } else { // Else use natural order (0,1,2,...,M-1)
        start_a = 0;
        std::iota(a.begin(), a.end(), 1); // Next element 1,2,...,M-2
    }

    d_t d(M), e(M);
    if (given_d.size()) { // If d is given copy it
        std::copy(given_d.begin(), given_d.end(), d.begin());
    } else { // Else use a d filled with "start"
        std::fill(d.begin(), d.end(), start);
    }

    // Commented code below because a is now a linked list
    for (size_t k = start; k < effective_stop; ++k) {
        //if constexpr (REPORT_MATCHES) algorithm_4_ReportSetMaximalMatches(hap_map[k], N, k, a, d, report);
        algorithm_2_BuildPrefixAndDivergenceArrays_linked_list(hap_map[k], k, start_a, a, b, d, e);
    }
    //if constexpr (REPORT_MATCHES) {
    //    if (effective_stop == N) {
    //        algorithm_4_ReportSetMaximalMatches(hap_map[N-1], N, N, a, d, report); // Special case
    //    }
    //}

    // Return a and d at stop position
    ppa_t ppa(M);
    ppa[0] = start_a;
    size_t next = a[start_a];
    for (size_t i = 1; i < M; ++i) {
        ppa[i] = next;
        next = a[next];
    }
    return {effective_stop, ppa, d};
}
//////
#endif /* EXPERIMENTAL */

std::vector<size_t> generate_positions_to_collect(const size_t N, const size_t THREADS) {
    if (THREADS == 0) {std::cerr << "THREADS should be at least 1" << std::endl; return {};} // Edge case, should not happen

    const size_t step = N / THREADS;

    std::vector<size_t> positions_to_collect(THREADS-1, 0);
    for (size_t i = 1; i < THREADS; ++i) {
        positions_to_collect[i-1] = step * i;
    }

    return positions_to_collect;
}

std::vector<a_d_arrays_at_pos> generate_a_d_arrays_for_positions_sequentially(const hap_map_t& hap_map, const std::vector<size_t>& positions) {
    std::vector<a_d_arrays_at_pos> a_d_arrays;

    if (positions.empty()) return a_d_arrays;

    const size_t N = hap_map.size(); // Number of variant sites to process
    const size_t M = hap_map[0].size(); // Number of haplotypes
    ppa_t a(M), b(M);
    std::iota(a.begin(), a.end(), 0);
    d_t d(M, 0), e(M);

    auto pos_iterator = positions.begin();

    // Traverse variant sites of hap_map (returns when all requested positions have been collected)
    for (size_t k = 0; k < N; ++k) {
        if (pos_iterator == positions.end()) { // If all positions have been generated break
            return a_d_arrays;
        } else {
            if (*pos_iterator == k) { // If we want to collect the a and d arrays at position k
                a_d_arrays.push_back({k, a, d});
                pos_iterator++;
            }
        }

        // Apply algorithm 2 at position k to generate the next a and d arrays
        algorithm_2_BuildPrefixAndDivergenceArrays(hap_map[k], k, a, b, d, e);
    }

    // If the final is requested (used for examples and test, in other cases the function returns earlier)
    if (*pos_iterator == N) {
        a_d_arrays.push_back({N, a, d});
        pos_iterator++; // Not necessary
    }

    return a_d_arrays;
}

std::vector<a_d_arrays_at_pos> generate_a_d_arrays_for_positions_in_parallel(const hap_map_t& hap_map, const std::vector<size_t>& positions) {
    const size_t M = hap_map[0].size(); // Number of haplotypes
    ppa_t a(M), b(M);
    std::iota(a.begin(), a.end(), 0);
    d_t d(M, 0), e(M);

    const size_t THREADS = positions.size();
    std::vector<std::thread> workers(THREADS);
    std::vector<a_d_arrays_at_pos> a_d_arrays(THREADS);

    if (positions.empty()) return a_d_arrays;

    for (size_t i = 0; i < THREADS; ++i) {
        workers[i] = std::thread([=, &hap_map, &a_d_arrays, &positions]{
            // If first thread start at 0 else start at ending position of last thread
            const size_t start = (i == 0) ? 0 : positions[i-1];
            const size_t stop = positions[i];
            a_d_arrays[i] = process_matrix_sequentially(hap_map, start, stop);
        });
    }
    for (auto& t : workers) {
        t.join();
    }

    // Fix possible errors sequentially
    fix_a_d(a_d_arrays);

    return a_d_arrays;
}

matches_t report_matches_sequentially(const hap_map_t& hap_map) {
    matches_t matches;
    auto f = [&](size_t ai, size_t bi, size_t start, size_t end){
        matches.push_back({
            .a = ai,
            .b = bi,
            .start = start,
            .end = end
        });
    };

    process_matrix_sequentially<true /* Report Matches */>(hap_map, 0, 0, {}, {}, f);

    return matches;
}

void report_matches_sequentially_to_file(const hap_map_t& hap_map, const std::string& ofname) {
    if (ofname.compare("-") == 0) {
        auto report = [](size_t ai, size_t bi, size_t start, size_t end){
            printf("MATCH\t%zu\t%zu\t%zu\t%zu\t%zu\n", ai, bi, start, end, end-start);
        };

        process_matrix_sequentially<true /* Report Matches */>(hap_map, 0, 0, {}, {}, report);
    } else {
        FILE* pFile = fopen(ofname.c_str(), "w"); // Old school file pointer

        // Reporting lambda function (copy captures the file pointer)
        auto report_to_file = [=](size_t ai, size_t bi, size_t start, size_t end){
            // Same reporting function as in Durbin2014
            fprintf(pFile, "MATCH\t%zu\t%zu\t%zu\t%zu\t%zu\n", ai, bi, start, end, end-start);
        };

        process_matrix_sequentially<true /* Report Matches */>(hap_map, 0, 0, {}, {}, report_to_file);
    }
}

std::vector<matches_t> report_matches_in_parallel_a_d_sequential(const hap_map_t& hap_map, const size_t THREADS = 1) {
    auto positions_to_collect = generate_positions_to_collect(hap_map.size(), THREADS);
    auto a_d_arrays = generate_a_d_arrays_for_positions_sequentially(hap_map, positions_to_collect);

    std::vector<std::thread> workers(THREADS);
    std::vector<matches_t> matches(THREADS);

    for (size_t i = 0; i < THREADS; ++i) {
        workers[i] = std::thread([=, &hap_map, &a_d_arrays, &matches]{
            auto f = [&](size_t ai, size_t bi, size_t start, size_t end){
                matches[i].push_back({
                    .a = ai,
                    .b = bi,
                    .start = start,
                    .end = end
                });
            };

            // If first thread start at 0 else start at ending position of last thread
            const size_t start = (i == 0) ? 0 : a_d_arrays[i-1].pos;
            // If last thread set stop to special value of 0 else stop at position
            const size_t stop = (i == THREADS-1) ? 0 : a_d_arrays[i].pos;
            if (i == 0) {
                // The first thread processes from natural order (empty a,d arrays given)
                process_matrix_sequentially<true /* Report */>(hap_map, start, stop, {/*a*/}, {/*d*/}, f);
            } else {
                process_matrix_sequentially<true /* Report */>(hap_map, start, stop, a_d_arrays[i-1].a, a_d_arrays[i-1].d, f);
            }
        });
    }
    for (auto& t : workers) {
        t.join();
    }

    return matches;
}

template<const bool TO_FILES = false>
std::vector<matches_t> report_matches_in_parallel(const hap_map_t& hap_map, const size_t THREADS = 1, const std::string ofname = "") {
    auto positions_to_collect = generate_positions_to_collect(hap_map.size(), THREADS);
    auto a_d_arrays = generate_a_d_arrays_for_positions_in_parallel(hap_map, positions_to_collect);

    std::vector<std::thread> workers(THREADS);
    std::vector<matches_t> matches(THREADS);

    for (size_t i = 0; i < THREADS; ++i) {
        workers[i] = std::thread([=, &hap_map, &a_d_arrays, &matches]{
            auto f = [&](size_t ai, size_t bi, size_t start, size_t end){
                matches[i].push_back({
                    .a = ai,
                    .b = bi,
                    .start = start,
                    .end = end
                });
            };

            std::stringstream filename;
            filename << ofname << "_" << i;
            FILE* pFile = nullptr;
            if constexpr (TO_FILES) {
                pFile = fopen(filename.str().c_str(), "w");
            }
            auto report_to_file = [=](size_t ai, size_t bi, size_t start, size_t end){
                fprintf(pFile, "MATCH\t%zu\t%zu\t%zu\t%zu\t%zu\n", ai, bi, start, end, end-start);
            };

            // If first thread start at 0 else start at ending position of last thread
            const size_t start = (i == 0) ? 0 : a_d_arrays[i-1].pos;
            // If last thread set stop to special value of 0 else stop at position
            const size_t stop = (i == THREADS-1) ? 0 : a_d_arrays[i].pos;
            if constexpr (TO_FILES) {
                if (i == 0) {
                    // The first thread processes from natural order (empty a,d arrays given)
                    process_matrix_sequentially<true /* Report */>(hap_map, start, stop, {/*a*/}, {/*d*/}, report_to_file);
                } else {
                    process_matrix_sequentially<true /* Report */>(hap_map, start, stop, a_d_arrays[i-1].a, a_d_arrays[i-1].d, report_to_file);
                }
                fclose(pFile);
            } else {
                if (i == 0) {
                    // The first thread processes from natural order (empty a,d arrays given)
                    process_matrix_sequentially<true /* Report */>(hap_map, start, stop, {/*a*/}, {/*d*/}, f);
                } else {
                    process_matrix_sequentially<true /* Report */>(hap_map, start, stop, a_d_arrays[i-1].a, a_d_arrays[i-1].d, f);
                }
            }
        });
    }
    //std::cerr << "Running with " << THREADS << " threads " << std::endl;
    for (auto& t : workers) {
        t.join();
    }

    return matches;
}

#endif /* __PBWT_EXP_HPP__ */