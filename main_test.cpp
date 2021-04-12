#include <gtest/gtest.h>
#include "pbwt_exp.hpp"

#include <algorithm>

namespace {

    // Example from paper
    const std::vector<std::vector<bool> > TEST_MATRIX = {
        {0,1,0,0,1,1,0,0,0,1},
        {1,1,1,0,0,0,1,1,1,0},
        {1,0,0,1,1,1,1,0,0,1},
        {1,1,1,1,1,1,1,0,0,1},
        {1,0,0,0,0,0,1,1,0,0},
        {0,0,0,1,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,1,0,0},
        {0,1,1,1,0,0,0,0,0,0},
        {0,0,0,0,1,1,0,0,1,1},
        {0,1,1,0,0,1,1,1,0,0}
    };

    template <typename T>
    static void check_vectors(const std::vector<T>& v1, const std::vector<T>& v2) {
        ASSERT_EQ(v1.size(), v2.size());
        for (size_t i = 0; i < v1.size(); ++i) {
            EXPECT_EQ(v1[i], v2[i]);
        }
    }

    TEST(PBWT, EmptyTest) {
    }

    // Example from paper
    TEST(PBWT, CheckIfSequentialAndParallelADAreSameSmall) {
        std::vector<size_t> positions_to_collect = {5, 10};

        auto result_seq = generate_a_d_arrays_for_positions_sequentially(TEST_MATRIX, positions_to_collect);
        auto result_par = generate_a_d_arrays_for_positions_in_parallel(TEST_MATRIX, positions_to_collect);

        ASSERT_NE(result_seq.size(), 0);

        // Check if correct nubmer of results
        ASSERT_EQ(result_seq.size(), result_par.size());

        // Check if positions are the same
        EXPECT_TRUE(std::equal(result_seq.begin(), result_seq.end(),
                               result_par.begin(),
                               [](const a_d_arrays_at_pos& a, const a_d_arrays_at_pos& b) -> bool {
                                   return a.pos == b.pos;
                               }));

        // Check if a arrays are same
        for(size_t i = 0; i < result_seq.size(); ++i) {
            //print_vector(result_seq[i].a);
            //print_vector(result_seq[i].d);
            //print_vector(result_par[i].a);
            //print_vector(result_par[i].d);
            check_vectors(result_seq[i].a, result_par[i].a);
        }

        // Check if d arrays are same
        for(size_t i = 0; i < result_seq.size(); ++i) {
            check_vectors(result_seq[i].d, result_par[i].d);
        }
    }

    // Bigger example
    TEST(PBWT, CheckIfFirstPositionIsOK) {
        const size_t THREADS = 4;

        auto hap_map = read_from_macs_file<bool>("11k.macs"); // TODO CHANGE FILE
        const size_t N = hap_map.size(); // Number of variant sites
        auto positions_to_collect = generate_positions_to_collect(N, THREADS);

        positions_to_collect.resize(1);

        std::vector<a_d_arrays_at_pos> result_seq;
        std::vector<a_d_arrays_at_pos> result_par;

        result_seq = generate_a_d_arrays_for_positions_sequentially(hap_map, positions_to_collect);
        result_par = generate_a_d_arrays_for_positions_in_parallel(hap_map, positions_to_collect);

        auto first_res_seq = result_seq.front();
        auto first_res_par = result_par.front();
        EXPECT_EQ(first_res_seq.pos, first_res_par.pos);
        EXPECT_EQ(first_res_seq.a, first_res_par.a);
        print_differences(first_res_seq.a, first_res_par.a);
        EXPECT_EQ(first_res_seq.d, first_res_par.d);
        print_differences(first_res_seq.d, first_res_par.d);
    }

    // Bigger example
    TEST(PBWT, CheckIfSequentialAndParallelADAreSame) {
        const size_t THREADS = 4; /// @todo range

        auto hap_map = read_from_macs_file<bool>("11k.macs"); // TODO CHANGE FILE
        auto positions_to_collect = generate_positions_to_collect(hap_map.size(), THREADS);

        std::vector<a_d_arrays_at_pos> result_seq;
        std::vector<a_d_arrays_at_pos> result_par;

        result_seq = generate_a_d_arrays_for_positions_sequentially(hap_map, positions_to_collect);
        result_par = generate_a_d_arrays_for_positions_in_parallel(hap_map, positions_to_collect);

        ASSERT_NE(result_seq.size(), 0);

        // Check if positions are the same
        EXPECT_TRUE(std::equal(result_seq.begin(), result_seq.end(),
                               result_par.begin(),
                               [](const a_d_arrays_at_pos& a, const a_d_arrays_at_pos& b) -> bool {
                                   return a.pos == b.pos;
                               }));

        // Check if a arrays are same
        EXPECT_TRUE(std::equal(result_seq.begin(), result_seq.end(),
                               result_par.begin(),
                               [](const a_d_arrays_at_pos& a, const a_d_arrays_at_pos& b) -> bool {
                                   //print_differences(a.a, b.a);
                                   return std::equal(a.a.begin(), a.a.end(), b.a.begin());
                               }));

        // Check if d arrays are same
        EXPECT_TRUE(std::equal(result_seq.begin(), result_seq.end(),
                               result_par.begin(),
                               [](const a_d_arrays_at_pos& a, const a_d_arrays_at_pos& b) -> bool {
                                   return std::equal(a.d.begin(), a.d.end(), b.d.begin());
                               }));
    }
}