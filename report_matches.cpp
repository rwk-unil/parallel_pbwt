/*******************************************************************************
 * Copyright (C) 2021 Rick Wertenbroek, University of Lausanne
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include <iostream>
#include <filesystem>
#include <thread>
#include <chrono>
namespace fs = std::filesystem;

#include "CLI11.hpp"

#include "pbwt_exp.hpp"

void printElapsedTime(const std::chrono::steady_clock::time_point& begin,
                      const std::chrono::steady_clock::time_point& end);
void report_long_matches(const std::string& filename, const std::string& ofname, const size_t THREADS, const size_t L);
void report_set_maximal_matches(const std::string& filename, const std::string& ofname, const size_t THREADS);

int main(int argc, const char *argv[]) {

    CLI::App app{"Report matches example application"};

    std::string filename = "-";
    app.add_option("-f,--file", filename, "Input file name, default is stdio");
    std::string ofname = "-";
    app.add_option("-o,--output", ofname, "Output file name, default is stdio");
    //char O = 'u';
    //app.add_option("-O, --output-type", O, "output type b|u|z|v");
    size_t threads = 2;
    app.add_option("-t,--threads", threads, "Number of threads");
    size_t L = 0;
    app.add_option("-L,--length", L, "Length for long matches, 0 or undefined reports set maximal matches");

    CLI11_PARSE(app, argc, argv);
    if (threads == 0) {return -1;}

    std::cerr << "Running with " << threads << " threads" << std::endl;

    auto begin = std::chrono::steady_clock::now();
    if (L) {
        std::cerr << "Reporting long matches L=" << L << std::endl;
        report_long_matches(filename, ofname, threads, L);
    } else {
        std::cerr << "Reporting set maximal matches" << std::endl;
        report_set_maximal_matches(filename, ofname, threads);
    }
    auto end = std::chrono::steady_clock::now();

    std::cerr << "Total time : ";
    printElapsedTime(begin, end);

    return 0;
}

void printElapsedTime(const std::chrono::steady_clock::time_point& begin,
                      const std::chrono::steady_clock::time_point& end) {
    std::cerr << "Time elapsed = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s] "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms] "
              << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[us] "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "[ns]" << std::endl;
}

void report_long_matches(const std::string& filename, const std::string& ofname, const size_t THREADS, const size_t L) {
    if(filename.compare("-") and !fs::exists(filename)) {
        std::cerr << "File " << filename << " does not exist" << std::endl;
        exit(-1);
    }
    
    auto begin = std::chrono::steady_clock::now();
    auto hap_map = read_from_bcf_file(filename);
    auto end = std::chrono::steady_clock::now();
    std::cerr << "Reading BCF file : ";
    printElapsedTime(begin, end);

    if (filename.compare("-") == 0) {
        std::cerr << "Cannot output to stdout with multi-threaded version" << std::endl;
        exit(-1);
    }
    report_matches_in_parallel<true /* TO FILE */, 3>(hap_map, THREADS, ofname, L);
}

void report_set_maximal_matches(const std::string& filename, const std::string& ofname, const size_t THREADS) {
    if(filename.compare("-") and !fs::exists(filename)) {
        std::cerr << "File " << filename << " does not exist" << std::endl;
        exit(-1);
    }
    
    auto begin = std::chrono::steady_clock::now();
    auto hap_map = read_from_bcf_file(filename);
    auto end = std::chrono::steady_clock::now();
    std::cerr << "Reading BCF file : ";
    printElapsedTime(begin, end);

    if (THREADS == 1) {
        report_matches_sequentially_to_file(hap_map, ofname);
    } else {
        if (filename.compare("-") == 0) {
            std::cerr << "Cannot output to stdout with multi-threaded version" << std::endl;
            exit(-1);
        }
        report_matches_in_parallel<true /* TO FILE */>(hap_map, THREADS, ofname);
    }
}