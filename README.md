# Parallel Positional Burrows-Wheeler Transform Algorithms

## Dependencies

- This software depends on [htslib](https://github.com/samtools/htslib) for reading VCF/BCF files (version 1.12).
- This software depends on [Google Benchmark](https://github.com/google/benchmark) for the benchmarks (version 1.5.2).

If the dependencies are already installed on the system the paths can be set accordingly in the Makefile. Else build the dependencies.

### Build the dependencies

```shell
# htslib
git submodule update --init --recursive htslib
cd htslib
autoreconf -i
./configure
make
cd ..

# benchmark
git submodule update --init benchmark
cd benchmark
git clone https://github.com/google/googletest.git # benchmark depends on Google Test
cmake -E make_directory "build"
cmake -E chdir "build" cmake -DCMAKE_BUILD_TYPE=Release ../
cmake --build "build" --config Release
cd ..
```

## Build and Run

```shell
make
make benchmark # Runs benchmark
make test # Runs unit tests
```