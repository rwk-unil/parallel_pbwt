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
make app # Generates application
make benchmark # Runs benchmark
make test # Runs unit tests
```

```shell
./app -h
Report matches example application
Usage: ./app [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -f,--file TEXT              Input file name, default is stdio
  -o,--output TEXT            Output file name, default is stdio
  -t,--threads UINT           Number of threads
  -L,--length UINT            Length for long matches, 0 or undefined reports set maximal matches
```

### Comparison with PBWT (Durbin 2014)

https://github.com/richarddurbin/pbwt

Note: For demonstration purposes this software only reads bi-allelic sites from VCF/BCF.

#### Report long matches

```shell
./pbwt -readVcfGT <input_file.vcf/bcf> -longWithin <length> > results.txt
```

```shell
./app -f <input_file.vcf/bcf> -L <length> -o results.txt -t <number of threads>
```

#### Report set maximal matches

```shell
./pbwt -readVcfGT <input_file.vcf/bcf> -maxWithin > results.txt
```

```shell
./app -f <input_file.vcf/bcf> -o results.txt -t <number of threads>
```