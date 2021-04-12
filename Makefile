BENCHMARK_PATH := ./benchmark
HTSLIB_PATH := ./htslib

# C++ Compiler
CXX=g++
INCLUDE_DIRS=-I . -I $(HTSLIB_PATH)/htslib -I $(BENCHMARK_PATH)/include -I $(BENCHMARK_PATH)/googletest/googletest/include
CXXFLAGS=-O3 -g -Wall -std=c++17 $(INCLUDE_DIRS) $(CXXEXTRAFLAGS)
# Linker
LD=g++
LIBS=-lpthread -lhts -lbenchmark
#LDFLAGS=-L benchmark/build/src/
LDFLAGS=-O3 -L $(HTSLIB_PATH) -L $(BENCHMARK_PATH)/build/src/

BENCHMARK_EXECUTABLE := pbwt_bench
BENCHMARK_SOURCE := main_bench.cpp
BENCHMARK_OBJ := $(BENCHMARK_SOURCE:.cpp=.o)

TEST_EXECUTABLE := pbwt_test
TEST_SOURCES := main_test.cpp $(BENCHMARK_PATH)/googletest/googletest/src/gtest_main.cc
TEST_OBJS := $(TEST_SOURCES:.cpp=.o)
TEST_OBJS := $(TEST_OBJS:.cc=.o)

CPP_SOURCES := $(wildcard *.cpp)
DEPENDENCIES := $(CPP_SOURCES:.cpp=.d)

# Rules
all : $(BENCHMARK_EXECUTABLE) $(TEST_EXECUTABLE)

benchmark : $(BENCHMARK_EXECUTABLE)
	./$(BENCHMARK_EXECUTABLE)

test : $(TEST_EXECUTABLE)
	./$(TEST_EXECUTABLE)

$(BENCHMARK_EXECUTABLE) : $(BENCHMARK_OBJ)
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

$(TEST_EXECUTABLE) : $(TEST_OBJS)
	$(LD) $(LDFLAGS) -L $(BENCHMARK_PATH)/build/lib -I $(BENCHMARK_PATH)/googletest/googletest/include $^ $(LIBS) -lgtest -o $@

# Do not include the depency rules for "clean"
ifneq ($(MAKECMDGOALS),clean)
-include $(DEPENDENCIES)
endif

# Compile
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to generate the dependency files
%.d : %.cpp
	$(CXX) $(INCLUDE_DIRS) -MG -MP -MM -MT '$(@:.d=.o)' $< -MF $@

# Remove artifacts
clean :
	rm -f *.o $(BENCHMARK_EXECUTABLE) $(BENCHMARK_OBJS)

# Rules that don't generate artifacts
.PHONY :
	all clean benchmark