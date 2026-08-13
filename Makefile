CXX      ?= clang++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -pthread
INCLUDES := -Isrc

BIN := validate gcsa bench_repetition simulate_repeats compare_algos ilp_baseline

all: $(BIN)

validate: src/validate.cpp src/shape.hpp src/gapped_sa.hpp src/sais.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) src/validate.cpp -o $@

gcsa: src/main.cpp src/shape.hpp src/gapped_sa.hpp src/sais.hpp src/compress.hpp src/serialize.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) src/main.cpp -o $@

bench_repetition: src/bench_repetition.cpp src/shape.hpp src/gapped_sa.hpp src/sais.hpp src/compress.hpp src/serialize.hpp src/random_dna.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) src/bench_repetition.cpp -o $@

simulate_repeats: src/simulate_repeats.cpp src/random_dna.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) src/simulate_repeats.cpp -o $@

compare_algos: src/compare_algos.cpp src/shape.hpp src/gapped_sa.hpp src/sais.hpp src/compress.hpp src/serialize.hpp src/random_dna.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) src/compare_algos.cpp -o $@

ilp_baseline: src/ilp_baseline.cpp src/shape.hpp src/gapped_sa.hpp src/sais.hpp src/compress.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) src/ilp_baseline.cpp -o $@

clean:
	rm -f $(BIN)

.PHONY: all clean
