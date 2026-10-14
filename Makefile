# Build with OpenCoarrays:  brew install opencoarrays   (macOS)
ifeq ($(origin FC),default)
FC := caf
endif
FFLAGS  ?= -O2 -fcheck=bounds
NP      ?= 4
NPMAX   ?= 10       # for test-np
CONFIG  ?= config/default.toml

# BLAS (cgemm, cdotc, scnrm2): Accelerate on macOS, -lblas elsewhere
ifeq ($(shell uname),Darwin)
LIBS ?= -framework Accelerate
else
LIBS ?= -lblas
endif
BUILD   := build
BIN     := $(BUILD)/su3

# order matters (module dependencies)
LIB_SRC := mod_su3_matrix mod_utils mod_toml mod_config mod_logical_field_var \
           mod_real_field_var mod_su3_field_var mod_complex_field_var \
           mod_gauge_field mod_interfaces mod_metropolis mod_measure
LIB_OBJ := $(addprefix $(BUILD)/,$(addsuffix .o,$(LIB_SRC)))

# tests: one program per file in tests/ (the lattice is a global singleton,
# so each test needs its own process)
TESTS    := test_toml test_decomposition test_cold_start test_hot_start test_metropolis test_distribution
TEST_BIN := $(addprefix $(BUILD)/,$(TESTS))

all: $(BIN)

$(BIN): $(LIB_OBJ) $(BUILD)/main.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

$(BUILD)/%.o: src/%.f95 | $(BUILD)
	$(FC) $(FFLAGS) -J$(BUILD) -c $< -o $@

$(BUILD)/mod_testing.o: tests/mod_testing.f95 | $(BUILD)
	$(FC) $(FFLAGS) -J$(BUILD) -c $< -o $@

$(BUILD)/test_%: tests/test_%.f95 $(LIB_OBJ) $(BUILD)/mod_testing.o
	$(FC) $(FFLAGS) -J$(BUILD) -o $@ $< $(LIB_OBJ) $(BUILD)/mod_testing.o $(LIBS)

$(BUILD):
	mkdir -p $@

run: $(BIN)
	cafrun -np $(NP) ./$(BIN) $(CONFIG)

test: $(TEST_BIN)
	@for t in $(TESTS); do \
	    echo "== $$t"; \
	    cafrun -np $(NP) ./$(BUILD)/$$t || exit 1; \
	done
	@echo "All tests passed"

# run the distribution tests for every image count 1..NPMAX
test-np: $(BUILD)/test_distribution $(BUILD)/test_metropolis
	@for n in $$(seq 1 $(NPMAX)); do \
	    for t in test_distribution test_metropolis; do \
	        echo "== $$t, $$n image(s)"; \
	        cafrun -np $$n ./$(BUILD)/$$t > $(BUILD)/test-np.log 2>&1 || { cat $(BUILD)/test-np.log; exit 1; }; \
	        grep -E "FAIL|ok" $(BUILD)/test-np.log | grep -c ok | sed 's/^/   checks ok: /'; \
	    done; \
	done
	@echo "All image counts passed"

clean:
	rm -rf $(BUILD)

.PHONY: all run test test-np clean
