# Build with OpenCoarrays:  brew install opencoarrays   (macOS)
ifeq ($(origin FC),default)
FC := caf
endif
FFLAGS  ?= -O2 -fcheck=bounds
NP      ?= 4
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
TESTS    := test_toml test_cold_start test_hot_start test_metropolis
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

clean:
	rm -rf $(BUILD)

.PHONY: all run test clean
