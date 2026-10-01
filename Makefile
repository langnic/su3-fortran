# Build with OpenCoarrays:  brew install opencoarrays   (macOS)
ifeq ($(origin FC),default)
FC := caf
endif
FFLAGS  ?= -O2 -fcheck=bounds
NP      ?= 4

# BLAS (cgemm, cdotc, scnrm2): Accelerate on macOS, -lblas elsewhere
ifeq ($(shell uname),Darwin)
LIBS ?= -framework Accelerate
else
LIBS ?= -lblas
endif
BUILD   := build
BIN     := $(BUILD)/su3

# order matters (module dependencies)
SRC := mod_su3_matrix mod_utils mod_logical_field_var mod_real_field_var \
       mod_su3_field_var mod_complex_field_var mod_gauge_field \
       mod_interfaces mod_metropolis mod_measure main
OBJ := $(addprefix $(BUILD)/,$(addsuffix .o,$(SRC)))

all: $(BIN)

$(BIN): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

$(BUILD)/%.o: src/%.f95 | $(BUILD)
	$(FC) $(FFLAGS) -J$(BUILD) -c $< -o $@

$(BUILD):
	mkdir -p $@

run: $(BIN)
	cafrun -np $(NP) ./$(BIN)

clean:
	rm -rf $(BUILD)

.PHONY: all run clean
