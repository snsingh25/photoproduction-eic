# Single top-level Makefile. No per-directory Makefiles.
#
#   make <name>           build src/**/<name>.cc or plots/**/<name>.cc -> bin/<name>
#   make all              build every .cc under src/ and plots/
#   make list             list all buildable targets
#   make clean            remove bin/
#
# Example:  make jetreco   -> bin/jetreco

include build.mk

BIN_DIR := bin

# Discover every .cc under src/ and plots/, keyed by basename.
SRC_CCS   := $(shell find src plots -name '*.cc' 2>/dev/null)
ALL_NAMES := $(sort $(basename $(notdir $(SRC_CCS))))
ALL_BINS  := $(addprefix $(BIN_DIR)/,$(ALL_NAMES))

# evtgen is the only script that actually links PYTHIA8 (it generates events).
# Everything else reads ROOT files and clusters with FastJet.
PYTHIA_NAMES := evtgen

.PHONY: all list clean $(ALL_NAMES)

all: $(ALL_BINS)

list:
	@echo "Buildable targets (make <name>):"
	@for n in $(ALL_NAMES); do \
	  src=$$(find src plots -name "$$n.cc" | head -1); \
	  printf "  %-35s  %s\n" "$$n" "$$src"; \
	done

clean:
	rm -rf $(BIN_DIR)

# Phony shortcut: `make jetreco` -> build bin/jetreco.
$(ALL_NAMES): %: $(BIN_DIR)/%

# Find the .cc for a given basename and compile it.
# Uses PYTHIA link flags only for names in PYTHIA_NAMES.
$(BIN_DIR)/%: $(shell :)
	@mkdir -p $(BIN_DIR)
	@src=$$(find src plots -name "$*.cc" | head -1); \
	if [ -z "$$src" ]; then echo "No .cc found for target '$*'" >&2; exit 1; fi; \
	if echo " $(PYTHIA_NAMES) " | grep -q " $* "; then \
	  link="$(LINK_PYTHIA_ROOT_FASTJET)"; \
	else \
	  link="$(LINK_ROOT_FASTJET)"; \
	fi; \
	echo "[CXX] $$src -> $@"; \
	$(CXX) $(CXX_COMMON) "$$src" -o "$@" $$link
