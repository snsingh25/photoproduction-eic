# Shared build configuration for ROOT + FastJet + PYTHIA8 programs.
# Adjust these paths to match your installation, then all scripts build against them.
#
# ROOT is sourced from the pyRoot-env conda environment, which also has a
# working Python (numpy, uproot, matplotlib) — one environment covers both
# the C++ build/run side and the Python analysis side.
# Activate with:  conda activate pyRoot-env

PYTHIA_PREFIX := /Users/siddharthsingh/HEPTOOLS/PYTHIA/pythia8312/examples/pythia8311
PYTHIA_INCLUDE := -I$(PYTHIA_PREFIX)/include
PYTHIA_LIB     := -L$(PYTHIA_PREFIX)/lib -Wl,-rpath,$(PYTHIA_PREFIX)/lib -lpythia8 -ldl

FASTJET_INCLUDE := -I/usr/local/include
FASTJET_LIB     := -L/usr/local/lib -Wl,-rpath,/usr/local/lib -lfastjet

# Pin ROOT to the pyRoot-env conda env so the build is independent of whether
# a user activated it — the binary will still run fine either way because the
# rpath is baked in by `root-config --glibs`.
ROOT_CONFIG := /Users/siddharthsingh/miniforge3/envs/pyRoot-env/bin/root-config
ROOT_FLAGS  := $(shell $(ROOT_CONFIG) --cflags --glibs)

CXX        := clang++
CXX_COMMON := -O2 -std=c++17 -fPIC -pthread -w

# Preset link combos used by the top-level Makefile.
LINK_ROOT_FASTJET := $(FASTJET_INCLUDE) $(FASTJET_LIB) $(ROOT_FLAGS)
LINK_ROOT_ONLY    := $(ROOT_FLAGS)
LINK_PYTHIA_ROOT_FASTJET := $(PYTHIA_INCLUDE) $(PYTHIA_LIB) $(FASTJET_INCLUDE) $(FASTJET_LIB) $(ROOT_FLAGS)
