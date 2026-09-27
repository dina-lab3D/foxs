# Build the standalone FoXS command-line executable.
#
# The minimal GAMB and DockingLib sources needed by FoXS are vendored under
# lib/. Boost remains an external dependency; override its location with:
#   make BOOST_PREFIX=/path/to/boost
#
# On macOS, Homebrew Boost is detected automatically in its standard Apple
# Silicon and Intel locations. On Linux, system compiler/library paths are
# used unless BOOST_PREFIX is provided.

CXX ?= c++
NVCC ?= nvcc
TARGET ?= foxs
GPU ?= 0
BUILD_DIR ?= $(if $(filter 1,$(GPU)),build-cuda,build)
HOST_OS := $(shell uname -s)

LIB_DIR := lib
BOOST_PREFIX ?=

# Find standard Homebrew installations without requiring brew to be on PATH.
ifeq ($(HOST_OS),Darwin)
  ifeq ($(strip $(BOOST_PREFIX)),)
    ifneq ($(wildcard /opt/homebrew/opt/boost/include/boost),)
      BOOST_PREFIX := /opt/homebrew/opt/boost
    else
      ifneq ($(wildcard /usr/local/opt/boost/include/boost),)
        BOOST_PREFIX := /usr/local/opt/boost
      endif
    endif
  endif
endif

BOOST_CPPFLAGS := $(if $(strip $(BOOST_PREFIX)),-isystem $(BOOST_PREFIX)/include)
BOOST_LDFLAGS := $(if $(strip $(BOOST_PREFIX)),-L$(BOOST_PREFIX)/lib)

CPPFLAGS += -I. -I$(LIB_DIR) $(BOOST_CPPFLAGS)
CXXFLAGS ?= -O3
CXXFLAGS += -std=c++17 -Wall -Wextra -MMD -MP
NVCCFLAGS ?= -O3 -std=c++17
ifneq ($(strip $(CUDA_ARCH)),)
  NVCCFLAGS += -arch=$(CUDA_ARCH)
endif
LDFLAGS += $(BOOST_LDFLAGS)
BOOST_LIBS ?= -lboost_program_options
LDLIBS += $(BOOST_LIBS)

# macOS architecture override. Set it only when Boost was built for the same
# architecture:
#   make ARCH=x86_64
ifeq ($(HOST_OS),Darwin)
  ifneq ($(strip $(ARCH)),)
    CXXFLAGS += -arch $(ARCH)
    LDFLAGS += -arch $(ARCH)
  endif
endif

FOXS_SOURCES := \
	foxs.cpp \
	Profile.cpp \
	Distribution.cpp \
	FormFactorTable.cpp \
	utility.cpp \
	ChiScore.cpp \
	ChiScoreLog.cpp \
	ChiFreeScore.cpp \
	RatioVolatilityScore.cpp \
	ProfileClustering.cpp \
	Gnuplot.cpp \
	JmolWriter.cpp \
	ColorCoder.cpp

DOCKING_LIB_SOURCES := \
	$(LIB_DIR)/ChemAtom.cc \
	$(LIB_DIR)/DotSphere.cc \
	$(LIB_DIR)/SASurface.cc

GAMB_SOURCES := \
	$(LIB_DIR)/Atom.cc \
	$(LIB_DIR)/CIF.cc \
	$(LIB_DIR)/Interface.cc \
	$(LIB_DIR)/PDB.cc \
	$(LIB_DIR)/SurfacePoint.cc

DEPENDENCY_SOURCES := $(DOCKING_LIB_SOURCES) $(GAMB_SOURCES)
OBJECTS := \
	$(FOXS_SOURCES:%.cpp=$(BUILD_DIR)/%.o) \
	$(DEPENDENCY_SOURCES:%.cc=$(BUILD_DIR)/%.o)
CUDA_OBJECTS :=

ifeq ($(GPU),1)
  CPPFLAGS += -DFOXS_SAXS_CUDA_LIB
  CUDA_OBJECTS += $(BUILD_DIR)/internal/cuda_helpers.o
  LINKER := $(NVCC)
else
  LINKER := $(CXX)
endif

OBJECTS += $(CUDA_OBJECTS)
DEPS := $(OBJECTS:.o=.d)

.PHONY: all clean check-deps test

all: $(TARGET)

# Standalone adaptation of IMP's modules/foxs/test/test_saxs.py::test_simple.
test: $(TARGET)
	python3 tests/test_foxs.py ./$(TARGET)

$(TARGET): $(OBJECTS) | check-deps
	$(LINKER) $(LDFLAGS) $(OBJECTS) $(LDLIBS) -o $@

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: %.cc
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/internal/cuda_helpers.o: internal/cuda_helpers.cu
	@mkdir -p $(dir $@)
	$(NVCC) $(CPPFLAGS) $(NVCCFLAGS) -c $< -o $@

check-deps:
	@test -f "$(LIB_DIR)/Atom.cc" || { echo "Missing vendored GAMB sources under $(LIB_DIR)"; exit 1; }
	@test -f "$(LIB_DIR)/SASurface.cc" || { echo "Missing vendored DockingLib sources under $(LIB_DIR)"; exit 1; }
	@if [ -n "$(BOOST_PREFIX)" ]; then \
		test -d "$(BOOST_PREFIX)/include/boost" || { echo "Missing Boost headers under $(BOOST_PREFIX)/include"; exit 1; }; \
	fi
	@if [ "$(GPU)" = "1" ]; then command -v "$(NVCC)" >/dev/null || { echo "Missing CUDA compiler: $(NVCC)"; exit 1; }; fi

clean:
	rm -rf "build" "build-cuda" "$(TARGET)"

-include $(DEPS)
