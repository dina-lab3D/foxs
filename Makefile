# Build the GAMB- and DockingLib-backed FoXS command-line executable.
#
# Override dependency locations when needed, for example:
#   make GAMB_DIR=/path/to/gamb DOCKING_LIB_DIR=/path/to/DockingLib \
#        BOOST_PREFIX=/path/to/boost
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

GAMB_DIR ?= ../gamb
GAMB_LIB ?= $(GAMB_DIR)/libgamb++.a
DOCKING_LIB_DIR ?= ../DockingLib
DOCKING_LIB ?= $(DOCKING_LIB_DIR)/libdockingLib.a
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

BOOST_CPPFLAGS := $(if $(strip $(BOOST_PREFIX)),-I$(BOOST_PREFIX)/include)
BOOST_LDFLAGS := $(if $(strip $(BOOST_PREFIX)),-L$(BOOST_PREFIX)/lib)

CPPFLAGS += -I. -I$(DOCKING_LIB_DIR) -I$(GAMB_DIR) $(BOOST_CPPFLAGS)
CXXFLAGS ?= -O3
CXXFLAGS += -std=c++17 -Wall -Wextra -MMD -MP
NVCCFLAGS ?= -O3 -std=c++17
ifneq ($(strip $(CUDA_ARCH)),)
  NVCCFLAGS += -arch=$(CUDA_ARCH)
endif
LDFLAGS += $(BOOST_LDFLAGS)
BOOST_LIBS ?= -lboost_program_options
LDLIBS += $(DOCKING_LIB) $(GAMB_LIB) $(BOOST_LIBS)

# macOS universal/c architecture override. Set it only when GAMB,
# DockingLib, and Boost were built for the same architecture:
#   make ARCH=x86_64
ifeq ($(HOST_OS),Darwin)
  ifneq ($(strip $(ARCH)),)
    CXXFLAGS += -arch $(ARCH)
    LDFLAGS += -arch $(ARCH)
  endif
endif

SOURCES := \
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

OBJECTS := $(SOURCES:%.cpp=$(BUILD_DIR)/%.o)
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

$(TARGET): check-deps $(DOCKING_LIB) $(OBJECTS)
	$(LINKER) $(LDFLAGS) $(OBJECTS) $(LDLIBS) -o $@

$(DOCKING_LIB):
	$(MAKE) -C "$(DOCKING_LIB_DIR)" CXX="$(CXX)" \
		GAMB_DIR="$(abspath $(GAMB_DIR))" BOOST_PREFIX="$(BOOST_PREFIX)"

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/internal/cuda_helpers.o: internal/cuda_helpers.cu
	@mkdir -p $(dir $@)
	$(NVCC) $(CPPFLAGS) $(NVCCFLAGS) -c $< -o $@

check-deps:
	@test -f "$(GAMB_LIB)" || { echo "Missing GAMB library: $(GAMB_LIB)"; exit 1; }
	@test -f "$(DOCKING_LIB_DIR)/Makefile" || { echo "Missing DockingLib source directory: $(DOCKING_LIB_DIR)"; exit 1; }
	@if [ -n "$(BOOST_PREFIX)" ]; then \
		test -d "$(BOOST_PREFIX)/include/boost" || { echo "Missing Boost headers under $(BOOST_PREFIX)/include"; exit 1; }; \
	fi
	@if [ "$(GPU)" = "1" ]; then command -v "$(NVCC)" >/dev/null || { echo "Missing CUDA compiler: $(NVCC)"; exit 1; }; fi

clean:
	rm -rf "build" "build-cuda" "$(TARGET)"

-include $(DEPS)
