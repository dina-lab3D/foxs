# Build the GAMB-backed FoXS command-line executable.
#
# Override dependency locations when needed, for example:
#   make GAMB_DIR=/path/to/gamb BOOST_PREFIX=/path/to/boost

CXX ?= c++
NVCC ?= nvcc
TARGET ?= foxs
GPU ?= 0
BUILD_DIR ?= $(if $(filter 1,$(GPU)),build-cuda,build)

GAMB_DIR ?= ../gamb
GAMB_LIB ?= $(GAMB_DIR)/libgamb++.a
BOOST_PREFIX ?=

# Leave BOOST_PREFIX empty to use Boost from the system include and library
# paths. Set it only for a non-standard installation, e.g.
# BOOST_PREFIX=/opt/homebrew/opt/boost.
BOOST_CPPFLAGS := $(if $(strip $(BOOST_PREFIX)),-I$(BOOST_PREFIX)/include)
BOOST_LDFLAGS := $(if $(strip $(BOOST_PREFIX)),-L$(BOOST_PREFIX)/lib)

CPPFLAGS += -I. -I$(GAMB_DIR) $(BOOST_CPPFLAGS)
CXXFLAGS ?= -O3
CXXFLAGS += -std=c++17 -Wall -Wextra -MMD -MP
NVCCFLAGS ?= -O3 -std=c++17
ifneq ($(strip $(CUDA_ARCH)),)
  NVCCFLAGS += -arch=$(CUDA_ARCH)
endif
LDFLAGS += $(BOOST_LDFLAGS)
BOOST_LIBS ?= -lboost_program_options
LDLIBS += $(GAMB_LIB) $(BOOST_LIBS)

# Set ARCH explicitly only when all three dependencies were built for it:
#   make ARCH=x86_64
ifneq ($(strip $(ARCH)),)
  CXXFLAGS += -arch $(ARCH)
  LDFLAGS += -arch $(ARCH)
endif

SOURCES := \
	foxs.cpp \
	Profile.cpp \
	Distribution.cpp \
	FormFactorTable.cpp \
	SolventAccessibleSurface.cpp \
	utility.cpp \
	ChiScore.cpp \
	ChiScoreLog.cpp \
	ChiFreeScore.cpp \
	RatioVolatilityScore.cpp \
	ProfileClustering.cpp \
	RadiusOfGyrationRestraint.cpp \
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

$(TARGET): check-deps $(OBJECTS)
	$(LINKER) $(LDFLAGS) $(OBJECTS) $(LDLIBS) -o $@

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/internal/cuda_helpers.o: internal/cuda_helpers.cu
	@mkdir -p $(dir $@)
	$(NVCC) $(CPPFLAGS) $(NVCCFLAGS) -c $< -o $@

check-deps:
	@test -f "$(GAMB_LIB)" || { echo "Missing GAMB library: $(GAMB_LIB)"; exit 1; }
	@if [ -n "$(BOOST_PREFIX)" ]; then \
		test -d "$(BOOST_PREFIX)/include/boost" || { echo "Missing Boost headers under $(BOOST_PREFIX)/include"; exit 1; }; \
	else \
		test -d /usr/include/boost || { echo "Missing system Boost headers; install Boost or set BOOST_PREFIX=/path/to/boost"; exit 1; }; \
	fi
	@if [ "$(GPU)" = "1" ]; then command -v "$(NVCC)" >/dev/null || { echo "Missing CUDA compiler: $(NVCC)"; exit 1; }; fi

clean:
	rm -rf "build" "build-cuda" "$(TARGET)"

-include $(DEPS)
