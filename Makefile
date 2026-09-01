# Build the GAMB-backed FoXS command-line executable.
#
# Override dependency locations when needed, for example:
#   make GAMB_DIR=/path/to/gamb BOOST_PREFIX=/path/to/boost

CXX ?= c++
TARGET ?= foxs
BUILD_DIR ?= build

GAMB_DIR ?= ../gamb
GAMB_LIB ?= $(GAMB_DIR)/libgamb++.a
BOOST_PREFIX ?= /opt/homebrew/opt/boost

CPPFLAGS += -I. -I$(GAMB_DIR) -I$(BOOST_PREFIX)/include
CXXFLAGS ?= -O3
CXXFLAGS += -std=c++17 -Wall -Wextra -MMD -MP
LDFLAGS += -L$(BOOST_PREFIX)/lib
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
DEPS := $(OBJECTS:.o=.d)

.PHONY: all clean check-deps test

all: $(TARGET)

# Standalone adaptation of IMP's modules/foxs/test/test_saxs.py::test_simple.
test: $(TARGET)
	python3 tests/test_foxs.py ./$(TARGET)

$(TARGET): check-deps $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) $(LDLIBS) -o $@

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

check-deps:
	@test -f "$(GAMB_LIB)" || { echo "Missing GAMB library: $(GAMB_LIB)"; exit 1; }
	@test -d "$(BOOST_PREFIX)/include/boost" || { echo "Missing Boost headers under $(BOOST_PREFIX)/include"; exit 1; }

clean:
	rm -rf "$(BUILD_DIR)" "$(TARGET)"

-include $(DEPS)
