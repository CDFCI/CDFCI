# A simple Makefile for CDFCI

# compiler path (g++ or icpx)
CC = g++
# compiler flag
CCFLAG = -O3 -std=c++17
# compiler debug flag
DBGFLAG = -g
# Optional target-specific optimization flags, for example ARCH_FLAGS=-mavx.
# No x86 instruction set is enabled by default so ARM builds remain portable.
ARCH_FLAGS ?=
CCFLAG += $(ARCH_FLAGS)

# include path for header files
IFLAG = -I./src -I./include -I./src/lib/eigen -I./src/lib/pybind11/include
CCFLAG += $(IFLAG)

# compiler OpenMP flag
OMP_FLAG = -fopenmp
OMP_LFLAG =
# Apple Clang needs the separately installed libomp runtime.
ifneq ($(OS), Windows_NT)
ifeq ($(shell uname -s), Darwin)
        LIBOMP_PREFIX ?= $(shell brew --prefix libomp 2>/dev/null)
ifneq ($(LIBOMP_PREFIX),)
        OMP_FLAG = -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
        OMP_LFLAG = -L$(LIBOMP_PREFIX)/lib -Wl,-rpath,$(LIBOMP_PREFIX)/lib -lomp
else
        OMP_FLAG = -Xpreprocessor -fopenmp
        OMP_LFLAG = -lomp
endif
endif
endif
# Intel compiler uses -qopenmp
ifeq (ic, $(findstring ic, $(CC)))
        OMP_FLAG = -qopenmp
endif

# Ask config.h whether this compiler actually selected __float128/quadmath.
CDFCI_USE_FLOAT128 := $(shell $(CC) $(IFLAG) -dM -E -x c++ include/config.h 2>/dev/null | grep -q '^\#define CDFCI_USE_FLOAT128' && echo 1)
ifeq ($(CDFCI_USE_FLOAT128),1)
	LFLAG += -lquadmath
endif

# disable some warnings
ifeq (ic, $(findstring ic, $(CC)))
	CCFLAG += -Wno-deprecated-builtins -Wno-tautological-constant-compare
endif

DIR_BIN = bin
DIR_DEBUG = debug
PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin

# executable name
# OpenMP disabled, single thread
TARGET_CDFCI = cdfci
TARGET_TOOLS = cdfci_tools
TARGET_OPTORBFCI = optorbfci
TARGET_XCDFCI = xcdfci

# OpenMP enabled
TARGET_CDFCI_OMP = $(TARGET_CDFCI)_omp
# TARGET_OPTORBFCI_OMP = $(TARGET_OPTORBFCI)_omp
TARGET_XCDFCI_OMP = $(TARGET_XCDFCI)_omp

# source file containing the main() function
SOURCE_CDFCI = src/main.cpp
SOURCE_TOOLS = src/tools.cpp
SOURCE_OPTORBFCI = src/optorbfci.cpp
SOURCE_XCDFCI = src/xcdfci.cpp

# Python module
PYTHON        := python3
PY_INC        := $(shell $(PYTHON)-config --includes)
PY_LDFLAGS    := $(shell $(PYTHON)-config --ldflags)
PY_EXT_SUFFIX := $(shell $(PYTHON)-config --extension-suffix)

PYMOD_SRC     := src/pycdfci.cpp
PYMOD_OUT     := python/_cdfci$(PY_EXT_SUFFIX)

# build test
# compiler
# Use g++ for CI integration
CC_TEST = $(CC)
OMP_TEST_FLAG = $(OMP_FLAG)
# test directory
DIR_TEST = test
SOURCE_TEST = test_system.cpp
TARGET_TEST = cdfci_test
TARGET_OMP_TEST = $(TARGET_TEST)_omp

# examples directory
DIR_EXAMPLES = examples

# default rule
default: all

.PHONY: cdfci
cdfci:
	mkdir -p $(DIR_BIN)
	$(CC) $(CCFLAG) $(SOURCE_CDFCI) -o $(DIR_BIN)/$(TARGET_CDFCI) $(LFLAG)
	$(CC) $(CCFLAG) $(OMP_FLAG) $(SOURCE_CDFCI) -o $(DIR_BIN)/$(TARGET_CDFCI_OMP) $(LFLAG) $(OMP_LFLAG)

.PHONY: tools
tools:
	mkdir -p $(DIR_BIN)
	$(CC) $(CCFLAG) $(SOURCE_TOOLS) -o $(DIR_BIN)/$(TARGET_TOOLS) $(LFLAG)

.PHONY: optorbfci
optorbfci:
	mkdir -p $(DIR_BIN)
	$(CC) $(CCFLAG) $(SOURCE_OPTORBFCI) -o $(DIR_BIN)/$(TARGET_OPTORBFCI) $(LFLAG)
# 	$(CC) $(CCFLAG) $(OMP_FLAG) $(SOURCE_OPTORBFCI) -o $(DIR_BIN)/$(TARGET_OPTORBFCI_OMP) $(LFLAG)

.PHONY: xcdfci
xcdfci:
	mkdir -p $(DIR_BIN)
	$(CC) $(CCFLAG) $(SOURCE_XCDFCI) -o $(DIR_BIN)/$(TARGET_XCDFCI) $(LFLAG)
	$(CC) $(CCFLAG) $(OMP_FLAG) $(SOURCE_XCDFCI) -o $(DIR_BIN)/$(TARGET_XCDFCI_OMP) $(LFLAG) $(OMP_LFLAG)

.PHONY: all
all: cdfci tools xcdfci optorbfci

.PHONY: check
check: test

$(PYMOD_OUT): $(PYMOD_SRC)
	$(CC) $(CCFLAG) -shared -fPIC $(OMP_FLAG) $(PY_INC) $(PYMOD_SRC) -o $@ $(PY_LDFLAGS) $(LFLAG) $(OMP_LFLAG)

python-module: $(PYMOD_OUT)

clean-python:
	rm -f $(PYMOD_OUT)

.PHONY: debug
debug:
	mkdir -p $(DIR_DEBUG)
	$(CC) $(CCFLAG) $(DBGFLAG) $(SOURCE_CDFCI) -o $(DIR_DEBUG)/$(TARGET_CDFCI) $(LFLAG)
# 	$(CC) $(CCFLAG) $(DBGFLAG) $(OMP_FLAG) $(SOURCE_CDFCI) -o $(DIR_DEBUG)/$(TARGET_CDFCI_OMP) $(LFLAG)

.PHONY: build_test
build_test:
	$(CC_TEST) $(CCFLAG) $(DIR_TEST)/$(SOURCE_TEST) -o $(DIR_TEST)/$(TARGET_TEST) $(LFLAG)
	$(CC_TEST) $(CCFLAG) $(OMP_TEST_FLAG) $(DIR_TEST)/$(SOURCE_TEST) -o $(DIR_TEST)/$(TARGET_OMP_TEST) $(LFLAG) $(OMP_LFLAG)

.PHONY: test
test: build_test
	@echo "Running regression tests..."
	cd $(DIR_TEST) && ./$(TARGET_TEST) -D
ifeq ($(OS), Windows_NT)
	cd $(DIR_TEST) && set OMP_NUM_THREADS=2 && ./$(TARGET_OMP_TEST) -D
else
	cd $(DIR_TEST) && OMP_NUM_THREADS=2 ./$(TARGET_OMP_TEST) -D
endif

.PHONY: install
install: all
	install -d $(DESTDIR)$(BINDIR)
	install -m 0755 $(DIR_BIN)/$(TARGET_CDFCI) $(DESTDIR)$(BINDIR)/$(TARGET_CDFCI)
	install -m 0755 $(DIR_BIN)/$(TARGET_CDFCI_OMP) $(DESTDIR)$(BINDIR)/$(TARGET_CDFCI_OMP)
	install -m 0755 $(DIR_BIN)/$(TARGET_TOOLS) $(DESTDIR)$(BINDIR)/$(TARGET_TOOLS)
	install -m 0755 $(DIR_BIN)/$(TARGET_XCDFCI) $(DESTDIR)$(BINDIR)/$(TARGET_XCDFCI)
	install -m 0755 $(DIR_BIN)/$(TARGET_XCDFCI_OMP) $(DESTDIR)$(BINDIR)/$(TARGET_XCDFCI_OMP)
	install -m 0755 $(DIR_BIN)/$(TARGET_OPTORBFCI) $(DESTDIR)$(BINDIR)/$(TARGET_OPTORBFCI)

.PHONY: uninstall
uninstall:
	rm -f $(DESTDIR)$(BINDIR)/$(TARGET_CDFCI)
	rm -f $(DESTDIR)$(BINDIR)/$(TARGET_CDFCI_OMP)
	rm -f $(DESTDIR)$(BINDIR)/$(TARGET_TOOLS)
	rm -f $(DESTDIR)$(BINDIR)/$(TARGET_XCDFCI)
	rm -f $(DESTDIR)$(BINDIR)/$(TARGET_XCDFCI_OMP)
	rm -f $(DESTDIR)$(BINDIR)/$(TARGET_OPTORBFCI)

.PHONY: release-check
release-check: clean all test examples

examples: all
	@echo "Running all demo examples..."
	cd $(DIR_EXAMPLES) && chmod +x ./run_all.sh && ./run_all.sh

.PHONY: format

format:
	clang-format -i src/*.h
	clang-format -i src/*.cpp

.PHONY: clean
clean:
	rm -f $(DIR_BIN)/$(TARGET_CDFCI)
	rm -f $(DIR_BIN)/$(TARGET_CDFCI_OMP)
	rm -f $(DIR_BIN)/$(TARGET_TOOLS)
	rm -f $(DIR_BIN)/$(TARGET_XCDFCI)
	rm -f $(DIR_BIN)/$(TARGET_XCDFCI_OMP)
	rm -f $(DIR_BIN)/$(TARGET_OPTORBFCI)
	rm -f $(DIR_DEBUG)/$(TARGET_CDFCI)
	rm -f $(DIR_TEST)/$(TARGET_TEST)
	rm -f $(DIR_TEST)/$(TARGET_OMP_TEST)
