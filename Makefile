PROG = cellcoal
DASH = -
VERSION = 2.0.0
SRCDIR = src
OBJDIR = out
BINDIR = bin

CC ?= gcc

_OBJS = cellcoal.o eigen.o signatures.o
OBJS=$(patsubst %,$(OBJDIR)/%,$(_OBJS))

# -ffp-contract=off keeps floating-point results bit-identical to the original code
# (do NOT add -ffast-math or -march=native). -fno-math-errno only avoids setting errno in libm calls (faster, same results).
#
# OpenMP (parallel read counts and VCF) is optional and detected automatically:
#   1. the compiler supports -fopenmp (gcc, clang on Linux);
#   2. otherwise, on macOS, Homebrew's libomp is used if installed (brew install libomp);
#   3. otherwise CellCoal is compiled without OpenMP (single thread, same results).
# Override with: make OMPFLAG=            (no OpenMP)
#                make OMPFLAG=-fopenmp    (force the flag for both compiling and linking)
ifeq ($(origin OMPFLAG),undefined)
  OMP_TEST_SRC = '\#include <omp.h>\nint main(void){return omp_get_max_threads()>0?0:1;}\n'
  ifeq ($(shell printf $(OMP_TEST_SRC) | $(CC) -fopenmp -x c - -o /dev/null >/dev/null 2>&1 && echo yes),yes)
    OMP_CFLAGS  = -fopenmp
    OMP_LDFLAGS = -fopenmp
  else
    LIBOMP_PREFIX ?= $(shell brew --prefix libomp 2>/dev/null)
    ifneq ($(LIBOMP_PREFIX),)
      ifeq ($(shell printf $(OMP_TEST_SRC) | $(CC) -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include -x c - -L$(LIBOMP_PREFIX)/lib -lomp -o /dev/null >/dev/null 2>&1 && echo yes),yes)
        OMP_CFLAGS  = -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
        OMP_LDFLAGS = -L$(LIBOMP_PREFIX)/lib -lomp
      endif
    endif
  endif
else
  OMP_CFLAGS  = $(OMPFLAG)
  OMP_LDFLAGS = $(OMPFLAG)
endif
ifeq ($(strip $(OMP_CFLAGS)),)
  OMP_CFLAGS = -Wno-unknown-pragmas
  $(info CellCoal: OpenMP not available, building single-threaded)
endif

CFLAGS ?= -Wall $(PROFILING) -O3 -ffp-contract=off -fno-math-errno $(OMP_CFLAGS) $(WARN) -DUSE_COLORS

LIBS = -lm

CPPFLAGS += -D_XOPEN_SOURCE=500 -D_GNU_SOURCE
ifeq ($(shell uname -s),Darwin)
CPPFLAGS += -D_DARWIN_C_SOURCE
endif

all: $(BINDIR)/$(PROG)-$(VERSION)

$(BINDIR)/$(PROG)-$(VERSION) : $(OBJS)
	mkdir -p $(BINDIR)
	$(CC) $(LDFLAGS) $(OMP_LDFLAGS) -o $@ $+ $(LIBS)

$(OBJS): $(OBJDIR)/%.o : $(SRCDIR)/%.c
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

.PHONY: clean 
.PHONY: remove

clean:
	rm -f $(OBJS)
	rm -rf $(OBJDIR)

remove:
	rm -rf $(BINDIR)