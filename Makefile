PROG = cellcoal
DASH = -
VERSION = 2.2.7
SRCDIR = src
OBJDIR = out
BINDIR = bin

CC ?= gcc

_OBJS = cellcoal.o globals.o tree.o mutation.o observation.o reads.o reads_output.o output.o rng.o options.o config.o sweep.o winexe.o bcf.o bgzf.o eigen.o signatures.o
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
# On Windows (MSYS2 / MinGW-w64) the test below cannot write to /dev/null, and libgomp always ships with the compiler
ifneq ($(findstring MINGW,$(shell uname -s))$(findstring MSYS,$(shell uname -s)),)
  ifeq ($(origin OMPFLAG),undefined)
    OMPFLAG = -fopenmp
  endif
endif
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

# zlib compresses the BGZF blocks of the BCF output (--bcf); without it the BCF is written uncompressed
ifneq ($(shell printf '\043include <zlib.h>\nint main(void){return 0;}\n' | $(CC) -x c - -lz -o /dev/null 2>/dev/null && echo yes),)
CPPFLAGS += -DHAVE_ZLIB
LIBS += -lz
endif

CPPFLAGS += -D_XOPEN_SOURCE=500 -D_GNU_SOURCE
ifeq ($(shell uname -s),Darwin)
CPPFLAGS += -D_DARWIN_C_SOURCE
endif

all: $(BINDIR)/$(PROG)-$(VERSION)

$(BINDIR)/$(PROG)-$(VERSION) : $(OBJS)
	mkdir -p $(BINDIR)
	$(CC) $(LDFLAGS) $(OMP_LDFLAGS) -o $@ $+ $(LIBS)

$(OBJS): $(OBJDIR)/%.o : $(SRCDIR)/%.c $(wildcard $(SRCDIR)/*.h)
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

.PHONY: clean test install

PREFIX ?= /usr/local

# make install [PREFIX=/usr/local] [DESTDIR=...]: installs the program as "cellcoal", the manual and the examples
install: all
	install -d $(DESTDIR)$(PREFIX)/bin $(DESTDIR)$(PREFIX)/share/doc/cellcoal $(DESTDIR)$(PREFIX)/share/cellcoal
	install -m 755 $(BINDIR)/$(PROG)-$(VERSION) $(DESTDIR)$(PREFIX)/bin/$(PROG)
	install -m 755 tools/bulk_vaf.py $(DESTDIR)$(PREFIX)/bin/cellcoal-bulk-vaf
	install -m 755 tools/cellcoal_plots.py $(DESTDIR)$(PREFIX)/bin/cellcoal-plots
	install -m 644 README.md $(wildcard docs/CHANGES-$(VERSION).md) $(DESTDIR)$(PREFIX)/share/doc/cellcoal
	install -m 644 docs/cellcoal.manual.v2.0.html $(DESTDIR)$(PREFIX)/share/doc/cellcoal
	cp -R examples parameters $(DESTDIR)$(PREFIX)/share/cellcoal/

test: all
	./tests/run_tests.sh $(BINDIR)/$(PROG)$(DASH)$(VERSION)

.PHONY: remove

clean:
	rm -f $(OBJS)
	rm -rf $(OBJDIR)

remove:
	rm -rf $(BINDIR)