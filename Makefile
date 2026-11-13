PROG = cellcoal
DASH = -
VERSION = 1.2.0
SRCDIR = src
OBJDIR = out
BINDIR = bin

CC ?= gcc

_OBJS = cellcoal.o eigen.o signatures.o
OBJS=$(patsubst %,$(OBJDIR)/%,$(_OBJS))

# -ffp-contract=off keeps floating-point results bit-identical to the original code
# (do NOT add -ffast-math or -march=native). OpenMP is optional: make OMPFLAG=
OMPFLAG ?= -fopenmp
CFLAGS ?= -Wall $(PROFILING) -O3 -ffp-contract=off $(OMPFLAG) $(WARN) -DUSE_COLORS

LIBS = -lm

CPPFLAGS += -D_XOPEN_SOURCE=500 -D_GNU_SOURCE

all: $(BINDIR)/$(PROG)-$(VERSION)

$(BINDIR)/$(PROG)-$(VERSION) : $(OBJS)
	mkdir -p $(BINDIR)
	$(CC) $(LDFLAGS) $(OMPFLAG) -o $@ $+ $(LIBS)

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