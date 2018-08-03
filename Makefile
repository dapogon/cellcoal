PROG = cellcoal
SRCDIR = src
OBJDIR= out
BINDIR= bin

CC ?= gcc

_OBJS = cellcoal.o eigen.o signatures.o
OBJS=$(patsubst %,$(OBJDIR)/%,$(_OBJS))

CFLAGS ?= -Wall $(PROFILING) -O3 $(WARN) -DUSE_COLORS

LIBS = -lm

CPPFLAGS += -D_XOPEN_SOURCE=500

all: $(BINDIR)/$(PROG)

$(BINDIR)/$(PROG) : $(OBJS)
	mkdir -p $(BINDIR)
	$(CC) $(LDFLAGS) -o $@ $+ $(LIBS)

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