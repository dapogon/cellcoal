PROG = coaltum
SRCDIR = CoalTumor

CC ?= gcc

OBJS = coaltumor.o eigen.o
CFLAGS ?= -Wall $(PROFILING) -O3 $(WARN) -DUSE_COLORS

LIBS = -lm

CPPFLAGS += -D_XOPEN_SOURCE=500

all: $(PROG)

$(PROG) : $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $+ $(LIBS)

%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(PROG)
