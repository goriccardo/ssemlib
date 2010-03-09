LIBNAME = ssemlib.so

LIBSRCS = sseasinh.c ssefrexp.c sselog.c sseatan.c
SUFFIXES = .c

OBJS = sseasinh.o ssefrexp.o sselog.o sseatan.o

LIBS = 

CC = gcc
CCFLAGS = -fPIC -O2 -msse2 -g -Wall
LDFLAGS = -shared

all: ssemlib

ssemlib: $(OBJS)
	$(CC) $(LDFLAGS) -o $(LIBNAME) $(OBJS)

clean:
	rm -f $(OBJS)

.SUFFIXES: $(SUFFIXES)

.c.o:
	$(CC) $(CCFLAGS) -c $<
