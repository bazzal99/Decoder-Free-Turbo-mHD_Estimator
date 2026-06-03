# Makefile for turbo-mhd-estimator
#
# Targets:
#   make          — build the estimator (default)
#   make clean    — remove compiled objects and binary
#   make run      — build and run immediately

CC      = gcc
CFLAGS  = -O2 -std=c99 -Wall -Wextra
LDFLAGS = -lm

TARGET  = turbo_mhd
SRCS    = turbo_mhd.c trellis.c distance.c rtz_search.c utils.c
OBJS    = $(SRCS:.c=.o)

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

run: all
	./$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)

# Header dependencies
turbo_mhd.o:  turbo_mhd.c  trellis.h distance.h rtz_search.h utils.h
trellis.o:    trellis.c    trellis.h
distance.o:   distance.c   distance.h trellis.h utils.h
rtz_search.o: rtz_search.c rtz_search.h trellis.h distance.h utils.h
utils.o:      utils.c      utils.h
