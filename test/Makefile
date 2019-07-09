INCLUDE := -I../include
CC      := cc
CFLAGS  := -Wall -Wextra -std=c99 -g -fno-omit-frame-pointer -fsanitize=address
LINK    := -lm -lcriterion -L./ -lsymspg

INCS    := $(wildcard ../include/*.h)
SOURCE  := $(wildcard ./*.c)
OBJS    := $(patsubst %.c, %.o, $(SOURCE))
EXE     := $(patsubst %.c, %.out, $(SOURCE))

LIBSRC  := $(wildcard ../src/*.c)
LIBOBJ	:= $(patsubst ../src/%.c, ../build/%.o, $(LIBSRC))

.SECONDARY:
	./*.o ../build/*

all: $(EXE)
	
%.out: %.o $(LIBOBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LINK)

%.o: %.c $(INCS)
	$(CC) -c -o $@ $< $(INCLUDE) $(CFLAGS)

../build/%.o: ../src/%.c $(INCS)
	$(CC) -c -o $@ $< $(INCLUDE) $(CFLAGS)

clean:
	@rm ./*.out ./*.o ../build/*
