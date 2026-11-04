#Simple MakeFile
CC      := gcc
CFLAGS  := -O2 -Wall -Wextra -std=c99
LDFLAGS := -lm

#Builds PolyMorph Binary
BIN := PolyMorph

#All necessary C files
SRC := \
	main.c \
	pdb.c \
	curve.c \
	curve_adapt.c \
	selection.c \
	expr.c \
	scaling.c \
	cli.c \
	curvature.c \
	curve_scale.c \
	place.c \
	centDiff.c \
	clash_sortgrid.c \
	steric_scale.c

OBJ := $(SRC:.c=.o)

.PHONY: all clean

all: $(BIN)

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(BIN) $(OBJ)
