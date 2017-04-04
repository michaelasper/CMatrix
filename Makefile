#
# Makefile for 'matrices'.
#
# Type 'make' or 'make matrices' to create the binary.
# Type 'make clean' or 'make clear' to delete all temporaries.
# Type 'make run' to execute the binary.
# Type 'make debug' to debug the binary using gdb(1).
#

# build target specs
CC = gcc
CFLAGS = -g -Wall -Werror -DDEBUG -D_DEBUG 
OUT_DIR = .
LIBS = -lm

# first target entry is the target invoked when typing 'make'
default: matrices

matrices: $(OUT_DIR)/matrices.c.o
	@echo -n 'Linking matrices... '
	@$(CC) $(CFLAGS) -o matrices $(OUT_DIR)/matrices.c.o $(LIBS)
	@echo Done.

$(OUT_DIR)/matrices.c.o: matrices.c
	@echo -n 'Compiling matrices.c... '
	@$(CC) $(CFLAGS) -o $(OUT_DIR)/matrices.c.o -c matrices.c
	@echo Done.

run:
	./matrices 

debug:
	gdb ./matrices

clean:
	@echo -n 'Removing all temporary binaries... '
	@rm -f matrices $(OUT_DIR)/*.o
	@echo Done.

clear:
	@echo -n 'Removing all temporary binaries... '
	@rm -f matrices $(OUT_DIR)/*.o
	@echo Done.

