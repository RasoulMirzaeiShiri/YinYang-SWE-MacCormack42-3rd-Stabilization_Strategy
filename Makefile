```makefile
# Compiler
FC = gfortran

# Compilation flags
FFLAGS = -O3

# Executable name
TARGET = strategy3

# Source files
SRC = $(wildcard *.f90)

# Object files
OBJ = $(SRC:.f90=.o)

# Build executable
$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJ)

# Compile source files
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Remove generated files
clean:
	rm -f *.o *.mod $(TARGET)

# Rebuild from scratch
rebuild: clean $(TARGET)

# Run the code
run: $(TARGET)
	./$(TARGET)

.PHONY: clean rebuild run
```
