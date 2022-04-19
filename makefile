PY_DIST = $(shell which python3)
MAKEFLAGS += --no-print-directory
MPI_NPROC = 4
MPI_FLAGS = 


MPI_minP = 4
MPI_maxP = 8
MPI_stepP = 2

# The default target is 'all'. It will be executed when the makefile is
# called with make.

# 1. make all - to compile, run a simulation and plot the results.
all:
	@$(MAKE) -C src/
	@make run
	@$(MAKE) plot

## 2. make compile - to compile the program.
compile:
	@$(MAKE) -C src/

## 3. make run - to run.
run:
	@echo "Running the simulation with $(MPI_NPROC) processors."
	mpirun -np $(MPI_NPROC) $(MPI_FLAGS) ./src/program.exe

## 4. make plot - to plot the graphics
plot:
	@echo "Python found in: ${PY_DIST}. Plotting now."
	@${PY_DIST} ./src/modules/plot_statistics.py

## 5. make clean - to remove build and program.exe
clean:
	rm -r src/build
	rm -r src/program.exe

## 6. check paralel - perform parallel benchmarking
check_paralel:
	@echo "Checking parallel performance."
	./src/modules/check_paralel.sh -m $(MPI_maxP) -d $(MPI_minP) -s $(MPI_stepP)

.PHONY : help
