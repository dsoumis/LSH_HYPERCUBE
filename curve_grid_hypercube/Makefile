OBJS_curve_grid_hypercube    = BruteForceImplementation.o HashFunctions.o CubeHashTable.o InputImplementation.o main.o utilitiesHypercubeCurves.o fFunction.o GridImplementation.o
SOURCE_curve_grid_hypercube  = BruteForceImplementation.cpp HashFunctions.cpp CubeHashTable.cpp InputImplementation.cpp main.cpp utilitiesHypercubeCurves.cpp fFunction.cpp GridImplementation.cpp 
HEADER_curve_grid_hypercube  = BruteForceImplementation.h HashFunctions.h CubeHashTable.h InputImplementation.h utilitiesHypercubeCurves.h fFunction.h GridImplementation.h
OUT_curve_grid_hypercube     = curve_grid_hypercube

CC      = g++
FLAGS   = -g -c -O3 -std=c++11
# -g option enables debugging mode
# -c flag generates object code for separate files

all: 	$(OBJS_curve_grid_hypercube)
		$(CC) -g $(OBJS_curve_grid_hypercube) -o $(OUT_curve_grid_hypercube)
		
BruteForceImplementation.o: BruteForceImplementation.cpp
	$(CC) $(FLAGS) BruteForceImplementation.cpp	
	
HashFunctions.o: HashFunctions.cpp
	$(CC) $(FLAGS) HashFunctions.cpp	
	
CubeHashTable.o: CubeHashTable.cpp
	$(CC) $(FLAGS) CubeHashTable.cpp	
	
InputImplementation.o: InputImplementation.cpp
	$(CC) $(FLAGS) InputImplementation.cpp	
	
main.o: main.cpp
	$(CC) $(FLAGS) main.cpp	
	
utilitiesHypercubeCurves.o: utilitiesHypercubeCurves.cpp
	$(CC) $(FLAGS) utilitiesHypercubeCurves.cpp	
	
fFunction.o: fFunction.cpp
	$(CC) $(FLAGS) fFunction.cpp	
	
GridImplementation.o: GridImplementation.cpp
	$(CC) $(FLAGS) GridImplementation.cpp
	
# clean house
clean:
	rm -f $(OBJS_curve_grid_hypercube) $(OUT_curve_grid_hypercube)