OBJS_curve_grid_lsh    = BruteForceImplementation.o HashFunctions.o HashTables.o InputImplementation.o main.o utilitiesLSHCurves.o GridImplementation.o
SOURCE_curve_grid_lsh  = BruteForceImplementation.cpp HashFunctions.cpp HashTables.cpp InputImplementation.cpp main.cpp utilitiesLSHCurves.cpp GridImplementation.cpp
HEADER_curve_grid_lsh  = BruteForceImplementation.h HashFunctions.h HashTables.h InputImplementation.h utilitiesLSHCurves.h GridImplementation.h
OUT_curve_grid_lsh     = curve_grid_lsh

CC      = g++
FLAGS   = -g -c -O3 -std=c++11
# -g option enables debugging mode
# -c flag generates object code for separate files

all: 	$(OBJS_curve_grid_lsh)
		$(CC) -g $(OBJS_curve_grid_lsh) -o $(OUT_curve_grid_lsh)
		
BruteForceImplementation.o: BruteForceImplementation.cpp
	$(CC) $(FLAGS) BruteForceImplementation.cpp

HashFunctions.o: HashFunctions.cpp
	$(CC) $(FLAGS) HashFunctions.cpp

HashTables.o: HashTables.cpp
	$(CC) $(FLAGS) HashTables.cpp

InputImplementation.o: InputImplementation.cpp
	$(CC) $(FLAGS) InputImplementation.cpp

main.o: main.cpp
	$(CC) $(FLAGS) main.cpp

utilitiesLSHCurves.o: utilitiesLSHCurves.cpp
	$(CC) $(FLAGS) utilitiesLSHCurves.cpp	
	
GridImplementation.o: GridImplementation.cpp
	$(CC) $(FLAGS) GridImplementation.cpp
	
# clean house
clean:
	rm -f $(OBJS_curve_grid_lsh) $(OUT_curve_grid_lsh)