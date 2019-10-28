OBJS_cube    = BruteForceImplementation.o HashFunctions.o CubeHashTable.o InputImplementation.o main.o utilitiesHypercube.o fFunction.o
SOURCE_cube  = BruteForceImplementation.cpp HashFunctions.cpp CubeHashTable.cpp InputImplementation.cpp main.cpp utilitiesHypercube.cpp fFunction.cpp 
HEADER_cube  = BruteForceImplementation.h HashFunctions.h CubeHashTable.h InputImplementation.h utilitiesHypercube.h fFunction.h
OUT_cube     = cube

CC      = g++
FLAGS   = -g -c -O3 -std=c++11
# -g option enables debugging mode
# -c flag generates object code for separate files

all: 	$(OBJS_cube)
		$(CC) -g $(OBJS_cube) -o $(OUT_cube)
		
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
	
utilitiesHypercube.o: utilitiesHypercube.cpp
	$(CC) $(FLAGS) utilitiesHypercube.cpp

fFunction.o: fFunction.cpp
	$(CC) $(FLAGS) fFunction.cpp
	
# clean house
clean:
	rm -f $(OBJS_cube) $(OUT_cube)