OBJS_lsh    = BruteForceImplementation.o HashFunctions.o HashTables.o InputImplementation.o main.o utilitiesLSH.o
SOURCE_lsh  = BruteForceImplementation.cpp HashFunctions.cpp HashTables.cpp InputImplementation.cpp main.cpp utilitiesLSH.cpp
HEADER_lsh  = BruteForceImplementation.h HashFunctions.h HashTables.h InputImplementation.h utilitiesLSH.h
OUT_lsh     = lsh

CC      = g++
FLAGS   = -g -c -O3 -std=c++11
# -g option enables debugging mode
# -c flag generates object code for separate files

all: 	$(OBJS_lsh)
		$(CC) -g $(OBJS_lsh) -o $(OUT_lsh)
		
# create/compile the individual files >>separately<<
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

utilitiesLSH.o: utilitiesLSH.cpp
	$(CC) $(FLAGS) utilitiesLSH.cpp
	
# clean house
clean:
	rm -f $(OBJS_lsh) $(OUT_lsh)