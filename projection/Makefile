OBJS_curve_projection    = main.o HashFunctions.o HashTables.o InputImplementation.o ProjectionArray.o
SOURCE_curve_projection   = main.cpp HashFunctions.cpp HashTables.cpp InputImplementation.cpp ProjectionArray.cpp
HEADER_curve_projection = HashFunctions.h HashTables.h InputImplementation.h ProjectionArray.h
OUT_curve_projection    = curve_projection

CC      = g++
FLAGS   = -g -c -O3 -std=c++11
# -g option enables debugging mode
# -c flag generates object code for separate files

all: 	$(OBJS_curve_projection)
		$(CC) -g $(OBJS_curve_projection) -o $(OUT_curve_projection)
		
	
HashTables.o: HashTables.cpp
	$(CC) $(FLAGS) HashTables.cpp	
	
HashFunctions.o: HashFunctions.cpp
	$(CC) $(FLAGS) HashFunctions.cpp	
	
InputImplementation.o: InputImplementation.cpp
	$(CC) $(FLAGS) InputImplementation.cpp	
	
main.o: main.cpp
	$(CC) $(FLAGS) main.cpp	
	
ProjectionArray.o: ProjectionArray.cpp
	$(CC) $(FLAGS) ProjectionArray.cpp	
	
# clean house
clean:
	rm -f $(OBJS_curve_projection) $(OUT_curve_projection)