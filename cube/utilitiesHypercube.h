#ifndef LSH_UTILITIES_H
#define LSH_UTILITIES_H
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <string.h>
#include "InputImplementation.h"
#include "HashFunctions.h"
#include "CubeHashTable.h"
#include "fFunction.h"
#include "BruteForceImplementation.h"
using namespace std;
void assign_values_from_arguments(int const& argc,char **argv,string& inputFile,string& queryFile,int& k,int& M,int& probes,string& outputFile);
void pointsHypercube(InputGenericVector<int>& pointsVector,CubeHashTable<int>& cubeHashTable,vector<HashFunctions<int>>& hashFunctions,unsigned int const& k,
        double const& w,int const& d,fFunctions<int>& functions);
void queriesHypercube(InputGenericVector<int>& queriesVector,CubeHashTable<int>& cubeHashTable,vector<HashFunctions<int>>& hashFunctions,unsigned int const& k, double const& w,int const& d,
                      fFunctions<int>& functions,string const& outputFile,int const& M,int const& probes,ExactNeighboursVector<int>& neighboursVector,double const& radius);
#endif //LSH_UTILITIES_H
