#ifndef LSH_UTILITIES_H
#define LSH_UTILITIES_H
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <string.h>
#include "InputImplementation.h"
#include "HashTables.h"
#include "BruteForceImplementation.h"
#include <chrono>
#include <tuple>
using namespace std;
void assign_values_from_arguments(int const& argc,char **argv,string& inputFile,string& queryFile,int& k,int& L,string& outputFile);
void queriesLSH(InputGenericVector<int>& queriesVector,HashTables<int>& hashTables,vector<HashFunctions<int>> const& hashFunctions,unsigned int const& k, double const& w,
        ExactNeighboursVector<int>& neighboursVector,double const& radius,string const& outputFile);
void pointsLSH(InputGenericVector<int>& pointsVector,HashTables<int>& hashTables,int const& L,vector<HashFunctions<int>>& hashFunctions,unsigned int const& k, double const& w);
#endif //LSH_UTILITIES_H
