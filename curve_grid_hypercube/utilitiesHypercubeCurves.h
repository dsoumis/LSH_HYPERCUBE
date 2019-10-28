#ifndef LSH_UTILITIES_H
#define LSH_UTILITIES_H
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <string.h>
#include "InputImplementation.h"
#include "HashFunctions.h"
#include "fFunction.h"
#include "CubeHashTable.h"
#include "GridImplementation.h"
#include "BruteForceImplementation.h"
using namespace std;
void assign_values_from_arguments(int const& argc,char **argv,string& inputFile,string& queryFile,int& k,int& L,string& outputFile,int& M,int& probes);
double deltaCalculator(InputGenericVector<pair<double,double>> const& point);
void curves_points_HyperCube(int const& L_grid,InputGenericVector<pair<double,double>>& pointsVector,vector<vector<vector<double>>> const& concatenated_gridCurve_vectors,int const& d,
                             int const& k_hypercube,vector<double> const& w,vector<vector<HashFunctions<double>>>& hashFunctions,fFunctions<double>& f,vector<CubeHashTable<pair<double,double>>>& cubeHashTables);
void curves_queries_HyperCube(int const& L_grid,InputGenericVector<pair<double,double>>& queriesVector,vector<vector<vector<double>>> const& concatenated_gridCurve_vectors,int const& d,
                              int const& k_hypercube,vector<double> const& w,vector<vector<HashFunctions<double>>>& hashFunctions,fFunctions<double>& f,
                              vector<CubeHashTable<pair<double,double>>>& cubeHashTables,vector<Grid<pair<double,double>>>& grids,int const& M,int const& probes,
                              unsigned int const& maxCurveSize,double const& max_coord,string const& outputFile,ExactNeighboursVector<pair<double,double>>& neighboursVector);
#endif //LSH_UTILITIES_H
