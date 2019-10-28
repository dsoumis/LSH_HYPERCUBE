#ifndef LSH_UTILITIES_H
#define LSH_UTILITIES_H
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <string.h>
#include "InputImplementation.h"
#include "HashFunctions.h"
#include "HashTables.h"
#include "GridImplementation.h"
#include "BruteForceImplementation.h"
using namespace std;
void assign_values_from_arguments(int const& argc,char **argv,string& inputFile,string& queryFile,int& k,int& L,string& outputFile);
double deltaCalculator(InputGenericVector<pair<double,double>> const& point);
void curves_points_LSH(int const& L_grid,InputGenericVector<pair<double,double>>& pointsVector,vector<HashFunctions<double>>& hashFunctions,vector<HashTables<pair<double,double>>>& hashTables,
               unsigned int const& k_vec,vector<double> const& w,vector<vector<vector<double>>> const& concatenated_gridCurve_vectors);
void curves_queries_LSH(int const& L_grid,InputGenericVector<pair<double,double>>& queriesVector,vector<Grid<pair<double,double>>>& grids,vector<HashFunctions<double>>& hashFunctions,
                        unsigned int const& maxCurveSize,double const& max_coord,unsigned int const& k_vec,vector<double> const& w,vector<HashTables<pair<double,double>>>& hashTables,string const& outputFile,
                        ExactNeighboursVector<pair<double,double>>& neighboursVector);
#endif //LSH_UTILITIES_H
