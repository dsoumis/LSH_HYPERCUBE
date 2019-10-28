#include <iostream>
#include "InputImplementation.h"
#include "GridImplementation.h"
#include "utilitiesLSHCurves.h"
#include "HashFunctions.h"
#include "HashTables.h"
#include "BruteForceImplementation.h"
#define dimensions 2

int main(int argc, char **argv) {

    string inputFile,queryFile,outputFile;
    int k_vec=4,L_grid=4;
    assign_values_from_arguments(argc,argv,inputFile,queryFile,k_vec,L_grid,outputFile);
    cout<<inputFile<<" "<<queryFile<<" "<<k_vec<<" "<<L_grid<<" "<<outputFile<<endl;
    while (inputFile.empty()) {
        cout << "Please give input file." << endl;
        getline(cin, inputFile);
    }
    while (outputFile.empty()) {
        cout << "Please give output file." << endl;
        getline(cin, outputFile);
    }
    while(true) {
        while (queryFile.empty()) {
            cout << "Please give query file." << endl;
            getline(cin, queryFile);
        }
        unsigned int maxCurveSize = 0, minCurveSize = 4294967295;
        InputGenericVector<pair<double, double>> pointsVector(inputFile, maxCurveSize, minCurveSize, true);

        InputGenericVector<pair<double, double>> queriesVector(inputFile, maxCurveSize, minCurveSize, false);
        maxCurveSize *= 2; //Because its vector will consist of x1,y1,x2,y2,x3,y3 not in pairs



        double delta = 0.00184445;//deltaCalculator(pointsVector);
        double max_coord = 0; //Consider max coordinate to be zero
        pointsVector.maxCoordFinder(max_coord); //Find  the max coordinate to use it in vector padding
        queriesVector.maxCoordFinder(max_coord);
        max_coord *= 100;
        max_coord /= delta;


        ExactNeighboursVector<pair<double, double>> neighboursVector(pointsVector, queriesVector, false);

        vector<Grid<pair<double, double>>> grids;
        grids.reserve((unsigned long) L_grid);

        vector<HashFunctions<double>> hashFunctions;
        hashFunctions.reserve((unsigned long) L_grid);

        vector<HashTables<pair<double, double>>> hashTables;
        hashTables.reserve((unsigned long) L_grid);

        vector<double> w; //One w for each grid
        w.resize((unsigned long) L_grid);

        vector<vector<vector<double>>> concatenated_gridCurve_vectors;
        concatenated_gridCurve_vectors.resize((unsigned long) L_grid);

        for (unsigned int v = 0; v < L_grid; v++) {
            concatenated_gridCurve_vectors[v].resize(pointsVector.itemValues.size());
        }


        for (unsigned int i = 0; i < L_grid; i++) {
            grids.emplace_back(Grid<pair<double, double>>(minCurveSize, dimensions, delta));//Create L_grid grids

            //Create the concatenated grid curve vector for each curve of each grid
            for (unsigned int curve = 0; curve < pointsVector.itemValues.size(); curve++) {
                concatenated_gridCurve_vectors[i][curve] = grids[i].concatenated_vector_from_gridCurve(
                        pointsVector.itemValues[curve].second, maxCurveSize, max_coord);
            }
            //Exact Neighbours here calculates only 5% of neighbors for speed
            ExactNeighboursVector<double> temp_concatenatedNeighbors(concatenated_gridCurve_vectors[i],
                                                                     concatenated_gridCurve_vectors[i], true);
            w[i] = temp_concatenatedNeighbors.wCalculator();
            //cout<<w[i]<<endl;

            hashFunctions.emplace_back(
                    HashFunctions<double>(k_vec, maxCurveSize, w[i]));//Create L_grid g functions. One for each grid.
            hashTables.emplace_back(HashTables<pair<double, double>>(1)); //Create L_grid LSHs with one hash table each.
        }


        curves_points_LSH(L_grid, pointsVector, hashFunctions, hashTables, (unsigned int) k_vec, w,
                          concatenated_gridCurve_vectors);

        curves_queries_LSH(L_grid, queriesVector, grids, hashFunctions, maxCurveSize, max_coord, (unsigned int) k_vec,
                           w, hashTables, outputFile, neighboursVector);
        cout << "Do you want to repeat the process with another query file?(Y/N)" << endl;
        string answer;
        while (answer != "N" && answer != "Y")
            getline(cin, answer);
        if (answer == "N")
            break;
        queryFile.erase();
    }
    return 5;
}