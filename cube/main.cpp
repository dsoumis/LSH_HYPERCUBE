#include <iostream>
#include "InputImplementation.h"
#include "CubeHashTable.h"
#include "utilitiesHypercube.h"
#include "fFunction.h"
#include <math.h>       /* log */
#include "BruteForceImplementation.h"
#include "HashFunctions.h"
int main(int argc,char **argv) {
    string inputFile,queryFile,outputFile;
    int d=-1,M=10,probes=2,k=4;
    assign_values_from_arguments(argc,argv,inputFile,queryFile,d,M,probes,outputFile);
    cout<<inputFile<<" "<<queryFile<<" "<<d<<" "<<M<<" "<<probes<<" "<<outputFile<<endl;

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


        InputGenericVector<int> pointsVector(inputFile);
        double radius;
        InputGenericVector<int> queriesVector(queryFile, radius);

        ExactNeighboursVector<int> neighboursVector(pointsVector, queriesVector);
        double w = neighboursVector.wCalculator();
        if(d==-1)
            d = (int) ceil(log(pointsVector.itemValues.size()));//d'=logn
        CubeHashTable<int> cubeHashTable;

        fFunctions<int> functions;
        vector<HashFunctions<int>> hashFunctions;
        hashFunctions.reserve((unsigned long) d);
        for (unsigned int i = 0; i < d; i++) {
            hashFunctions.emplace_back(HashFunctions<int>(k, pointsVector.itemValues[0].second.size(), w));
        }

        pointsHypercube(pointsVector, cubeHashTable, hashFunctions, (unsigned int) k, w, d, functions);

        queriesHypercube(queriesVector, cubeHashTable, hashFunctions, (unsigned int) k, w, d, functions, outputFile, M,
                         probes, neighboursVector, radius);
        cout << "Do you want to repeat the process with another query file?(Y/N)" << endl;
        string answer;
        while (answer!="N" && answer!="Y")
            getline(cin, answer);
        if (answer == "N")
            break;
        queryFile.erase();
    }
return 5;
}