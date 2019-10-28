#include "InputImplementation.h"
#include "BruteForceImplementation.h"
#include "HashTables.h"
#include <string.h>
#include "HashFunctions.h"
#include "utilitiesLSH.h"
//#define TYPE_OF_DATA int

int main(int argc, char **argv) {
    string inputFile,queryFile,outputFile;
    int k=4,L=5;
    assign_values_from_arguments(argc,argv,inputFile,queryFile,k,L,outputFile);
    cout<<inputFile<<" "<<queryFile<<" "<<k<<" "<<L<<" "<<outputFile<<endl;
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
        //ExactNeighboursVector<int> inputNeighboursVector(pointsVector,pointsVector, true);

        double w = 4656.37;//inputNeighboursVector.wCalculator();
        ExactNeighboursVector<int> neighboursVector(pointsVector, queriesVector, false);
        //neighboursVector.printNeighborsToFile();
        auto tableSize = (unsigned int) pointsVector.itemValues.size() / 8;
        HashTables<int> hashTables((unsigned int) L);


        vector<HashFunctions<int>> hashFunctions;
        hashFunctions.reserve((unsigned long) L);
        for (unsigned int i = 0; i < L; i++) { //Create L g functions
            hashFunctions.emplace_back(HashFunctions<int>(k, pointsVector.itemValues[0].second.size(), w));
        }


        pointsLSH(pointsVector, hashTables, L, hashFunctions, (unsigned int) k, w);


        queriesLSH(queriesVector, hashTables, hashFunctions, (unsigned int) k, w, neighboursVector, radius, outputFile);

        cout<<"Do you want to repeat the process with another query file?(Y/N)"<<endl;
        string answer;
        while (answer!="N" && answer!="Y")
            getline(cin, answer);
        if(answer=="N")
            break;
        queryFile.erase();
    }
    return 5;
}