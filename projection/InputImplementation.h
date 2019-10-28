#ifndef LSH_INPUTIMPLEMENTATION_H
#define LSH_INPUTIMPLEMENTATION_H
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <bits/stdc++.h>
#define MAX_POINTS_NUMBER 4
using namespace std;
template <class inputData>
class InputGenericVector{
private:

    vector<inputData> temp;
    string itemID;
    void push(inputData const&);  // push value to temp
    void save(); //save temp to inputValues
    void constructorFunction(InputGenericVector&,string const&);
public:
    //Vector of pairs. Each pair has as first item the itemID and as second item a vector with all the points of each itemID.
    vector< pair<string,vector<inputData> >> itemValues;
    //Function overloading
    explicit InputGenericVector(string path,double& radius); //for queries
    explicit InputGenericVector(string const& path); //for points
    explicit InputGenericVector(string const& path,unsigned int& maxCurveSize,unsigned int& minCurveSize,bool const& input); //for curves
    void printVector();
    void maxCoordFinder(double& max);
};
#endif
