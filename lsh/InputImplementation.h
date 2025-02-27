#ifndef LSH_INPUTIMPLEMENTATION_H
#define LSH_INPUTIMPLEMENTATION_H
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
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
    explicit InputGenericVector(string path,double& radius);
    explicit InputGenericVector(string path);
    void printVector();
};
//GenericVector<int> inputToVector();
#endif
