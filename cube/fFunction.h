#ifndef HYPERCUBE_FFUNCTION_H
#define HYPERCUBE_FFUNCTION_H

#include "HashFunctions.h"
#include <unordered_map>
template <class inputData>
class fFunctions{
private:
    //An unordered map of values which were already encountered
    unordered_map<int,string> alreadySeenValues;
public:
    string binaryStringCalculator(vector<inputData> const& p,int d,int const& k,double const& w,vector<HashFunctions<inputData>> hashFunctions);
};
#endif //HYPERCUBE_FFUNCTION_H
