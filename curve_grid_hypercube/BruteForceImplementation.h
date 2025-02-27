#ifndef LSH_INPUTNEIGHBOURS_H
#define LSH_INPUTNEIGHBOURS_H
#include "InputImplementation.h"
#include <cstdio>
#include <ctime>
#include<tuple> // for tuple
#include <chrono>
template <class inputData>
class ExactNeighboursVector{
private:
    //Vector of pairs. Each pair has an itemID and a tuple with its nearest neighbour id,distance and time to find it.
    vector< pair<string, tuple<string,double,double>> > neighboursVector;
    void save (string const& itemID,string const& neighborID,double const& value,double const& time);
    double manhattanDistance(vector<inputData> const& point, vector<inputData> const& query);
    inline double Dtw(vector<pair<double,double>> const& p,vector<pair<double,double>> const& q);
public:
    ExactNeighboursVector(InputGenericVector<inputData> const& pointsVector,InputGenericVector<inputData> const& queriesVector,bool input);
    ExactNeighboursVector(vector<vector<double>> const& pointsVector,vector<vector<double>> const& queriesVector,bool input);
    void printNeighborsToFile();
    double wCalculator();
    double getRealDistance(int index);
    double getRealTime(int index);
    string getRealNeighbor(int index);
};
#endif //LSH_INPUTNEIGHBOURS_H
