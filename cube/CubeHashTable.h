#ifndef HYPERCUBE_CUBEHASHTABLE_H
#define HYPERCUBE_CUBEHASHTABLE_H
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <cstdio>
#include <ctime>
#include <chrono>
using namespace std;
template <class inputData>
class CubeHashTable{
private:
    //Unordered multimap allows many values to 1 key,which means collisions in 1 bucket.
    //A bucket has an binary string key and a pair of itemID, points of item
    unordered_multimap<string, pair<string,vector<inputData>> > hashTable;
    double manhattanDistance(vector<inputData> const& point, vector<inputData> const& query);
public:
    void insertCube(string key,pair<string,vector<inputData >>& point);
    tuple<string,double,double> nearestNeighbor(pair<string,vector<inputData >>& query,string const& key,int const& M,int const& probes);
    vector<string> rangeSearch (pair<string,vector<inputData >>& query,string const& key,double radius);
};
#endif //HYPERCUBE_CUBEHASHTABLE_H
