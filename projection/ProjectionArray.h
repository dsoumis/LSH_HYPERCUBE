
#ifndef CURVES_PROJECTIONARRAY_H
#define CURVES_PROJECTIONARRAY_H

#include <iostream>
#include <vector>
#include <bits/stdc++.h>
#include <string>
#include <unordered_map>
using namespace std;

typedef pair<int, int> pairs;

class cell{
public:
    vector<vector<pairs>> RT;
    unordered_map<string, vector<pair<int,vector<double>>>> imap;        //key: "i-j-traversal"
    unordered_map<string, vector<pair<string,vector<double>>>> jmap;        //key: "i-j-traversal"
};

class ProjectionArray {
public:
    vector<vector<cell>>array;      //2-d vector
     ProjectionArray(unsigned long n);
     void  calculateArray(unsigned long size);
     void calculatePoints(unsigned long m,  unsigned long n,vector<vector<pairs>> &RT);
     void calculateAllPaths(int m, int n,  const set<pairs>& s, vector<vector<pairs>> &RT );
     void calculateAllPathsUtil(int i, int j, int m, int n, vector<pairs> &path, int pi,  const set<pairs>& s, vector<vector<pairs>> &RT);
};

void calculateG( vector<double> &G, int d, double e);
double generateNumberG(double const& range_from,double const& range_to);
void displayPairs(const set<pairs>& s);
void createVector(vector<pair<double,double>> const &curve, vector<double> &G, vector<pairs> &RT , string &curveId, int ipos, int jpos, cell & cell, int traversalpos, int flag);
void calculateVectors(vector<pair<double,double>> const &curve, vector<double> &G,cell &traversals, string &curveId, int ipos, int jpos,int flag);
#endif //CURVES_PROJECTIONARRAY_H
