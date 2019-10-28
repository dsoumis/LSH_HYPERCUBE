#include "HashTables.h"
template <class inputData>
void HashTables<inputData>::insertHashtable(unsigned int whichHashTable, unsigned int const& g,pair<string,vector<inputData >>& point){
    hashTables[whichHashTable].insert(make_pair(g, make_tuple(point.first,point.second)));
}

template <class inputData>
inline double HashTables<inputData>::manhattanDistance(vector<inputData> const& point, vector<inputData> const& query){
    double distance=0;
    for(unsigned int i=0; i<point.size(); i++){
        distance+=abs(point[i] - query[i]);
    }
    return distance;
}

template <class inputData>
inline tuple<string,double,double> HashTables<inputData>::nearestNeighbor (pair<string,vector<inputData >>& query,
        vector<HashFunctions<int>> hashFunctions,unsigned int const& k,double const& w) {
    auto start = std::chrono::system_clock::now();
    tuple<string,double,double> b; //Nearest Neighbor
    double db=INT32_MAX; //Min distance
    bool allBucketsEmpty= true;
    int L=(int)hashTables.size();
    for(unsigned int i=0; i<L; i++){ //for i from 1 to L do
        unsigned int g=hashFunctions[i].gCalculator(query.second,k,w);
        //cout<<"to bucket "<<bucket<<endl;
        for(auto p=hashTables[i].equal_range(g).first; p!=hashTables[i].equal_range(g).second; p++) { //for each item p in bucket gi(q) do
            if (hashTables[i].count(g) > 100 * L) //if large number of retrieved items (e.g.>10L) then Break
                break;

            allBucketsEmpty = false;
            //if dist(q,p)<db then b←p; db←dist(q,p)
            double dist = manhattanDistance(get<1>(p->second), query.second);
            if (dist < db) {
                get<0>(b) = get<0>(p->second);
                db = dist;
                get<1>(b) = dist;
            }

        }
    }
    if(allBucketsEmpty){
        get<0>(b)="None";
        get<1>(b)=0;
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> duration = end-start;
    get<2>(b)=duration.count();
    //cout<<"geitonas "<<get<0>(b)<<" apostasi "<<get<1>(b)<<" xronos "<<duration.count()<<endl;
    return b;
}
template <class inputData>
vector<string> HashTables<inputData>::rangeSearch (pair<string,vector<inputData >>& query,double radius,
                                                   vector<HashFunctions<int>> hashFunctions,unsigned int const& k,double const& w) {
    vector<string> b; //Range Neighbors
    bool allBucketsEmpty= true;
    for(unsigned int i=0; i<hashTables.size(); i++){ //for i from 1 to L do
        unsigned int g=hashFunctions[i].gCalculator(query.second,k,w);
        for(auto p=hashTables[i].equal_range(g).first; p!=hashTables[i].equal_range(g).second; p++) { //for each item p in bucket gi(q) do
            if (hashTables[i].count(g) >100 * hashTables.size()) //if large number of retrieved items (e.g.>10L) then Break
                break;

            allBucketsEmpty = false;
            //if dist(q,p) is within radius save the range neighbor
            double dist = manhattanDistance(get<1>(p->second), query.second);
            if (dist <= radius) {
                b.push_back(get<0>(p->second));
            }

        }
    }
    if(allBucketsEmpty){
        b.emplace_back("None");
    }
    return b;
}
template class HashTables<int>; //In order to not fail the compile as the compiler wants to see the data that the templated class will have.