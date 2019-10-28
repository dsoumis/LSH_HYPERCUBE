#include "CubeHashTable.h"
template <>
inline double CubeHashTable<int>::manhattanDistance(vector<int> const& point, vector<int> const& query){
    unsigned int distance=0;
    for(unsigned int i=0; i<point.size(); i++){
        distance+=abs(point[i] - query[i]);
    }
    return distance;
}
template <class inputData>
void CubeHashTable<inputData>::insertCube(string key, pair<string, vector<inputData>> &point){
    hashTable.insert(make_pair(key, make_pair(point.first,point.second)));
}


inline void all_combinations_of_strings_with_hammingDistance(string str,int strLen,int hammingDistance,vector<string>& keys) {
    if (hammingDistance == 0) {
        keys.push_back(str);
        return;
    }
    if (strLen < 0) return;
    // flip current bit
    str[strLen] = str[strLen] == '0' ? '1' : '0';
    all_combinations_of_strings_with_hammingDistance(str, strLen-1, hammingDistance-1,keys);
    // or don't flip it (flip it again to undo)
    str[strLen] = str[strLen] == '0' ? '1' : '0';
    all_combinations_of_strings_with_hammingDistance(str, strLen-1, hammingDistance,keys);
}
template <class inputData>
inline tuple<string,double,double> CubeHashTable<inputData>::nearestNeighbor (pair<string,vector<inputData >>& query,string const& key,int const& M,int const& probes) {
    vector<string> keys;
    keys.push_back(key);
    //In order to check all near vertices
    //Calculate vertices with hamming distance from 1 to x until vertices with distance h are equal or more to probes
    int hamming=1;
    while(probes>keys.size()) {
        all_combinations_of_strings_with_hammingDistance(key, (int) key.length() - 1, hamming, keys);
        hamming++;
    }
    auto start = std::chrono::system_clock::now();
    tuple<string,double,double> b; //Nearest Neighbor
    double db=INT32_MAX; //Min distance
    int vertices=0;
    for(unsigned int i=0; i<keys.size(); i++) {
        int points_checked=0;
        vertices++;
        if(vertices>probes)
            break;
        for (auto p = hashTable.equal_range(keys[i]).first; p != hashTable.equal_range(keys[i]).second; p++) { //for each item p in vertex with key
            points_checked++;
            if(M<points_checked)
                break;
            //if dist(q,p)<db then b←p; db←dist(q,p)
            double dist = manhattanDistance(p->second.second, query.second);
            if (dist < db) {
                get<0>(b) = p->second.first;
                db = dist;
                get<1>(b) = dist;
            }

        }
    }


    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> duration = end-start;
    get<2>(b)=duration.count();
    //cout<<"geitonas "<<get<0>(b)<<" apostasi "<<get<1>(b)<<" xronos "<<duration.count()<<endl;
    return b;
}


template <class inputData>
inline vector<string> CubeHashTable<inputData>::rangeSearch (pair<string,vector<inputData >>& query,string const& key,double radius) {
    vector<string> b; //Range Neighbors
    for(auto p=hashTable.equal_range(key).first; p!=hashTable.equal_range(key).second; p++) { //for each item p in bucket gi(q) do

        //if dist(q,p) is within radius save the range neighbor
        double dist = manhattanDistance(p->second.second, query.second);
        if (dist <= radius) {
            b.push_back(p->second.first);
        }

    }

    return b;
}
template class CubeHashTable<int>; //In order to not fail the compile as the compiler wants to see the data that the templated class will have.