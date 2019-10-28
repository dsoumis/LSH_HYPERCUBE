#include "HashTables.h"
template <class inputData>
void HashTables<inputData>::insertHashtable(unsigned int whichHashTable, unsigned int const& g,pair<string,vector<inputData >>& point){
    hashTables[whichHashTable].insert(make_pair(g, make_tuple(point.first,point.second)));
}

template <>
double HashTables<int>::manhattanDistance(vector<int> const& point, vector<int> const& query){
    double distance=0;
    for(unsigned int i=0; i<point.size(); i++){
        distance+=abs(point[i] - query[i]);
    }
    return distance;
}

inline double minf(const double a,const double b,const double c) {
    double min;
    min=a;
    if(b<min) min=b;
    if(c<min) min=c;
    return min;
}
inline double Dtw(vector<pair<double,double>> &p,vector<pair<double,double>> &q){
    double distance;
    double C[p.size()][q.size()];
    for(int i=0; i<p.size(); i++) {
        for(int j = 0; j<q.size(); j++) {
            if(i==0 && j==0) {
                C[0][0] = sqrt((q[j].first - p[i].first)*(q[j].first - p[i].first) + (q[j].second - p[i].second)*(q[j].second - p[i].second));             //initialize first value with ||p-q||
                continue;
            }
            distance = sqrt((q[j].first - p[i].first)*(q[j].first - p[i].first) + (q[j].second - p[i].second)*(q[j].second - p[i].second));   //Euclidean distance

            if (i > 0 && j > 0) {
                C[i][j] = minf(C[i - 1][j], C[i - 1][j - 1], C[i][j-1]) + distance;
            } else if (i > 0 && j==0) {
                C[i][0] = C[i-1][0] + distance;
            } else if (j > 0 && i==0) {
                C[0][j] =  C[0][j-1] + distance;
            }
        }
    }
    return C[p.size()-1][q.size()-1];
}







template <>
tuple<string,double,double> HashTables<int>::nearestNeighbor (pair<string,vector<int >>& query,
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
template <>
tuple<string,double,double> HashTables<pair<double,double>>::nearestNeighbor (pair<string,vector<pair<double,double> >>& query,unsigned int const& g) {
    auto start = std::chrono::system_clock::now();
    tuple<string,double,double> b; //Nearest Neighbor
    double db=INT32_MAX; //Min distance
    bool bucketEmpty= true;
    //cout<<"giauto to grid: "<<hashTables[0].count(g)<<endl;
    auto iterator=hashTables[0].equal_range(g); //iterate all values with key g
    for(auto p = iterator.first; p != iterator.second; ++p) { //for each item p in bucket gi(q) do
//        if(hashTables[0].count(g)>150)
//            break;
        bucketEmpty = false;
        //if dist(q,p)<db then b←p; db←dist(q,p)
        double dist = Dtw(get<1>(p->second), query.second);
        if (dist < db) {
            get<0>(b) = get<0>(p->second);
            db = dist;
            get<1>(b) = dist;
        }
    }
    //exit(1);
    if(bucketEmpty){
        //cout<<"brethike"<<endl;
        get<0>(b)="None";
        get<1>(b)=0;
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> duration = end-start;
    get<2>(b)=duration.count();
    //cout<<"geitonas "<<get<0>(b)<<" apostasi "<<get<1>(b)<<" xronos "<<duration.count()<<endl;
    return b;
}
template <>
vector<string> HashTables<int>::rangeSearch (pair<string,vector<int >>& query,double radius,
                                                   vector<HashFunctions<int>> hashFunctions,unsigned int const& k,double const& w) {
    vector<string> b; //Range Neighbors
    bool allBucketsEmpty= true;
    for(unsigned int i=0; i<hashTables.size(); i++){ //for i from 1 to L do
        unsigned int g=hashFunctions[i].gCalculator(query.second,k,w);
        for(auto p=hashTables[i].equal_range(g).first; p!=hashTables[i].equal_range(g).second; p++) { //for each item p in bucket gi(q) do
            if (hashTables[i].count(g) > 100 * hashTables.size()) //if large number of retrieved items (e.g.>100L) then Break
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
template class HashTables<pair<double,double>>; //In order to not fail the compile as the compiler wants to see the data that the templated class will have.