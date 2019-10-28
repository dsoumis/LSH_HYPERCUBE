#include "BruteForceImplementation.h"
#include <math.h>

template <class inputData>
void ExactNeighboursVector<inputData>::save (string const& itemID,string const& neighborID,double const& value,double const& time) {
    neighboursVector.push_back(make_pair(itemID,make_tuple(neighborID,value,time)));
}


template <>
double ExactNeighboursVector<double>::manhattanDistance(vector<double> const& point, vector<double> const& query){
    double distance=0;
    for(unsigned int i=0; i<point.size(); i++){
        distance+=abs(point[i] - query[i]);
    }
    return distance;
}
template <>
double ExactNeighboursVector<int>::manhattanDistance(vector<int> const& point, vector<int> const& query){
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
template <>
inline double ExactNeighboursVector<pair<double,double>>::Dtw(vector<pair<double,double>> const& p,vector<pair<double,double>> const& q){
    double distance;
    double C[p.size()][q.size()];
    for(int i=0; i<p.size(); i++) {
        for(int j = 0; j<q.size(); j++) {
            if(i==0 && j==0) {
                C[0][0] = sqrt(((q[j].first - p[i].first)*(q[j].first - p[i].first)) + ((q[j].second - p[i].second)*(q[j].second - p[i].second)));             //initialize first value with ||p-q||
                continue;
            }
            distance = sqrt(((q[j].first - p[i].first)*(q[j].first - p[i].first)) + ((q[j].second - p[i].second)*(q[j].second - p[i].second)));   //Euclidean distance

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
ExactNeighboursVector<int>::ExactNeighboursVector(InputGenericVector<int> const& pointsVector,InputGenericVector<int> const& queriesVector,bool input){
    auto total_start = std::chrono::system_clock::now();
    for (unsigned int query=0; query<queriesVector.itemValues.size(); query++){
        auto start = std::chrono::system_clock::now();
        double min=INT32_MAX;
        string neighborID;
        double neighborDistance=0;
        for(unsigned int point=0; point<pointsVector.itemValues.size(); point++){
            if(point==query && input)
                continue;
            double distance=manhattanDistance(pointsVector.itemValues[point].second,queriesVector.itemValues[query].second);

            //cout << "apostasi" << distance << endl;
            if (min > distance) {
                min = distance;
                neighborID=pointsVector.itemValues[point].first;
                neighborDistance=distance;
            }
        }
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> duration = end-start;
        this->save(queriesVector.itemValues[query].first,neighborID,neighborDistance,duration.count());
    }
    auto total_end = std::chrono::system_clock::now();
    std::chrono::duration<double> total_duration = total_end-total_start;
    cout<<"Time for exact neighbors through brute force: "<<total_duration.count()<<endl;
}
template <>
ExactNeighboursVector<double>::ExactNeighboursVector(vector<vector<double>> const& pointsVector,vector<vector<double>> const& queriesVector,bool input){ //Only 5% of neighbors

    for (unsigned int query=0; query<queriesVector.size()*0.05; query++){
        double min=INT32_MAX;
        double neighborDistance=0;
        for(unsigned int point=0; point<pointsVector.size()*0.05; point++){
            if(point==query && input)
                continue;
            double distance=manhattanDistance(pointsVector[point],queriesVector[query]);

            //cout << "apostasi" << distance << endl;
            if (min > distance) {
                min = distance;
                neighborDistance=distance;
            }
        }
        this->save("temp_query","temp",neighborDistance,0);
    }
}
template <>
ExactNeighboursVector<pair<double,double>>::ExactNeighboursVector(InputGenericVector<pair<double,double>> const& pointsVector,InputGenericVector<pair<double,double>> const& queriesVector,bool input){
    auto total_start = std::chrono::system_clock::now();
    for (unsigned int query=0; query<queriesVector.itemValues.size(); query++){
        cout<<"neo"<<query<<endl;
        auto start = std::chrono::system_clock::now();
        auto min=DBL_MAX;
        string neighborID;
        double neighborDistance=0;
        for(unsigned int point=0; point<pointsVector.itemValues.size(); point++){
            cout<<"neo"<<point<<endl;
            if(point==query && input)
                continue;
            double distance=Dtw(pointsVector.itemValues[point].second,queriesVector.itemValues[query].second);

            //cout << "apostasi" << distance << endl;
            if (min > distance) {
                min = distance;
                neighborID=pointsVector.itemValues[point].first;
                neighborDistance=distance;
            }
        }
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> duration = end-start;
        this->save(queriesVector.itemValues[query].first,neighborID,neighborDistance,duration.count());
    }
    auto total_end = std::chrono::system_clock::now();
    std::chrono::duration<double> total_duration = total_end-total_start;
    cout<<"Time for exact neighbors through brute force: "<<total_duration.count()<<endl;
}

template <class inputData>
void ExactNeighboursVector<inputData>::printNeighborsToFile(){
    streambuf *psbuf, *backup;
    ofstream filestr;
    filestr.open ("C:\\Users\\User\\Desktop\\neighbors.txt");

    backup = cout.rdbuf();     // back up cout's streambuf

    psbuf = filestr.rdbuf();        // get file's streambuf
    cout.rdbuf(psbuf);         // assign streambuf to cout
    cout<<neighboursVector.size()<<endl;
    for (unsigned int i=0; i<neighboursVector.size(); i++){
      cout<<"QueryID: "<<neighboursVector[i].first<<"    Neighbor: (Id: "<<get<0>(neighboursVector[i].second)<<",Distance: "<<get<1>(neighboursVector[i].second)<<",Time: "<<get<2>(neighboursVector[i].second)<<")"<<endl;
    }
    std::cout.rdbuf(backup);        // restore cout's original streambuf

    filestr.close();
}
template <class inputData>
double ExactNeighboursVector<inputData>::wCalculator(){
    double w=0;
    for (unsigned int i=0; i<neighboursVector.size(); i++){
        w+=get<1>(neighboursVector[i].second);
    }
    w=w/(neighboursVector.size());
    w=w*4;
    return w;
}
template <class inputData>
string ExactNeighboursVector<inputData>::getRealNeighbor(int index){
    return get<0>(neighboursVector[index].second);
}
template <class inputData>
double ExactNeighboursVector<inputData>::getRealDistance(int index){
    return get<1>(neighboursVector[index].second);
}
template <class inputData>
double ExactNeighboursVector<inputData>::getRealTime(int index){
    return get<2>(neighboursVector[index].second);
}
template class InputGenericVector<int>; //In order to not fail the compile as the compiler wants to see the data that the templated class will have.
template class ExactNeighboursVector<int>;

template class InputGenericVector<pair<double,double>>;
template class ExactNeighboursVector<pair<double,double>>;

template class InputGenericVector<double>;
template class ExactNeighboursVector<double>;