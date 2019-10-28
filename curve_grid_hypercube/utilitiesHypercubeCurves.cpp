#include "utilitiesHypercubeCurves.h"
#include "InputImplementation.h"
void assign_values_from_arguments(int const& argc,char **argv,string& inputFile,string& queryFile,int& k,int& L,string& outputFile,int& M,int& probes){
    for(int i=1; i<argc; i+=2){
        if(strcmp(argv[i], "-d") == 0)
            inputFile=string(argv[i+1]);
        else if(strcmp(argv[i], "-q") == 0)
            queryFile=string(argv[i+1]);
        else if(strcmp(argv[i], "-k_hypercube") == 0)
            k=stoi(argv[i+1]);
        else if(strcmp(argv[i], "-L_grid") == 0)
            L=stoi(argv[i+1]);
        else if(strcmp(argv[i], "-o") == 0)
            outputFile=string(argv[i+1]);
        else if(strcmp(argv[i], "-M") == 0)
            M=stoi(argv[i+1]);
        else if(strcmp(argv[i], "-probes") == 0)
            probes=stoi(argv[i+1]);
    }
}

double deltaCalculator(InputGenericVector<pair<double,double>> const& point){
double delta=0;
for(unsigned int i=0; i<point.itemValues.size(); i++){
double temp=0;
for(unsigned int coord=0; coord<point.itemValues[i].second.size()-1; coord++){
temp+=sqrt((point.itemValues[i].second[coord].first - point.itemValues[i].second[coord+1].first)*(point.itemValues[i].second[coord].first - point.itemValues[i].second[coord+1].first)
+ (point.itemValues[i].second[coord].second - point.itemValues[i].second[coord+1].second)*(point.itemValues[i].second[coord].second - point.itemValues[i].second[coord+1].second));   //Euclidean distance
}
temp/=point.itemValues[i].second.size();
delta += temp;
}
delta/=point.itemValues.size();
return delta;
}

void curves_points_HyperCube(int const& L_grid,InputGenericVector<pair<double,double>>& pointsVector,vector<vector<vector<double>>> const& concatenated_gridCurve_vectors,int const& d,
        int const& k_hypercube,vector<double> const& w,vector<vector<HashFunctions<double>>>& hashFunctions,fFunctions<double>& f,vector<CubeHashTable<pair<double,double>>>& cubeHashTables){
    for(unsigned int i=0; i<L_grid;i++){
        for(unsigned int curve=0; curve<pointsVector.itemValues.size(); curve++){
            string key=f.binaryStringCalculator(concatenated_gridCurve_vectors[i][curve],(int)d,k_hypercube,w[i],hashFunctions[i]);
            cubeHashTables[i].insertCube(key,pointsVector.itemValues[curve]);
        }
    }
}

void curves_queries_HyperCube(int const& L_grid,InputGenericVector<pair<double,double>>& queriesVector,vector<vector<vector<double>>> const& concatenated_gridCurve_vectors,int const& d,
                             int const& k_hypercube,vector<double> const& w,vector<vector<HashFunctions<double>>>& hashFunctions,fFunctions<double>& f,
                             vector<CubeHashTable<pair<double,double>>>& cubeHashTables,vector<Grid<pair<double,double>>>& grids,int const& M,int const& probes,
                              unsigned int const& maxCurveSize,double const& max_coord,string const& outputFile,ExactNeighboursVector<pair<double,double>>& neighboursVector){
    streambuf *psbuf, *backup;
    ofstream filestr;
    filestr.open (outputFile);

    backup = cout.rdbuf();     // back up cout's streambuf

    psbuf = filestr.rdbuf();        // get file's streambuf
    cout.rdbuf(psbuf);         // assign streambuf to cout
    double total_duration=0;
    double max=-10;
    double averageAT=0;
    for(unsigned int curve=0; curve<queriesVector.itemValues.size(); curve++){
        tuple<string,double,double> neighbor;
        vector<double> concatenated_gridCurve_vector;

        auto start = std::chrono::system_clock::now();
        concatenated_gridCurve_vector=grids[0].concatenated_vector_from_gridCurve(queriesVector.itemValues[curve].second,maxCurveSize,max_coord);
        string key=f.binaryStringCalculator(concatenated_gridCurve_vector,(int)d,k_hypercube,w[0],hashFunctions[0]);
        neighbor = cubeHashTables[0].nearestNeighbor(queriesVector.itemValues[curve], key,M,probes);
        for(unsigned int grid=1; grid<L_grid; grid++){
            //cout<<"grid "<<grid<<endl;
            //Create the concatenated grid curve vector for each curve of each grid
            concatenated_gridCurve_vector=grids[grid].concatenated_vector_from_gridCurve(queriesVector.itemValues[curve].second,maxCurveSize,max_coord);


            key=f.binaryStringCalculator(concatenated_gridCurve_vector,(int)d,k_hypercube,w[grid],hashFunctions[grid]);
            //cout<<"komple"<<endl;
            tuple<string, double, double> temp;
            temp = cubeHashTables[grid].nearestNeighbor(queriesVector.itemValues[curve], key,M,probes);
            if((get<1>(neighbor)>get<1>(temp) && get<0>(temp)!="None") || get<0>(neighbor)=="None")
                neighbor=temp;

        }
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> duration = end-start;
        total_duration+=duration.count();
        if(get<0>(neighbor)!="None" && neighboursVector.getRealDistance(curve)!=0) {
            double tempAf=get<1>(neighbor)/neighboursVector.getRealDistance(curve);
            if (tempAf > max)
                max = tempAf;
            averageAT += get<2>(neighbor);
        }
        cout<<"Query: "<<queriesVector.itemValues[curve].first<<endl;
        cout<<"Method: LSH"<<endl;
        cout<<"HashFunction: Hypercube"<<endl;
        cout<<"Found Nearest neighbor: "<<get<0>(neighbor)<<endl;
        cout<<"True Nearest neighbor: "<<neighboursVector.getRealNeighbor(curve)<<endl;
        cout<<"distanceFound: "<<get<1>(neighbor)<<endl;
        cout<<"distanceTrue: "<<neighboursVector.getRealDistance(curve)<<endl;
    }
    std::cout.rdbuf(backup);        // restore cout's original streambuf
    filestr.close();
    cout<<"Time for approximate neighbors through lsh: "<<total_duration<<endl;
    averageAT=averageAT/queriesVector.itemValues.size();
    cout<<"MaxAF: "<<max<<endl;
    cout<<"AverageAT: "<<averageAT<<endl;
}