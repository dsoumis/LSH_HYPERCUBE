#include "utilitiesHypercube.h"
void assign_values_from_arguments(int const& argc,char **argv,string& inputFile,string& queryFile,int& k,int& M,int& probes,string& outputFile){
    for(int i=1; i<argc; i+=2){
        if(strcmp(argv[i], "-d") == 0)
            inputFile=string(argv[i+1]);
        else if(strcmp(argv[i], "-q") == 0)
            queryFile=string(argv[i+1]);
        else if(strcmp(argv[i], "-k") == 0)
            k=stoi(argv[i+1]);
        else if(strcmp(argv[i], "-M") == 0)
            M=stoi(argv[i+1]);
        else if(strcmp(argv[i], "-probes") == 0)
            probes=stoi(argv[i+1]);
        else if(strcmp(argv[i], "-o") == 0)
            outputFile=string(argv[i+1]);
    }
}

void pointsHypercube(InputGenericVector<int>& pointsVector,CubeHashTable<int>& cubeHashTable,vector<HashFunctions<int>>& hashFunctions,unsigned int const& k, double const& w,int const& d,
               fFunctions<int>& functions){
    for(unsigned int i=0; i<pointsVector.itemValues.size(); i++){
        string key=functions.binaryStringCalculator(pointsVector.itemValues[i].second,d,k,w,hashFunctions);
        cubeHashTable.insertCube(key,pointsVector.itemValues[i]);
    }
}

void queriesHypercube(InputGenericVector<int>& queriesVector,CubeHashTable<int>& cubeHashTable,vector<HashFunctions<int>>& hashFunctions,unsigned int const& k, double const& w,int const& d,
                     fFunctions<int>& functions,string const& outputFile,int const& M,int const& probes,ExactNeighboursVector<int>& neighboursVector,
                      double const& radius){
    streambuf *psbuf, *backup;
    ofstream filestr;
    filestr.open (outputFile);

    backup = cout.rdbuf();     // back up cout's streambuf

    psbuf = filestr.rdbuf();        // get file's streambuf
    cout.rdbuf(psbuf);         // assign streambuf to cout
    double total_duration=0;
    double max=-10;
    double averageAT=0;
    for(unsigned int i=0; i<queriesVector.itemValues.size(); i++){
        tuple<string,double,double> neighbor;
        auto start = std::chrono::system_clock::now();

        string key=functions.binaryStringCalculator(queriesVector.itemValues[i].second,d,k,w,hashFunctions);
        neighbor=cubeHashTable.nearestNeighbor(queriesVector.itemValues[i],key,M,probes);

        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> duration = end-start;
        total_duration+=duration.count();

        vector<string> rangeNeighbors;
        if(radius!=0) {
            rangeNeighbors = cubeHashTable.rangeSearch(queriesVector.itemValues[i], key, radius);
        }


        if(get<0>(neighbor)!="None"&& neighboursVector.getRealDistance(i)!=0) {
            double tempAf=get<1>(neighbor)/neighboursVector.getRealDistance(i);
            if (tempAf > max)
                max = tempAf;
            averageAT += get<2>(neighbor);
        }


        cout<<"Query: "<<queriesVector.itemValues[i].first<<endl;
        cout<<"Nearest neighbor: "<<get<0>(neighbor)<<endl;
        cout<<"distanceCube: "<<get<1>(neighbor)<<endl;
        cout<<"distanceTrue: "<<neighboursVector.getRealDistance(i)<<endl;
        cout<<"tCube: "<<get<2>(neighbor)<<endl;
        cout<<"tTrue: "<<neighboursVector.getRealTime(i)<<endl;

        if(radius!=0) {
            cout<<"R-near neighbors:"<<endl;
            for(unsigned int rangeNeighbor=0; rangeNeighbor<rangeNeighbors.size(); rangeNeighbor++)
                cout<<rangeNeighbors[rangeNeighbor]<<endl;
        }
    }
    std::cout.rdbuf(backup);        // restore cout's original streambuf
    filestr.close();
    cout<<"Time for approximate neighbors through lsh: "<<total_duration<<endl;
    averageAT=averageAT/queriesVector.itemValues.size();
    cout<<"MaxAF: "<<max<<endl; //Best 2.7
    cout<<"AverageAT: "<<averageAT<<endl;
}