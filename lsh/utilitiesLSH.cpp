#include "utilitiesLSH.h"
void assign_values_from_arguments(int const& argc,char **argv,string& inputFile,string& queryFile,int& k,int& L,string& outputFile){
    for(int i=1; i<argc; i+=2){
        if(strcmp(argv[i], "-d") == 0)
            inputFile=string(argv[i+1]);
        else if(strcmp(argv[i], "-q") == 0)
            queryFile=string(argv[i+1]);
        else if(strcmp(argv[i], "-k") == 0)
            k=stoi(argv[i+1]);
        else if(strcmp(argv[i], "-L") == 0)
            L=stoi(argv[i+1]);
        else if(strcmp(argv[i], "-o") == 0)
            outputFile=string(argv[i+1]);
    }
}

void pointsLSH(InputGenericVector<int>& pointsVector,HashTables<int>& hashTables,int const& L,vector<HashFunctions<int>>& hashFunctions,unsigned int const& k, double const& w){
    for(unsigned int i=0; i<pointsVector.itemValues.size(); i++){
        for(unsigned int numberHashtable=0; numberHashtable<L; numberHashtable++){
            unsigned int g=hashFunctions[numberHashtable].gCalculator(pointsVector.itemValues[i].second,k,w);
            hashTables.insertHashtable(numberHashtable,g,pointsVector.itemValues[i]);
        }
    }
}





void queriesLSH(InputGenericVector<int>& queriesVector,HashTables<int>& hashTables,vector<HashFunctions<int>> const& hashFunctions,unsigned int const& k, double const& w,
        ExactNeighboursVector<int>& neighboursVector,double const& radius,string const& outputFile){

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
        neighbor=hashTables.nearestNeighbor(queriesVector.itemValues[i],hashFunctions,k,w);
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> duration = end-start;
        total_duration+=duration.count();

        vector<string> rangeNeighbors;
        if(radius!=0) {
            rangeNeighbors = hashTables.rangeSearch(queriesVector.itemValues[i], radius, hashFunctions,(unsigned int) k, w);
        }

        if(get<0>(neighbor)!="None"&& neighboursVector.getRealDistance(i)!=0) {
            double tempAf=get<1>(neighbor)/neighboursVector.getRealDistance(i);
            if (tempAf > max)
                max = tempAf;
            averageAT += get<2>(neighbor);
        }
        cout<<"Query: "<<queriesVector.itemValues[i].first<<endl;
        cout<<"Nearest neighbor: "<<get<0>(neighbor)<<endl;
        cout<<"distanceLSH: "<<get<1>(neighbor)<<endl;
        cout<<"distanceTrue: "<<neighboursVector.getRealDistance(i)<<endl;
        cout<<"tLSH: "<<get<2>(neighbor)<<endl;
        cout<<"tTrue: "<<neighboursVector.getRealTime(i)<<endl;
        if(radius!=0) {
            cout<<"R-near neighbors:"<<endl;
            for (unsigned int rangeNeighbor = 0; rangeNeighbor < rangeNeighbors.size(); rangeNeighbor++)
                cout << rangeNeighbors[rangeNeighbor] << endl;
        }
    }

    std::cout.rdbuf(backup);        // restore cout's original streambuf
    filestr.close();
    cout<<"Time for approximate neighbors through lsh: "<<total_duration<<endl;
    averageAT=averageAT/queriesVector.itemValues.size();
    cout<<"MaxAF: "<<max<<endl;
    cout<<"AverageAT: "<<averageAT<<endl;
}