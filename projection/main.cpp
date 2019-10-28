#include "InputImplementation.h"
#include "ProjectionArray.h"
#include "HashFunctions.h"
#include "HashTables.h"
#define d 2
#define e 0.5

void assign_values_from_arguments(int const& argc,char **argv,string& inputFile,string& queryFile,int& k,int& L,string& outputFile){
    for(int i=1; i<argc; i+=2){
        if(strcmp(argv[i], "-d") == 0)
            inputFile=string(argv[i+1]);
        else if(strcmp(argv[i], "-q") == 0)
            queryFile=string(argv[i+1]);
        else if(strcmp(argv[i], "-k_vec") == 0)
            k=stoi(argv[i+1]);
        else if(strcmp(argv[i], "-L_vec") == 0)
            L=stoi(argv[i+1]);
        else if(strcmp(argv[i], "-o") == 0)
            outputFile=string(argv[i+1]);
    }
}
int main(int argc, char** argv) {
    string inputFile,queryFile,outputFile;
    int k_vec=4,L_vec=4;
    assign_values_from_arguments(argc,argv,inputFile,queryFile,k_vec,L_vec,outputFile);
    int maxCurveSize=637;
    unsigned int maxPointsSize=0,maxQueriesSize=0,minPointsSize=4294967295,minQueriesSize=4294967295;
    InputGenericVector<pair<double,double>> pointsVector(inputFile, maxPointsSize,minPointsSize, true);

    InputGenericVector<pair<double,double>> queriesVector(inputFile, maxQueriesSize, minQueriesSize, false);


    ProjectionArray A(maxPointsSize);
    A.calculateArray(maxPointsSize);
    cout<<"Number of input curves "<<pointsVector.itemValues.size()<<endl;
    cout<<"Number of queries curves "<<queriesVector.itemValues.size()<<endl<<endl;


    vector<double> G;
    calculateG(G, d, e);


    for(int i=0; i<(int)pointsVector.itemValues.size(); i++) { // for points
        int k = pointsVector.itemValues[i].second.size()-1;
        for (int j = 0; j < (int) maxPointsSize; j++) {
            calculateVectors(pointsVector.itemValues[i].second, G, A.array[k][j], pointsVector.itemValues[i].first,i, j, 0); //flag 0 to calculate points
        }
    }
    cout<<"Vectors for points created."<<endl;

    for(int j=0; j<(int)queriesVector.itemValues.size();j++){   //for queries
        int k = queriesVector.itemValues[j].second.size()-1;
        for(int i=0; i<(int)maxQueriesSize; i++){
            calculateVectors(queriesVector.itemValues[j].second, G, A.array[i][k], queriesVector.itemValues[j].first,i,j,1); //flag 1 to calculate queries
        }
    }
    cout<<"Vectors for queries created."<<endl<<endl;

    int sum = 0;
    for(int i=0; i<maxPointsSize; i++){
        for(int j=0; j<maxPointsSize; j++){
            sum +=  A.array[i][j].imap.size();
        }
    }

    vector<HashFunctions<double>> hashFunctions;
    hashFunctions.reserve((unsigned long) sum);
    vector<HashTables<pair<double, double>>> hashTables;
    hashTables.reserve((unsigned long) sum);

    vector<double> w; //One w for each traversal
    w.resize((unsigned long) sum);

    for (unsigned int i = 0; i < sum; i++) {
//        //Exact Neighbours here calculates only 5% of neighbors for speed
//        ExactNeighboursVector<double> temp_concatenatedNeighbors(concatenated_gridCurve_vectors[i],
//                                                                 concatenated_gridCurve_vectors[i], true);
        w[i] = 50000;//temp_concatenatedNeighbors.wCalculator();
        //cout<<w[i]<<endl;

        hashFunctions.emplace_back(HashFunctions<double>(k_vec, maxCurveSize, w[i]));//Create sum g functions. One for each sum.
        hashTables.emplace_back(HashTables<pair<double, double>>(1)); //Create L_grid LSHs with one hash table each.
    }

    unsigned int no_traversal=0;
    for(int i=0; i<maxPointsSize;i++){
        for(int j=0; j<maxPointsSize; j++ ){    //for every i,j in n*n array

            for (auto& it: A.array[i][j].imap) {

                for (auto& vec: it.second){     //for every vector in imap
                    unsigned int g=hashFunctions[no_traversal].gCalculator(vec.second,(unsigned int)k_vec,w[no_traversal]);
                    hashTables[no_traversal].insertHashtable(0,g,pointsVector.itemValues[vec.first]);
                    no_traversal++;
                }


            }
        }
    }

    cout<< "Input curves inserted with LSH."<<endl<<endl;

    no_traversal=0;
    for(unsigned int curve=0; curve<queriesVector.itemValues.size(); curve++){      //for every query curve
        int k = (int)queriesVector.itemValues[curve].second.size()-1;               //k =size of curve -1
        tuple<string,double,double> neighbor;
        int counter =0;
        for(int i=0; i<maxQueriesSize; i++){                //for every i in row k

            for(int traversal=0; traversal<(int)A.array[i][k].RT.size(); traversal++){  //for every RT
                string key =to_string(i)+"-"+to_string(k)+"-"+to_string(traversal);

                for(auto& vec: A.array[i][k].jmap[key]) { //for every i-j-traversal found in jmap

                    if(vec.first == queriesVector.itemValues[curve].first) {    // check if its the same curve
                        unsigned int g = hashFunctions[no_traversal].gCalculator(vec.second, (unsigned int) k_vec, w[no_traversal]);
                        tuple<string, double, double> temp;
                        if (counter == 0)
                            neighbor = hashTables[no_traversal].nearestNeighbor(queriesVector.itemValues[no_traversal],g);
                        else
                            temp = hashTables[no_traversal].nearestNeighbor(queriesVector.itemValues[no_traversal], g);
                        no_traversal++;
                        counter++;
                        if ((get<1>(neighbor) > get<1>(temp) && get<0>(temp) != "None") || get<0>(neighbor) =="None") //If the distance is less in another hashtable but it is not from fault
                            neighbor = temp;
                    }
                }
            }


        }
        cout<<"Query: "<<queriesVector.itemValues[curve].first<<endl;
        cout<<"Method: Projection"<<endl;
        cout<<"HashFunction: LSH"<<endl;
        cout<<"Found Nearest neighbor: "<<get<0>(neighbor)<<endl;
        //cout<<"True Nearest neighbor: "<<neighboursVector.getRealNeighbor(curve)<<endl;
        cout<<"distanceFound: "<<get<1>(neighbor)<<endl;
        //cout<<"distanceTrue: "<<neighboursVector.getRealDistance(curve)<<endl;
    }


    return 0;
}


