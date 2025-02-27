#include "InputImplementation.h"
template <class inputData>
void InputGenericVector<inputData>::push (inputData const& value) {
    temp.push_back(value);
}
template <class inputData>
void InputGenericVector<inputData>::save () {
    itemValues.push_back(make_pair(itemID,temp));
    temp.clear();
}
template <>
void InputGenericVector<int>::printVector(){
    for (unsigned int i=0; i<itemValues.size(); i++){
        //cout<<itemValues[i].second.size()<<endl;
        cout<<itemValues[i].first<<endl;
        for(unsigned int j=0; j<itemValues[i].second.size(); j++){
            cout<<itemValues[i].second.at(j)<<" ";
        }
        cout<<endl;
    }
}
template <>
void InputGenericVector<pair<double,double>>::printVector(){
    for (unsigned int i=0; i<itemValues.size(); i++){
        //cout<<itemValues[i].second.size()<<endl;
        cout<<itemValues[i].first<<endl;
        cout.precision(numeric_limits< double >::max_digits10);
        for(unsigned int j=0; j<itemValues[i].second.size(); j++){
            cout<<"("<<itemValues[i].second.at(j).first<<","<<itemValues[i].second.at(j).second<<") ";
        }
        cout<<endl;
    }
}
template <>
void InputGenericVector<pair<double,double>>::maxCoordFinder(double& max){
    for (unsigned int i=0; i<itemValues.size(); i++){
        for(unsigned int j=0; j<itemValues[i].second.size(); j++){
            if(max<itemValues[i].second.at(j).first)
                max=itemValues[i].second.at(j).first;
            if(max<itemValues[i].second.at(j).second)
                max=itemValues[i].second.at(j).second;
        }
    }
}
template <>
void InputGenericVector<int>::constructorFunction(InputGenericVector<int>& vector,string const& value){
    vector.push(stoi(value));
}
template <>
void InputGenericVector<pair<double,double>>::constructorFunction(InputGenericVector<pair<double,double>>& vector,string const& value){
    stringstream stream(value);
    pair<double ,double > vector_value;
    string temp;
    getline(stream, temp, ' ');
    vector_value.first=stod(temp);
    getline(stream, temp, ' ');
    vector_value.second=stod(temp);
    vector.push(vector_value);
}
template <class inputData>
InputGenericVector<inputData>::InputGenericVector(string const& path){ //For input file
    //The input file
    ifstream inputFile;
    try{
        inputFile.open(path);
        if (!inputFile.is_open())
            throw "Can't open file";
        while (!inputFile.eof()) {
            string sLine;
            //Read line by line
            getline(inputFile, sLine);
            if(sLine.length()==0)//Break if it's last line.
                break;
            size_t pos = 0;

            bool itemIDflag= true;
            while ((pos = sLine.find(' ')) != string::npos) {
                //Push every value delimited by space
                if(!itemIDflag) //If it is not the item's id
                    //this->push(stoi(sLine.substr(0, pos)));
                    constructorFunction(*this,sLine.substr(0, pos));
                else{
                    itemID=sLine.substr(0, pos);
                    itemIDflag= false;
                }
                sLine.erase(0, pos + 1);
            }
            //Save the pushed values
            this->save();
        }

        inputFile.close();
    }catch (const char* msg) {
        cerr << msg << endl;
    }
}
template <class inputData>
InputGenericVector<inputData>::InputGenericVector(string path,double& radius){ //For query file
    //The input file
    ifstream inputFile;
    try{
        inputFile.open(path);
        if (!inputFile.is_open())
            throw "Can't open file";
        while (!inputFile.eof()) {
            string sLine;
            //Read line by line
            getline(inputFile, sLine);
            if(sLine.length()==0)//Break if it's last line.
                break;
            size_t pos = 0;
            if(sLine.find("Radius:")!= string::npos){
                pos=sLine.find(' ');
                sLine.erase(0, pos + 1);
                pos=sLine.find(' ');
                radius=stod(sLine.substr(0, pos));
                continue;
            }


            bool itemIDflag= true;
            while ((pos = sLine.find(' ')) != string::npos) {
                //Push every value delimited by space
                if(!itemIDflag) //If it is not the item's id
                    //this->push(stoi(sLine.substr(0, pos)));
                    constructorFunction(*this,sLine.substr(0, pos));
                else{
                    itemID=sLine.substr(0, pos);
                    itemIDflag= false;
                }
                sLine.erase(0, pos + 1);
            }
            //Save the pushed values
            this->save();
        }

        inputFile.close();
    }catch (const char* msg) {
        cerr << msg << endl;
    }
}

template <>
InputGenericVector<pair<double,double>>::InputGenericVector(string const& path,unsigned int& maxCurveSize,unsigned int& minCurveSize,bool const& input){ //Specific for trajectories dataset. Wont work with other files due to line variable.
    //The input file
    ifstream inputFile;
    try{
        inputFile.open(path);
        if (!inputFile.is_open())
            throw "Can't open file";

        unsigned int line=0; //Number of line that is currently read
        unsigned int max_m=0; //Max size of all grids
        unsigned int min_m=4294967295;//Min size of all grids

        while (!inputFile.eof()) {
            string sLine;
            //Read line by line
            getline(inputFile, sLine);
            line++;
            if(line>7401 && input)
                break;
            if(line<7402 && !input)
                continue;
            if(sLine.length()==0)//Break if it's last line=empty line.
                break;
            size_t pos = 0;

            while ((pos = sLine.find('\t')) != string::npos) {
                itemID=sLine.substr(0, pos);

                sLine.erase(0, pos + 1);
                pos = sLine.find('\t');
                unsigned int curve_size=(unsigned int)stoi(sLine.substr(0, pos));
                if(curve_size>max_m)
                    max_m=curve_size;
                if(curve_size<min_m)
                    min_m=curve_size;
                sLine.erase(0, pos + 1);
                for(unsigned int i=0; i<curve_size; i++){
                    pos = sLine.find(',');
                    string temp=sLine.substr(0, pos);
                    temp.erase(0,1); //Remove first character of string which is (

                    sLine.erase(0, pos + 1);
                    pos = sLine.find(')');
                    temp+=sLine.substr(0, pos);
                    //cout<<"to temp pou stelnw "<<temp<<endl;
                    constructorFunction(*this, temp);
                    sLine.erase(0, pos + 2);
                }
            }

            //Save the pushed values
            this->save();
        }

        inputFile.close();
        if(max_m>maxCurveSize)
            maxCurveSize=max_m;
        if(min_m<minCurveSize)
            minCurveSize=min_m;
    }catch (const char* msg) {
        cerr << msg << endl;
    }
}
template class InputGenericVector<int>; //In order to not fail the compile as the compiler wants to see the data that the templated class will have.
template class InputGenericVector<pair<double,double>>;