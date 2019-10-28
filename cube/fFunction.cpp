#include "fFunction.h"

inline int generateNumber(int const& range_from,int const& range_to){      //https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution

    random_device                  rand_dev;
    mt19937                        generator(rand_dev());
    uniform_int_distribution<>  distr(range_from, range_to);

    return distr(generator);

}
inline bool check_key(unordered_map<int, string> const& map, int key)
{
    // Key is not present
    if (map.find(key) == map.end())
        return false;

    return true;
}
template <class inputData>
inline string fFunctions<inputData>::binaryStringCalculator(vector<inputData> const& p,int d,int const& k,double const& w,vector<HashFunctions<inputData>> hashFunctions){
    string binaryString;
    for(int i=0; i<d; i++){
        unsigned int g=hashFunctions[i].gCalculator(p,(unsigned int)k,w);
        string f;
        if(!check_key(alreadySeenValues,g)){ //Check if we have already encountered this g
            f=to_string(generateNumber(0,1)); //If not, label it 0 or 1
            alreadySeenValues.insert(make_pair(g, f)); //Save it for future use
        }else{
            auto binaryValue = alreadySeenValues.find(g); //if we have already encountered this g,then find the label 0 or 1 that it has
            f=binaryValue->second;
        }
        binaryString+=f;
    }
    return binaryString;
}
template class fFunctions<int>; //In order to not fail the compile as the compiler wants to see the data that the templated class will have.