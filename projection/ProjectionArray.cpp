
#include "ProjectionArray.h"

ProjectionArray::ProjectionArray(unsigned long n) {
    array.resize(n);
    for(int i=0; i<n; i++){
        array[i].resize(n);
    }

}

void ProjectionArray:: calculateArray(unsigned long size){

    for(int i=0; i<size; i++){
        for(int j=0; j<size; j++){
            if(abs(i-j)<4)
            calculatePoints((unsigned long)i,(unsigned long)j,array[i][j].RT);
            //displayPairs(array[i][j]);
        }
    }

}


void ProjectionArray:: calculatePoints(const unsigned long m, const unsigned long n,vector<vector<pairs>> &RT){

    set<pairs> points;
    for(int i=0; i<m; i++){     //calculate the diagonal pairs
        int j = (int)(n/m)*i;
        points.insert(pair<int,int>(i,j));

    }
    for(int j=0; j<n; j++){     //calculate the diagonal pairs
        int i = (int)(m/n)*j;
        points.insert(pair<int,int>(i,j));
    }

    for (auto it=points.begin(); it != points.end(); ++it) {
        if ((*it).first + 1 < m) {
            points.insert(pair<int, int>((*it).first + 1, (*it).second));         //add the upper points
        }
    }


    if(m==0&&n==0){
        return;
    }
    else if(m==0){
        vector<pairs> traversal;
        for(int j=0;j<n;j++){
            traversal.push_back(make_pair(0,j));
        }
        RT.push_back(traversal);
        return;
    }
    else if(n==0){
        vector<pairs> traversal;
        for(int i=0;i<m;i++){
            traversal.push_back(make_pair(i,0));
        }
        RT.push_back(traversal);
        return;
    }


    calculateAllPaths((int)m,(int)n,points,RT);



}
// The main function that calculates all paths from top left to bottom right
void ProjectionArray:: calculateAllPaths(int m, int n, const set<pairs>& s, vector<vector<pairs>> &RT ) {
    vector<pairs> path ;
    path.resize((unsigned long)m+n);

    calculateAllPathsUtil(0, 0, m, n, path, 0, s, RT);

}
void ProjectionArray:: calculateAllPathsUtil(int i, int j, int m, int n, vector<pairs> &path, int pi, const set<pairs>& s, vector<vector<pairs>> &RT) { // idea from: https://www.geeksforgeeks.org/print-all-possible-paths-from-top-left-to-bottom-right-of-a-mxn-matrix/
    // Reached the bottom of the matrix so we are left with
    // only option to move right
    if (i == m - 1) {
        for (int k = j; k < n; k++) {
            path[pi + k - j].first = i;
            path[pi + k - j].second = k;
        }
        int flag = 0;
        vector<pairs> traversal;
        for (int l = 0; l < pi + n - j; l++) {

            traversal.push_back(make_pair(path[l].first,path[l].second));

            if(s.count(make_pair(path[l].first,path[l].second))) {   //pair is in set, add it as Relevant Traversal

            }else{
                flag =1;
                break;
            }

        }
        if(flag == 0){
            RT.push_back(traversal);
        }
        traversal.clear();
        return;
    }

        // Reached the right corner of the matrix we are left with
        // only the downward movement.

    else if (j == n - 1) {
        vector<pairs> traversal;
        for (int k = i; k < m; k++) {
            path[pi + k - i].first = k;
            path[pi + k - i].second = j;
        }
        int flag = 0;
        for (int l = 0; l < pi + m - i; l++){

            traversal.push_back(make_pair(path[l].first,path[l].second));
            if(s.count(make_pair(path[l].first,path[l].second)) ){   //pair is in set, add it as Relevant Traversal

            }else{
                flag = 1;
                break;
            }

        }

        if(flag == 0){
            RT.push_back(traversal);
        }
        traversal.clear();
        return;
    }

    // Add the current cell to the path being generated
    path[pi].first=i;
    path[pi].second=j;

    // Print all the paths that are possible after moving down
    calculateAllPathsUtil(i+1, j, m, n, path, pi + 1, s, RT);

    // Print all the paths that are possible after moving right
    calculateAllPathsUtil(i, j+1, m, n, path, pi + 1, s, RT);

    // Print all the paths that are possible after moving diagonal
    calculateAllPathsUtil(i+1, j+1, m, n, path, pi + 1, s, RT);
}

double generateNumberG(double const& range_from,double const& range_to){      //https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution

    unsigned seed =
            chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);

    normal_distribution<double> distribution (range_from,range_to);
    return distribution(generator);

}

void calculateG( vector<double> &G, const int d, const double e){
    double K = (-d*log2(e))/(e*e);


    G.resize((unsigned long)K*d);

    for(int i=0; i<K*d; i++){
        G[i] = generateNumberG(0.0,1.0);
    }

}

void calculateVectors(vector<pair<double,double>> const &curve, vector<double> &G,cell &cell, string &curveId,int ipos, int jpos, int flag){

    int i=0;
    for(auto& rel_tr:cell.RT) {
        createVector(curve, G, rel_tr, curveId,ipos,jpos, cell, i,flag);
        i++;
    }

}

void createVector(vector<pair<double,double>> const &curve, vector<double> &G, vector<pairs> &RT, string &curveId, int ipos, int jpos, cell & cell, int traversalpos, int flag){
    vector<double> x;
    x.resize(RT.size());

    for(int i=0; i<(int)RT.size(); i++){
        if(flag==0) {
            int pos = RT[i].first;          //find the position of the first int of the traversal
            x[i] = curve[pos].first;
        }else{
            int pos = RT[i].second;          //find the position of the second int of the traversal
            x[i] = curve[pos].second;
        }
    }
    vector<vector<double>> t;
    t.resize(G.size());
    for(int i=0; i<(int)G.size();i++){
        t[i].resize(x.size());
    }

    for(int i=0; i<(int)G.size(); i++) {            //multiply G with every x point
        for (int j = 0; j < (int) x.size(); j++) {
            t[i][j] = G[i] * x[j];
        }
    }

    vector<double> concat;
    concat.resize(G.size());
    for(int i=0; i<(int)t.size();i++){
        double sum=0;
        for(int j=0; j<(int)t[i].size(); j++){
            sum += t[i][j];
        }
        concat[i] = sum;
    }



    string key =to_string(ipos)+"-"+to_string(jpos)+"-"+to_string(traversalpos);
    if(flag == 0) {
        cell.imap[key].push_back(make_pair(ipos, concat));   //for points
    }else{
        cell.jmap[key].push_back(make_pair(curveId, concat)); //for queries
    }


}