#include <bits/stdc++.h>
using namespace std; 

double getMean(vector <double> vec) {
    double vecSize = vec.size() * 1.0;
    double sum = accumulate(vec.begin(), vec.end(), 0.0);

    return sum / vecSize;
}

double getMedian(vector <double> vec){
    sort(vec.begin(), vec.end());
    int medianPosition = vec.size() / 2;

    return vec.size() % 2 ? vec[medianPosition] : (vec[medianPosition] + vec[medianPosition - 1]) / 2.0 ;
}

double getMode(vector <double> vec){
    unordered_map<double, int> umap; 
    
    //frequency of a value in vec
    for (auto x: vec) umap[x]++; 
    
    double mode = 0.0;
    int maxCnt = 0;

    // Traversing an unordered map 
    for (auto x : umap) {
        if (x.second >= maxCnt) {
            maxCnt = x.second;
            mode = x.first;
        }
    }

    return mode;
}

double getFromND(vector <double> vec, double sigma){
    double mean = getMean(vec);
    double vecSize = vec.size() * 1.0;

    double squared_diff = 0.0;
    for (int i = 0; i < vec.size(); i++) {
        squared_diff += (vec[i] - mean) * (vec[i] - mean);
    }

    double standard_deviation = sqrt(squared_diff / vecSize);
    
    double low = mean - sigma * standard_deviation;
    double high = mean + sigma * standard_deviation;

    double pairCnt = 0.0;
    double sum = 0.0;

    for (int i = 0; i < vec.size(); i++) {
        if (vec[i] >= low && vec[i] <= high) {
            sum += vec[i];
            pairCnt += 1.0;
        }
    }

    return sum / pairCnt;
}
