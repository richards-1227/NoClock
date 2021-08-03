//
//  main.cpp
//  noClock
//
//  Created by Andrew Richards on 11/25/20.
//  Copyright © 2020 Andrew Richards. All rights reserved.
//

#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <map>
#include <random>
#include <algorithm>
#include <cstdlib> // for exit()


struct collectPattern{
    int patternID;
    int gIndex;
    int integralID;
    int hIndex[6];
    int count;
};

std::vector<collectPattern> getProbSym(double t1, double t2, double t3, double gammaA, double gammaB, double gammaC, double gammaD, double gammaAB, double gammaCD, double theta, std::vector<collectPattern> s);
std::vector<collectPattern> getProbAsymm(double t1, double t2, double t3, double gammaA, double gammaB, double gammaC, double gammaD, double gammaCD, double gammaBCD, double theta, std::vector<collectPattern> s);
std::vector<collectPattern> consolidatePattern(std::vector<collectPattern> t);

std::vector<collectPattern> consolidatePattern(std::vector<collectPattern> t){
    std::vector<collectPattern> s{};
    std::vector<int> consolidatedInd {};
    for(int i=0; i<t.size(); ++i){
        consolidatedInd.push_back(0);
    }
    
    for(int i=0; i<t.size(); ++i){
        if(consolidatedInd[i]==0){

            for(int j=(i+1); j<t.size(); ++j){
                if(consolidatedInd[j]==0){
                    if(t[i].patternID==t[j].patternID && t[i].integralID==t[j].integralID
                       && t[i].gIndex==t[j].gIndex && t[i].hIndex[0]==t[j].hIndex[0] && t[i].hIndex[1]==t[j].hIndex[1] && t[i].hIndex[2]==t[j].hIndex[2]
                       && t[i].hIndex[3]==t[j].hIndex[3] && t[i].hIndex[4]==t[j].hIndex[4] && t[i].hIndex[5]==t[j].hIndex[5]){
                        consolidatedInd[j]=1;
                        t[i].count += t[j].count;
                    }
                }
            }
                        s.push_back(t[i]);
        }
    }
    
    return s;
}

std::vector<collectPattern> getProbSym(double t1, double t2, double t3, double gammaA, double gammaB, double gammaC, double gammaD, double gammaAB, double gammaCD, double theta, std::vector<collectPattern> s){
    std::vector<double> y{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double a{4.0/3};
     std::vector<std::vector<int>> coefMat{
        {1,3,3,3,3,3,3,9,6,6,6,6,12},
        {1,3,-1,3,-1,3,-1,-3,6,-2,-2,-2,-4},
        {1,3,-1,-1,3,-1,3,-3,-2,6,-2,-2,-4},
        {1,-1,3,3,3,-1,-1,-3,-2,-2,6,-2,-4},
        {1,-1,3,-1,-1,3,3,-3,-2,-2,-2,6,-4},
        {1,3,3,-1,-1,-1,-1,9,-2,-2,-2,-2,-4},
        {1,-1,-1,3,-1,-1,3,1,-2,-2,-2,-2,4},
        {1,-1,-1,-1,3,3,-1,1,-2,-2,-2,-2,4},
        {1,3,-1,-1,-1,-1,-1,-3,-2,-2,2,2,4},
        {1,-1,-1,3,-1,-1,-1,1,-2,2,-2,2,0},
        {1,-1,-1,-1,3,-1,-1,1,2,-2,-2,2,0},
        {1,-1,-1,-1,-1,3,-1,1,-2,2,2,-2,0},
        {1,-1,-1,-1,-1,-1,3,1,2,-2,2,-2,0},
        {1,-1,3,-1,-1,-1,-1,-3,2,2,-2,-2,4},
        {1,-1,-1,-1,-1,-1,-1,1,2,2,2,2,-4}
    };
    
    std::vector<std::vector<int>> permMat{
        {1,1,1,1,1,1,1,1,1,1,1,1},
        {2,2,4,2,2,4,2,2,4,5,5,5},
        {3,4,2,3,4,2,5,5,5,2,2,4},
        {4,3,3,5,5,5,3,4,2,3,4,2},
        {5,5,5,4,3,3,4,3,3,4,3,3},
        {6,7,7,6,7,7,8,8,6,8,8,6},
        {7,6,8,8,8,6,6,7,7,7,6,8},
        {8,8,6,7,6,8,7,6,8,6,7,7},
        {9,10,10,9,10,10,12,12,14,12,12,14},
        {10,9,11,12,12,14,9,10,10,13,14,12},
        {11,11,9,13,14,12,13,14,12,9,10,10},
        {12,12,14,10,9,11,10,9,11,14,13,13},
        {13,14,12,11,11,9,14,13,13,10,9,11},
        {14,13,13,14,13,13,11,11,9,11,11,9},
        {15,15,15,15,15,15,15,15,15,15,15,15}
    };
    
    std::vector<std::vector<double>> g3Mat { {0,0,0,2,2,2,2,0, 2,2,2,2,2},
    {0,2,0,2,2,0,0,2, 2,2,2,0,2}
    };
    
     std::vector<int> g3Index {0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    std::vector<std::vector<double>> g2Mat {
        {0, 2, 0,0,0,0,0, 2, 1, 1, 0,0, 1 },
        {0, 2, 0,0,0,0,0, 2, 1, 1, 0,0, 1},
        {0,0,0,0,0, 2, 2, 0, 1, 1, 0, 2, 1}
    };
    
     std::vector<int> g2Index {0,1,0,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
    
    
     std::vector<std::vector<double>> g1Mat{
        {0,0, 2, 0,0,0,0, 2, 0,0, 1, 1,1},
        {0,0, 2, 0,0,0,0, 2, 0,0, 1, 1, 1},
        {0,0, 2, 0,0,0,0, 2, 0,0, 1, 1,1}
    };
    std::vector<int> g1Index {0,0,1,1,1,1,1,1,1,0,1,1,1,0,1,1,1,1,1,2,1,1,1,2,1};

    
    std::vector<std::vector<std::vector<int>>> g0Mat{
        { {0,0,0,0,0,0},
            {0,0,1,1,0,0}, {1,1,0,0,0,0}, {1,0,1,0,1,1}, {0,1,1,0,1,1}, {1,0,0,1,1,1}, {0,1,0,1,1,1},
            {1,1,1,1,0,0}, {1,0,1,1,1,1}, {0,1,1,1,1,1}, {1,1,1,0,1,1}, {1,1,0,1,1,1}, {1,1,1,1,1,1}
        },
        {  {0,0,0,0,0,0},
            {0,0,1,1,2,0}, {1,1,0,0,0,0}, {1,0,1,0,1,1}, {0,1,1,0,1,1}, {1,0,0,1,1,1}, {0,1,0,1,1,1},
            {1,1,1,1,2,0}, {1,0,1,1,2,1}, {0,1,1,1,2,1}, {1,1,1,0,1,1}, {1,1,0,1,1,1}, {1,1,1,1,2,1}
        },
        {  {0,0,0,0,0,0},
            {0,0,1,1,0,0}, {1,1,0,0,0,2}, {1,0,1,0,1,1}, {0,1,1,0,1,1}, {1,0,0,1,1,1}, {0,1,0,1,1,1},
            {1,1,1,1,0,2}, {1,0,1,1,1,1}, {0,1,1,1,1,1}, {1,1,1,0,1,2}, {1,1,0,1,1,2}, {1,1,1,1,1,2}
        },
        {  {0,0,0,0,0,0},
            {0,0,1,1,2,0}, {1,1,0,0,0,2}, {1,0,1,0,1,1}, {0,1,1,0,1,1}, {1,0,0,1,1,1}, {0,1,0,1,1,1},
            {1,1,1,1,2,2}, {1,0,1,1,2,1}, {0,1,1,1,2,1}, {1,1,1,0,1,2}, {1,1,0,1,1,2}, {1,1,1,1,2,2}
        },
        {  {0,0,0,0,0,0},
            {0,0,1,1,2,0}, {1,1,0,0,0,2}, {1,0,1,0,1,1}, {0,1,1,0,1,1}, {1,0,0,1,1,1}, {0,1,0,1,1,1},
            {1,1,1,1,2,2}, {1,0,1,1,2,1}, {0,1,1,1,2,1}, {1,1,1,0,1,2}, {1,1,0,1,1,2}, {1,1,1,1,2,2}
        },
        {  {0,0,0,0,0,0},
            {1,0,1,0,1,1}, {0,1,0,1,1,1}, {0,0,1,1,2,0}, {0,1,1,0,1,1}, {1,0,0,1,1,1}, {1,1,0,0,0,2},
            {1,1,1,1,2,2}, {1,0,1,1,2,1}, {1,1,1,0,1,2}, {0,1,1,1,2,1}, {1,1,0,1,1,2}, {1,1,1,1,2,2}
        },
        {  {0,0,0,0,0,0},
            {1,0,1,0,1,1}, {0,1,0,1,1,1}, {0,0,1,1,2,0}, {0,1,1,0,1,1}, {1,0,0,1,1,1}, {1,1,0,0,0,2},
            {1,1,1,1,2,2}, {1,0,1,1,2,1}, {1,1,1,0,1,2}, {0,1,1,1,2,1}, {1,1,0,1,1,2}, {1,1,1,1,2,2}
        },
        {  {0,0,0,0,0,0},
            {0,1,1,0,1,1}, {1,0,0,1,1,1}, {0,0,1,1,2,0}, {1,0,1,0,1,1}, {0,1,0,1,1,1}, {1,1,0,0,0,2},
            {1,1,1,1,2,2}, {0,1,1,1,2,1}, {1,1,1,0,1,2}, {1,0,1,1,2,1}, {1,1,0,1,1,2}, {1,1,1,1,2,2}
        },
        {  {0,0,0,0,0,0},
            {0,1,1,0,1,1}, {1,0,0,1,1,1}, {0,0,1,1,2,0}, {1,0,1,0,1,1}, {0,1,0,1,1,1}, {1,1,0,0,0,2},
            {1,1,1,1,2,2}, {0,1,1,1,2,1}, {1,1,1,0,1,2}, {1,0,1,1,2,1}, {1,1,0,1,1,2}, {1,1,1,1,2,2}
        },
        {  {0,0,0,0,0,0},
            {0,0,1,1,2,0}, {1,1,0,0,0,0}, {1,0,1,0,1,1}, {0,1,1,0,1,1}, {1,0,0,1,1,1}, {0,1,0,1,1,1},
            {1,1,1,1,2,0}, {1,0,1,1,2,1}, {0,1,1,1,2,1}, {1,1,1,0,1,1}, {1,1,0,1,1,1}, {1,1,1,1,2,1}
        },
        {  {0,0,0,0,0,0},
            {0,0,1,1,2,0}, {1,1,0,0,0,2}, {1,0,1,0,1,1}, {0,1,1,0,1,1}, {1,0,0,1,1,1}, {0,1,0,1,1,1},
            {1,1,1,1,2,2}, {1,0,1,1,2,1}, {0,1,1,1,2,1}, {1,1,1,0,1,2}, {1,1,0,1,1,2}, {1,1,1,1,2,2}
        },
        {  {0,0,0,0,0,0},
            {1,0,1,0,1,1}, {0,1,0,1,1,1}, {0,0,1,1,2,0}, {0,1,1,0,1,1}, {1,0,0,1,1,1}, {1,1,0,0,0,2},
            {1,1,1,1,2,2}, {1,0,1,1,2,1}, {1,1,1,0,1,2}, {0,1,1,1,2,1}, {1,1,0,1,1,2}, {1,1,1,1,2,2}
        },
        {  {0,0,0,0,0,0},
            {0,1,1,0,1,1}, {1,0,0,1,1,1}, {0,0,1,1,2,0}, {1,0,1,0,1,1}, {0,1,0,1,1,1}, {1,1,0,0,0,2},
            {1,1,1,1,2,2}, {0,1,1,1,2,1}, {1,1,1,0,1,2}, {1,0,1,1,2,1}, {1,1,0,1,1,2}, {1,1,1,1,2,2}
        },
        {  {0,0,0,0,0,0},
            {0,0,1,1,2,0}, {1,1,0,0,0,0}, {1,0,0,1,1,1}, {0,1,0,1,1,1}, {1,0,1,0,1,1}, {0,1,1,0,1,1},
            {1,1,1,1,2,0}, {1,0,1,1,2,1}, {0,1,1,1,2,1}, {1,1,0,1,1,1}, {1,1,1,0,1,1}, {1,1,1,1,2,1}
        },
        {  {0,0,0,0,0,0},
            {0,0,1,1,2,0}, {1,1,0,0,0,2}, {1,0,0,1,1,1}, {0,1,0,1,1,1}, {1,0,1,0,1,1}, {0,1,1,0,1,1},
            {1,1,1,1,2,2}, {1,0,1,1,2,1}, {0,1,1,1,2,1}, {1,1,0,1,1,2}, {1,1,1,0,1,2}, {1,1,1,1,2,2}
        },
        {  {0,0,0,0,0,0},
            {1,0,0,1,1,1}, {0,1,1,0,1,1}, {0,0,1,1,2,0}, {0,1,0,1,1,1}, {1,0,1,0,1,1}, {1,1,0,0,0,2},
            {1,1,1,1,2,2}, {1,0,1,1,2,1}, {1,1,0,1,1,2}, {0,1,1,1,2,1}, {1,1,1,0,1,2}, {1,1,1,1,2,2}
        },
        {  {0,0,0,0,0,0},
            {0,1,0,1,1,1}, {1,0,1,0,1,1}, {0,0,1,1,2,0}, {1,0,0,1,1,1}, {0,1,1,0,1,1}, {1,1,0,0,0,2},
            {1,1,1,1,2,2}, {0,1,1,1,2,1}, {1,1,0,1,1,2}, {1,0,1,1,2,1}, {1,1,1,0,1,2}, {1,1,1,1,2,2}
        },
        {  {0,0,0,0,0,0},
            {1,0,1,0,1,1}, {0,1,0,1,1,1}, {1,0,0,1,1,1}, {1,1,0,0,0,2}, {0,0,1,1,2,0}, {0,1,1,0,1,1},
            {1,1,1,1,2,2}, {1,0,1,1,2,1}, {1,1,1,0,1,2}, {1,1,0,1,1,2}, {0,1,1,1,2,1}, {1,1,1,1,2,2}
        },
        {  {0,0,0,0,0,0},
            {1,0,0,1,1,1}, {0,1,1,0,1,1}, {1,0,1,0,1,1}, {1,1,0,0,0,2}, {0,0,1,1,2,0}, {0,1,0,1,1,1},
            {1,1,1,1,2,2}, {1,0,1,1,2,1}, {1,1,0,1,1,2}, {1,1,1,0,1,2}, {0,1,1,1,2,1}, {1,1,1,1,2,2}
        },
        {  {0,0,0,0,0,0},
            {1,1,0,0,0,2}, {0,0,1,1,0,0}, {1,0,1,0,1,1}, {1,0,0,1,1,1}, {0,1,1,0,1,1}, {0,1,0,1,1,1},
            {1,1,1,1,0,2}, {1,1,1,0,1,2}, {1,1,0,1,1,2}, {1,0,1,1,1,1}, {0,1,1,1,1,1}, {1,1,1,1,1,2}
        },
        {  {0,0,0,0,0,0},
            {1,1,0,0,0,2}, {0,0,1,1,2,0}, {1,0,1,0,1,1}, {1,0,0,1,1,1}, {0,1,1,0,1,1}, {0,1,0,1,1,1},
            {1,1,1,1,2,2}, {1,1,1,0,1,2}, {1,1,0,1,1,2}, {1,0,1,1,2,1}, {0,1,1,1,2,1}, {1,1,1,1,2,2}
        },
        {  {0,0,0,0,0,0},
            {0,1,1,0,1,1}, {1,0,0,1,1,1}, {0,1,0,1,1,1}, {1,1,0,0,0,2}, {0,0,1,1,2,0}, {1,0,1,0,1,1},
            {1,1,1,1,2,2}, {0,1,1,1,2,1}, {1,1,1,0,1,2}, {1,1,0,1,1,2}, {1,0,1,1,2,1}, {1,1,1,1,2,2}
        },
        {  {0,0,0,0,0,0},
            {0,1,0,1,1,1}, {1,0,1,0,1,1}, {0,1,1,0,1,1}, {1,1,0,0,0,2}, {0,0,1,1,2,0}, {1,0,0,1,1,1},
            {1,1,1,1,2,2}, {0,1,1,1,2,1}, {1,1,0,1,1,2}, {1,1,1,0,1,2}, {1,0,1,1,2,1}, {1,1,1,1,2,2}
        },
        {  {0,0,0,0,0,0},
            {1,1,0,0,0,2}, {0,0,1,1,0,0}, {0,1,1,0,1,1}, {0,1,0,1,1,1}, {1,0,1,0,1,1}, {1,0,0,1,1,1},
            {1,1,1,1,0,2}, {1,1,1,0,1,2}, {1,1,0,1,1,2}, {0,1,1,1,1,1}, {1,0,1,1,1,1}, {1,1,1,1,1,2}
        },
        {  {0,0,0,0,0,0},
            {1,1,0,0,0,2}, {0,0,1,1,2,0}, {0,1,1,0,1,1}, {0,1,0,1,1,1}, {1,0,1,0,1,1}, {1,0,0,1,1,1},
            {1,1,1,1,2,2}, {1,1,1,0,1,2}, {1,1,0,1,1,2}, {0,1,1,1,2,1}, {1,0,1,1,2,1}, {1,1,1,1,2,2}
        }
    };
    
    std::vector<int> integralIndex{1,2,3,4,5,4,5,4,5,2,4,4,4,2,4,4,4,4,4,6,4,4,4,6,4};
    std::vector<int> permIndex{0,0,0,0,0,1,1,2,2,0,0,1,2,3,3,4,5,6,7,8,8,9,10,11,11};
    
    for(int m=0; m<15; ++m){
        for(int i=0; i<25; ++i){
            for(int j=0; j<13; ++j){
                int coefMult {};
                int index1 {permMat[m][permIndex[i]]-1};
                coefMult = coefMat[index1][j];
                collectPattern r{};
                r.patternID=m;
                r.integralID=integralIndex[i];
                int g1{};
                g1=g1Mat[g1Index[i]][j];
                int g2{};
                g2=g2Mat[g2Index[i]][j];
                int g3{};
                g3=g3Mat[g3Index[i]][j];
                if(g1==2){
                    r.gIndex=1;
                }
                if(g2==2){
                    r.gIndex=2;
                }
                if(g3==2){
                    r.gIndex=3;
                }
                if(g1==1 && g2==2){
                    r.gIndex=4;
                }
                if(g1==1 && g3==2){
                    r.gIndex=5;
                }
                if(g2==1 && g3==2){
                    r.gIndex=6;
                }
                if(g1==1 && g2==1 && g3==2){
                    r.gIndex=7;
                }
                if(g1==2 && g2==2 ){
                    r.gIndex=8;
                }
                if(g1==2 && g3==2){
                    r.gIndex=9;
                }
                for(int k=0; k<6; ++k){
                    r.hIndex[k]=g0Mat[i][j][k];
                }
                r.count=coefMult;
                s.push_back(r);
        }
    }
}
    
    return s;
}

std::vector<collectPattern> getProbAsymm(double t1, double t2, double t3, double gammaA, double gammaB, double gammaC, double gammaD, double gammaCD, double gammaBCD, double theta, std::vector<collectPattern> s){
    std::vector<double> z{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double a{4.0/3};
    std::vector<std::vector<int>> coefMat{
        {1,3,3,3,3,3,3,9,6,6,6,6,12},
        {1,3,-1,3,-1,3,-1,-3,6,-2,-2,-2,-4},
        {1,3,-1,-1,3,-1,3,-3,-2,6,-2,-2,-4},
        {1,-1,3,3,3,-1,-1,-3,-2,-2,6,-2,-4},
        {1,-1,3,-1,-1,3,3,-3,-2,-2,-2,6,-4},
        {1,3,3,-1,-1,-1,-1,9,-2,-2,-2,-2,-4},
        {1,-1,-1,3,-1,-1,3,1,-2,-2,-2,-2,4},
        {1,-1,-1,-1,3,3,-1,1,-2,-2,-2,-2,4},
        {1,3,-1,-1,-1,-1,-1,-3,-2,-2,2,2,4},
        {1,-1,-1,3,-1,-1,-1,1,-2,2,-2,2,0},
        {1,-1,-1,-1,3,-1,-1,1,2,-2,-2,2,0},
        {1,-1,-1,-1,-1,3,-1,1,-2,2,2,-2,0},
        {1,-1,-1,-1,-1,-1,3,1,2,-2,2,-2,0},
        {1,-1,3,-1,-1,-1,-1,-3,2,2,-2,-2,4},
        {1,-1,-1,-1,-1,-1,-1,1,2,2,2,2,-4}
    };
    
    std::vector<std::vector<int>> permMat{
        {1,1,1,1,1,1,1,1,1,1,1,1},
        {2,2,4,2,2,4,2,2,4,5,5,5},
        {3,4,2,3,4,2,5,5,5,2,2,4},
        {4,3,3,5,5,5,3,4,2,3,4,2},
        {5,5,5,4,3,3,4,3,3,4,3,3},
        {6,7,7,6,7,7,8,8,6,8,8,6},
        {7,6,8,8,8,6,6,7,7,7,6,8},
        {8,8,6,7,6,8,7,6,8,6,7,7},
        {9,10,10,9,10,10,12,12,14,12,12,14},
        {10,9,11,12,12,14,9,10,10,13,14,12},
        {11,11,9,13,14,12,13,14,12,9,10,10},
        {12,12,14,10,9,11,10,9,11,14,13,13},
        {13,14,12,11,11,9,14,13,13,10,9,11},
        {14,13,13,14,13,13,11,11,9,11,11,9},
        {15,15,15,15,15,15,15,15,15,15,15,15}
    };
    
   std::vector<std::vector<double>> g3Mat { {0,0,0,2.0*a,2.0*a,2.0*a,2.0*a,0, 2.0*a,2.0*a,2.0*a,2.0*a,2.0*a},
    {0,2.0*a,0,2.0*a,2.0*a,0,0,2.0*a, 2.0*a,2.0*a,2.0*a,0,2.0*a}
    };
    
    std::vector<int> g3Index {0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    
     std::vector<std::vector<double>> g2Mat {
        {0,0,0,0,0, 2.0*a*gammaBCD, 2.0*a*gammaBCD, 0, a*gammaBCD, a*gammaBCD, 0, 2.0*a*gammaBCD, a*gammaBCD},
        {0, 2.0*a, 0,0,0,0,0, 2.0*a, a, a, 0,0, a},
        {0,0,0,0,0, 2.0*a, 2.0*a, 0,a, a, 0, 2.0*a, a}
    };
    
    std::vector<int> g2Index{1,1,1,1,1,1,1,1,1,1,0,2,0,2,2,0,2,2,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
    
     std::vector<std::vector<double>> g1Mat{
        {0,0, 2.0*a*gammaCD, 0,0,0,0, 2.0*a*gammaCD, 0,0, a*gammaCD, a*gammaCD,a*gammaCD},
        {0,0, 2.0*a, 0,0,0,0, 2.0*a, 0,0, a, a, a},
        {0,0, 2.0*a*gammaBCD, 0,0,0,0, 2.0*a*gammaBCD, 0,0, a*gammaBCD, a*gammaBCD,a*gammaBCD}
    };
    
    std::vector<int> g1Index{0,2,1,1,2,1,1,2,1,1,0,0,2,2,1,2,2,1,2,2,1,0,2,1,1,1,2,1,1,1,2,1,1,1};
    
    std::vector<double> g0Vec{exp(-a*t1*(gammaC-gammaCD)), exp(-a*t1*(gammaD-gammaCD)), exp(-a*t2*(gammaB-gammaBCD)), exp(-a*t2*(gammaCD-gammaBCD)), exp(-a*t3*(gammaA-1.0)), exp(-a*t3*(gammaBCD-1.0)) };
    
       std::vector<std::vector<std::vector<int>>> g0Mat{
        { {0,0,0,0,0,0},
           {0,0,1,0,1,1},{1,1,0,0,0,0},{1,0,0,1,1,1},{0,1,0,1,1,1},{1,0,1,1,0,2},{0,1,1,1,0,2},
           {1,1,1,0,1,1},{1,0,1,1,1,2},{0,1,1,1,1,2},{1,1,0,1,1,1},{1,1,1,1,0,2},{1,1,1,1,1,2},
        },//1
        { {0,0,0,0,0,0},
           {0,0,1,0,1,1},{1,1,0,2,0,0},{1,0,0,1,1,1},{0,1,0,1,1,1},{1,0,1,1,0,2},{0,1,1,1,0,2},
           {1,1,1,2,1,1},{1,0,1,1,1,2},{0,1,1,1,1,2},{1,1,0,2,1,1},{1,1,1,2,0,2},{1,1,1,2,1,2},
        },//2
        { {0,0,0,0,0,0},
           {0,0,1,0,1,1},{1,1,0,2,0,2},{1,0,0,1,1,1},{0,1,0,1,1,1},{1,0,1,1,0,2},{0,1,1,1,0,2},
           {1,1,1,2,1,3},{1,0,1,1,1,2},{0,1,1,1,1,2},{1,1,0,2,1,2},{1,1,1,2,0,3},{1,1,1,2,1,3},
        },//3
        { {0,0,0,0,0,0},
           {0,0,1,0,1,1},{1,1,0,2,0,2},{1,0,0,1,1,1},{0,1,0,1,1,1},{1,0,1,1,0,2},{0,1,1,1,0,2},
           {1,1,1,2,1,3},{1,0,1,1,1,2},{0,1,1,1,1,2},{1,1,0,2,1,2},{1,1,1,2,0,3},{1,1,1,2,1,3},
        },//4
        { {0,0,0,0,0,0},
           {1,0,0,1,1,1},{0,1,1,1,0,0},{0,0,1,0,1,1},{0,1,0,1,1,1},{1,0,1,1,0,2},{1,1,0,2,0,2},
           {1,1,1,2,1,1},{1,0,1,1,1,2},{1,1,0,2,1,2},{0,1,1,1,1,1},{1,1,1,2,0,2},{1,1,1,2,1,2},
        },//5
        { {0,0,0,0,0,0},
           {1,0,0,1,1,1},{0,1,1,1,0,2},{0,0,1,0,1,1},{0,1,0,1,1,1},{1,0,1,1,0,2},{1,1,0,2,0,2},
           {1,1,1,2,1,3},{1,0,1,1,1,2},{1,1,0,2,1,2},{0,1,1,1,1,2},{1,1,1,2,0,3},{1,1,1,2,1,3},
        },//6
        { {0,0,0,0,0,0},
           {1,0,0,1,1,1},{0,1,1,1,0,2},{0,0,1,0,1,1},{0,1,0,1,1,1},{1,0,1,1,0,2},{1,1,0,2,0,2},
           {1,1,1,2,1,3},{1,0,1,1,1,2},{1,1,0,2,1,2},{0,1,1,1,1,2},{1,1,1,2,0,3},{1,1,1,2,1,3},
        },//7
        { {0,0,0,0,0,0},
           {0,1,0,1,1,1},{1,0,1,1,0,0},{0,0,1,0,1,1},{1,0,0,1,1,1},{0,1,1,1,0,2},{1,1,0,2,0,2},
           {1,1,1,2,1,1},{0,1,1,1,1,2},{1,1,0,2,1,2},{1,0,1,1,1,1},{1,1,1,2,0,2},{1,1,1,2,1,2},
        },//8
        { {0,0,0,0,0,0},
           {0,1,0,1,1,1},{1,0,1,1,0,2},{0,0,1,0,1,1},{1,0,0,1,1,1},{0,1,1,1,0,2},{1,1,0,2,0,2},
           {1,1,1,2,1,3},{0,1,1,1,1,2},{1,1,0,2,1,2},{1,0,1,1,1,2},{1,1,1,2,0,3},{1,1,1,2,1,3},
        },//9
        { {0,0,0,0,0,0},
           {0,1,0,1,1,1},{1,0,1,1,0,2},{0,0,1,0,1,1},{1,0,0,1,1,1},{0,1,1,1,0,2},{1,1,0,2,0,2},
           {1,1,1,2,1,3},{0,1,1,1,1,2},{1,1,0,2,1,2},{1,0,1,1,1,2},{1,1,1,2,0,3},{1,1,1,2,1,3},
        },//10
        { {0,0,0,0,0,0},
           {0,0,1,0,1,1},{1,1,0,0,0,0},{1,0,0,1,1,1},{0,1,0,1,1,1},{1,0,1,1,0,0},{0,1,1,1,0,0},
           {1,1,1,0,1,1},{1,0,1,1,1,1},{0,1,1,1,1,1},{1,1,0,1,1,1},{1,1,1,1,0,0},{1,1,1,1,1,1},
        },//11
        { {0,0,0,0,0,0},
           {0,0,1,0,1,1},{1,1,0,0,0,0},{1,0,0,1,1,1},{0,1,0,1,1,1},{1,0,1,1,0,2},{0,1,1,1,0,2},
           {1,1,1,0,1,1},{1,0,1,1,1,2},{0,1,1,1,1,2},{1,1,0,1,1,1},{1,1,1,1,0,2},{1,1,1,1,1,2},
        },//12
        { {0,0,0,0,0,0},
           {0,0,1,0,1,1},{1,1,0,2,0,0},{1,0,0,1,1,1},{0,1,0,1,1,1},{1,0,1,1,0,0},{0,1,1,1,0,0},
           {1,1,1,2,1,1},{1,0,1,1,1,1},{0,1,1,1,1,1},{1,1,0,2,1,1},{1,1,1,2,0,0},{1,1,1,2,1,1},
        },//13
        { {0,0,0,0,0,0},
           {0,0,1,0,1,1},{1,1,0,2,0,0},{1,0,0,1,1,1},{0,1,0,1,1,1},{1,0,1,1,0,2},{0,1,1,1,0,2},
           {1,1,1,2,1,1},{1,0,1,1,1,2},{0,1,1,1,1,2},{1,1,0,2,1,1},{1,1,1,2,0,2},{1,1,1,2,1,2},
        },//14
        { {0,0,0,0,0,0},
           {0,0,1,0,1,1},{1,1,0,2,0,2},{1,0,0,1,1,1},{0,1,0,1,1,1},{1,0,1,1,0,2},{0,1,1,1,0,2},
           {1,1,1,2,1,3},{1,0,1,1,1,2},{0,1,1,1,1,2},{1,1,0,2,1,2},{1,1,1,2,0,3},{1,1,1,2,1,3},
        },//15
        { {0,0,0,0,0,0},
           {1,0,0,1,1,1},{0,1,1,1,0,0},{0,0,1,0,1,1},{0,1,0,1,1,1},{1,0,1,1,0,0},{1,1,0,2,0,0},
           {1,1,1,2,1,1},{1,0,1,1,1,1},{1,1,0,2,1,1},{0,1,1,1,1,1},{1,1,1,2,0,0},{1,1,1,2,1,1},
        },//16
        { {0,0,0,0,0,0},
           {1,0,0,1,1,1},{0,1,1,1,0,0},{0,0,1,0,1,1},{0,1,0,1,1,1},{1,0,1,1,0,2},{1,1,0,2,0,2},
           {1,1,1,2,1,1},{1,0,1,1,1,2},{1,1,0,2,1,2},{0,1,1,1,1,1},{1,1,1,2,0,2},{1,1,1,2,1,2},
        },//17
        { {0,0,0,0,0,0},
           {1,0,0,1,1,1},{0,1,1,1,0,2},{0,0,1,0,1,1},{0,1,0,1,1,1},{1,0,1,1,0,2},{1,1,0,2,0,2},
           {1,1,1,2,1,3},{1,0,1,1,1,2},{1,1,0,2,1,2},{0,1,1,1,1,2},{1,1,1,2,0,3},{1,1,1,2,1,3},
        },//18
        { {0,0,0,0,0,0},
           {0,1,0,1,1,1},{1,0,1,1,0,0},{0,0,1,0,1,1},{1,0,0,1,1,1},{0,1,1,1,0,0},{1,1,0,2,0,0},
           {1,1,1,2,1,1},{0,1,1,1,1,1},{1,1,0,2,1,1},{1,0,1,1,1,1},{1,1,1,2,0,0},{1,1,1,2,1,1},
        },//19
        { {0,0,0,0,0,0},
           {0,1,0,1,1,1},{1,0,1,1,0,0},{0,0,1,0,1,1},{1,0,0,1,1,1},{0,1,1,1,0,2},{1,1,0,2,0,2},
           {1,1,1,2,1,1},{0,1,1,1,1,2},{1,1,0,2,1,2},{1,0,1,1,1,1},{1,1,1,2,0,2},{1,1,1,2,1,2},
        },//20
        { {0,0,0,0,0,0},
           {0,1,0,1,1,1},{1,0,1,1,0,2},{0,0,1,0,1,1},{1,0,0,1,1,1},{0,1,1,1,0,2},{1,1,0,2,0,2},
           {1,1,1,2,1,3},{0,1,1,1,1,2},{1,1,0,2,1,2},{1,0,1,1,1,2},{1,1,1,2,0,3},{1,1,1,2,1,3},
        },//21
        { {0,0,0,0,0,0},
           {0,0,1,0,1,1},{1,1,0,0,0,0},{1,0,1,1,0,2},{0,1,1,1,0,2},{1,0,0,1,1,1},{0,1,0,1,1,1},
           {1,1,1,0,1,1},{1,0,1,1,1,2},{0,1,1,1,1,2},{1,1,1,1,0,2},{1,1,0,1,1,1},{1,1,1,1,1,2},
        },//22
        { {0,0,0,0,0,0},
           {0,0,1,0,1,1},{1,1,0,2,0,0},{1,0,1,1,0,2},{0,1,1,1,0,2},{1,0,0,1,1,1},{0,1,0,1,1,1},
           {1,1,1,2,1,1},{1,0,1,1,1,2},{0,1,1,1,1,2},{1,1,1,2,0,2},{1,1,0,2,1,1},{1,1,1,2,1,2},
        },//23
        { {0,0,0,0,0,0},
           {0,0,1,0,1,1},{1,1,0,2,0,2},{1,0,1,1,0,2},{0,1,1,1,0,2},{1,0,0,1,1,1},{0,1,0,1,1,1},
           {1,1,1,2,1,3},{1,0,1,1,1,2},{0,1,1,1,1,2},{1,1,1,2,0,3},{1,1,0,2,1,2},{1,1,1,2,1,3},
        },//24
        { {0,0,0,0,0,0},
           {1,0,1,1,0,2},{0,1,0,1,1,1},{0,0,1,0,1,1},{0,1,1,1,0,2},{1,0,0,1,1,1},{1,1,0,2,0,2},
           {1,1,1,2,1,3},{1,0,1,1,1,2},{1,1,1,2,0,3},{0,1,1,1,1,2},{1,1,0,2,1,2},{1,1,1,2,1,3},
        },//25
        { {0,0,0,0,0,0},
           {0,1,1,1,0,2},{1,0,0,1,1,1},{0,0,1,0,1,1},{1,0,1,1,0,2},{0,1,0,1,1,1},{1,1,0,2,0,2},
           {1,1,1,2,1,3},{0,1,1,1,1,2},{1,1,1,2,0,3},{1,0,1,1,1,2},{1,1,0,2,1,2},{1,1,1,2,1,3},
        },//26
        { {0,0,0,0,0,0},
           {1,0,0,1,1,1},{0,1,1,1,0,0},{1,0,1,1,0,2},{1,1,0,2,0,2},{0,0,1,0,1,1},{0,1,0,1,1,1},
           {1,1,1,2,1,1},{1,0,1,1,1,2},{1,1,0,2,1,2},{1,1,1,2,0,2},{0,1,1,1,1,1},{1,1,1,2,1,2},
        },//27
        { {0,0,0,0,0,0},
           {1,0,0,1,1,1},{0,1,1,1,0,2},{1,0,1,1,0,2},{1,1,0,2,0,2},{0,0,1,0,1,1},{0,1,0,1,1,1},
           {1,1,1,2,1,3},{1,0,1,1,1,2},{1,1,0,2,1,2},{1,1,1,2,0,3},{0,1,1,1,1,2},{1,1,1,2,1,3},
        },//28
        { {0,0,0,0,0,0},
           {1,0,1,1,0,2},{0,1,0,1,1,1},{1,0,0,1,1,1},{1,1,0,2,0,2},{0,0,1,0,1,1},{0,1,1,1,0,2},
           {1,1,1,2,1,3},{1,0,1,1,1,2},{1,1,1,2,0,3},{1,1,0,2,1,2},{0,1,1,1,1,2},{1,1,1,2,1,3},
        },//29
        { {0,0,0,0,0,0},
           {1,1,0,2,0,2},{0,0,1,0,1,1},{1,0,0,1,1,1},{1,0,1,1,0,2},{0,1,0,1,1,1},{0,1,1,1,0,2},
           {1,1,1,2,1,3},{1,1,0,2,1,2},{1,1,1,2,0,3},{1,0,1,1,1,2},{0,1,1,1,1,2},{1,1,1,2,1,3},
        },//30
        { {0,0,0,0,0,0},
           {0,1,0,1,1,1},{1,0,1,1,0,0},{0,1,1,1,0,2},{1,1,0,2,0,2},{0,0,1,0,1,1},{1,0,0,1,1,1},
           {1,1,1,2,1,1},{0,1,1,1,1,2},{1,1,0,2,1,2},{1,1,1,2,0,2},{1,0,1,1,1,1},{1,1,1,2,1,2},
        },//31
        { {0,0,0,0,0,0},
           {0,1,0,1,1,1},{1,0,1,1,0,2},{0,1,1,1,0,2},{1,1,0,2,0,2},{0,0,1,0,1,1},{1,0,0,1,1,1},
           {1,1,1,2,1,3},{0,1,1,1,1,2},{1,1,0,2,1,2},{1,1,1,2,0,3},{1,0,1,1,1,2},{1,1,1,2,1,3},
        },//32
        { {0,0,0,0,0,0},
           {0,1,1,1,0,2},{1,0,0,1,1,1},{0,1,0,1,1,1},{1,1,0,2,0,2},{0,0,1,0,1,1},{1,0,1,1,0,2},
           {1,1,1,2,1,3},{0,1,1,1,1,2},{1,1,1,2,0,3},{1,1,0,2,1,2},{1,0,1,1,1,2},{1,1,1,2,1,3},
        },//33
        { {0,0,0,0,0,0},
           {1,1,0,2,0,2},{0,0,1,0,1,1},{0,1,0,1,1,1},{0,1,1,1,0,2},{1,0,0,1,1,1},{1,0,1,1,0,2},
           {1,1,1,2,1,3},{1,1,0,2,1,2},{1,1,1,2,0,3},{0,1,1,1,1,2},{1,0,1,1,1,2},{1,1,1,2,1,3},
        } //34
    };
    
   std::vector<int> integralIndex{1,2,3,4,2,3,4,2,3,4,5,1,6,2,3,6,2,3,6,2,3,1,2,3,3,3,2,3,3,3,2,3,3,3};
   std::vector<int> permIndex{0,0,0,0,1,1,1,2,2,2,0,0,0,0,0,1,1,1,2,2,2,3,3,3,4,5,6,6,7,8,9,9,10,11};
 
    for(int m=0; m<15; ++m){
        for(int i=0; i<34; ++i){
            for(int j=0; j<13; ++j){
                int coefMult {};
                int index1 {permMat[m][permIndex[i]]-1};
                coefMult = coefMat[index1][j];
                collectPattern r{};
                r.patternID=m;
                r.integralID=integralIndex[i];
                int g1{};
                g1=g1Mat[g1Index[i]][j];
                int g2{};
                g2=g2Mat[g2Index[i]][j];
                int g3{};
                g3=g3Mat[g3Index[i]][j];
                if(g1==2){
                    r.gIndex=1;
                }
                if(g2==2){
                    r.gIndex=2;
                }
                if(g3==2){
                    r.gIndex=3;
                }
                if(g1==1 && g2==2){
                    r.gIndex=4;
                }
                if(g1==1 && g3==2){
                    r.gIndex=5;
                }
                if(g2==1 && g3==2){
                    r.gIndex=6;
                }
                if(g1==1 && g2==1 && g3==2){
                    r.gIndex=7;
                }
                if(g1==2 && g2==2 ){
                    r.gIndex=8;
                }
                if(g1==2 && g3==2){
                    r.gIndex=9;
                }
                for(int k=0; k<6; ++k){
                    r.hIndex[k]=g0Mat[i][j][k];
                }
                r.count=coefMult;
                s.push_back(r);
            }
        }
    }
    
    
    return s;
}



int main(int argc, const char * argv[]) {

    std::vector<collectPattern> a{};
    
  //  std::vector<collectPattern> symProbVec{getProbSym(1.0, 2.0, 3.0, 1.1, 1.0, 1.0, 1.0, 1.0, 1.0, 0.003,a)};
     std::vector<collectPattern> asymmProbVec{getProbAsymm(1.0, 2.00, 3.0, 1.1, 1.0, 1.0, 1.2, 1.0, 1.1, 0.003, a)};
    
        std::cout << asymmProbVec.size() << '\n';
    
    std::vector<collectPattern> b{consolidatePattern(asymmProbVec)};
    
    for(int i=0; i<146; ++i){
        std::cout << b[i].patternID << " " << b[i].integralID << " " << b[i].gIndex <<
        " ";
        for(int j=0; j<6; ++j){
            std::cout << b[i].hIndex[j] << " ";
        }
        std::cout << b[i].count << '\n';
    }
    std::cout << b.size() << '\n';
    
    int kMatrix[133][15]{};
    for(int i=0; i<b.size(); ++i){
        if(b[i].integralID==1 && b[i].gIndex==1){
            kMatrix[0][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==1 && b[i].hIndex[0]==1 && b[i].hIndex[1]==1){
            kMatrix[1][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==1 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[2][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==1 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[3][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==1 && b[i].hIndex[0]==1 && b[i].hIndex[1]==1){
            kMatrix[4][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==1 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[5][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==1 && b[i].hIndex[0]==1 && b[i].hIndex[4]==1){
            kMatrix[6][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==1 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[7][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==1 && b[i].hIndex[1]==1 && b[i].hIndex[4]==1){
            kMatrix[8][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==1 && b[i].hIndex[2]==1 && b[i].hIndex[4]==1){
            kMatrix[9][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==1 && b[i].hIndex[0]==1 && b[i].hIndex[1]==1){
            kMatrix[10][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==1 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[11][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==1 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[12][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==1){
            kMatrix[13][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==1 && b[i].hIndex[0]==1 && b[i].hIndex[1]==1){
            kMatrix[14][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==1 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[15][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==1 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[16][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[17][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[4]==1){
            kMatrix[18][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==2 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[19][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==2 && b[i].hIndex[1]==1 && b[i].hIndex[4]==1){
            kMatrix[20][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==2 && b[i].hIndex[2]==1 && b[i].hIndex[4]==1){
            kMatrix[21][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[1]==1){
            kMatrix[22][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[23][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[4]==1){
            kMatrix[24][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==2 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[25][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==2 && b[i].hIndex[1]==1 && b[i].hIndex[4]==1){
            kMatrix[26][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==2 && b[i].hIndex[2]==1 && b[i].hIndex[4]==1){
            kMatrix[27][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[1]==1){
            kMatrix[28][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[29][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[4]==1){
            kMatrix[30][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==2 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[31][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==2 && b[i].hIndex[1]==1 && b[i].hIndex[4]==1){
            kMatrix[32][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==2 && b[i].hIndex[2]==1 && b[i].hIndex[4]==1){
            kMatrix[33][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[4]==1){
            kMatrix[34][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==2 && b[i].hIndex[1]==1 && b[i].hIndex[4]==1){
            kMatrix[35][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==2 && b[i].hIndex[2]==1 && b[i].hIndex[4]==1){
            kMatrix[36][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[37][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==2 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[38][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[1]==1){
            kMatrix[39][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[40][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==2 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[41][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[42][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[4]==1){
            kMatrix[43][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[44][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[4]==1){
            kMatrix[45][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==3 && b[i].hIndex[2]==1 && b[i].hIndex[4]==1){
            kMatrix[46][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[1]==1){
            kMatrix[47][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[48][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[4]==1){
            kMatrix[49][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[50][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[4]==1){
            kMatrix[51][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==3 && b[i].hIndex[2]==1 && b[i].hIndex[4]==1){
            kMatrix[52][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[1]==1){
            kMatrix[53][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[54][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[4]==1){
            kMatrix[55][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[56][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[4]==1){
            kMatrix[57][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==3 && b[i].hIndex[2]==1 && b[i].hIndex[4]==1){
            kMatrix[58][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4&& b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[1]==1){
            kMatrix[59][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[60][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[4]==1){
            kMatrix[61][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[62][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[4]==1){
            kMatrix[63][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==3 && b[i].hIndex[2]==1 && b[i].hIndex[4]==1){
            kMatrix[64][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[4]==1){
            kMatrix[65][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[4]==1){
            kMatrix[66][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==3 && b[i].hIndex[2]==1 && b[i].hIndex[4]==1){
            kMatrix[67][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[4]==1){
            kMatrix[68][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[4]==1){
            kMatrix[69][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==3 && b[i].hIndex[2]==1 && b[i].hIndex[4]==1){
            kMatrix[70][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==4 && b[i].hIndex[4]==0 ){
            kMatrix[71][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==4 && b[i].hIndex[2]==0 ){
            kMatrix[72][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==4 && b[i].hIndex[4]==0 ){
            kMatrix[73][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==4 && b[i].hIndex[2]==0 ){
            kMatrix[74][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==4 && b[i].hIndex[1]==0 ){
            kMatrix[75][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==4 && b[i].hIndex[0]==0 ){
            kMatrix[76][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==4 && b[i].hIndex[4]==0 ){
            kMatrix[77][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==4 && b[i].hIndex[2]==0 ){
            kMatrix[78][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==4 && b[i].hIndex[1]==0 ){
            kMatrix[79][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==4 && b[i].hIndex[0]==0 ){
            kMatrix[80][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==4  ){
            kMatrix[81][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==4  ){
            kMatrix[82][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==5 && b[i].hIndex[4]==0 ){
            kMatrix[83][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==5 && b[i].hIndex[2]==0 ){
            kMatrix[84][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==5 && b[i].hIndex[4]==0 ){
            kMatrix[85][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==5 && b[i].hIndex[2]==0 ){
            kMatrix[86][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==5 && b[i].hIndex[1]==0 ){
            kMatrix[87][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==5 && b[i].hIndex[0]==0 ){
            kMatrix[88][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==5 && b[i].hIndex[4]==0 ){
            kMatrix[89][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==5 && b[i].hIndex[2]==0 ){
            kMatrix[90][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==5 && b[i].hIndex[1]==0 ){
            kMatrix[91][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==5 && b[i].hIndex[0]==0 ){
            kMatrix[92][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==5 && b[i].hIndex[4]==0 ){
            kMatrix[93][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==5 && b[i].hIndex[2]==0 ){
            kMatrix[94][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==5 && b[i].hIndex[1]==0 ){
            kMatrix[95][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==5 && b[i].hIndex[0]==0 ){
            kMatrix[96][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==5  ){
            kMatrix[97][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==5 && b[i].hIndex[2]==0 ){
            kMatrix[98][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==5 && b[i].hIndex[1]==0 ){
            kMatrix[99][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==5 && b[i].hIndex[0]==0 ){
            kMatrix[100][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==6 && b[i].hIndex[1]==0 ){
            kMatrix[101][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==6 && b[i].hIndex[0]==0 ){
            kMatrix[102][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==6 && b[i].hIndex[2]==0 ){
            kMatrix[103][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==6 && b[i].hIndex[1]==0 ){
            kMatrix[104][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==6 && b[i].hIndex[0]==0 ){
            kMatrix[105][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==6 && b[i].hIndex[4]==0 ){
            kMatrix[106][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==6 && b[i].hIndex[2]==0 ){
            kMatrix[107][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==6 && b[i].hIndex[1]==0 ){
            kMatrix[108][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==6 && b[i].hIndex[0]==0 ){
            kMatrix[109][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==6 && b[i].hIndex[2]==0 ){
            kMatrix[110][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==6 && b[i].hIndex[1]==0 ){
            kMatrix[111][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==6 && b[i].hIndex[0]==0 ){
            kMatrix[112][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==6 && b[i].hIndex[1]==0 ){
            kMatrix[113][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==6 && b[i].hIndex[0]==0 ){
            kMatrix[114][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==6 && b[i].hIndex[2]==0 ){
            kMatrix[115][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==6 && b[i].hIndex[1]==0 ){
            kMatrix[116][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==6 && b[i].hIndex[0]==0 ){
            kMatrix[117][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==7  ){
            kMatrix[118][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==7  ){
            kMatrix[119][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==7  ){
            kMatrix[120][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==7  ){
            kMatrix[121][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==7  ){
            kMatrix[122][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==7  ){
            kMatrix[123][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==8  ){
            kMatrix[124][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==8  ){
            kMatrix[125][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==8  ){
            kMatrix[126][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==8  ){
            kMatrix[127][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==9  ){
            kMatrix[128][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==9  ){
            kMatrix[129][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==9  ){
            kMatrix[130][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==9  ){
            kMatrix[131][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==9  ){
            kMatrix[132][b[i].patternID]=b[i].count;
        }
    }
    
    
   /* int kMatrix[116][15]{};
    for(int i=0; i<b.size(); ++i){
        if(b[i].integralID==1 && b[i].gIndex==1){
            kMatrix[0][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==1){
            kMatrix[1][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==1){
            kMatrix[2][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==1 && b[i].hIndex[0]==1 && b[i].hIndex[1]==1){
            kMatrix[3][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==1 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[4][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==1 && b[i].hIndex[0]==1 && b[i].hIndex[3]==1){
            kMatrix[5][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==1 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[6][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==1 && b[i].hIndex[1]==1 && b[i].hIndex[3]==1){
            kMatrix[7][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==1 && b[i].hIndex[2]==1 && b[i].hIndex[3]==1){
            kMatrix[8][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==1 && b[i].hIndex[0]==1 && b[i].hIndex[1]==1){
            kMatrix[9][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==1 && b[i].hIndex[0]==1 && b[i].hIndex[3]==1){
            kMatrix[10][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==1 && b[i].hIndex[1]==1 && b[i].hIndex[3]==1){
            kMatrix[11][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==1){
            kMatrix[12][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==2){
            kMatrix[13][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==2 && b[i].hIndex[2]==1 && b[i].hIndex[3]==1){
            kMatrix[14][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[15][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[3]==1){
            kMatrix[16][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==2 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[17][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==2 && b[i].hIndex[1]==1 && b[i].hIndex[3]==1){
            kMatrix[18][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==2){
            kMatrix[19][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[1]==1){
            kMatrix[20][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[21][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[3]==1){
            kMatrix[22][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==2 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[23][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==2 && b[i].hIndex[1]==1 && b[i].hIndex[3]==1){
            kMatrix[24][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==2 && b[i].hIndex[2]==1 && b[i].hIndex[3]==1){
            kMatrix[25][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==2 && b[i].hIndex[2]==1 && b[i].hIndex[3]==1){
            kMatrix[26][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[27][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==2 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[28][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[29][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==2 && b[i].hIndex[0]==1 && b[i].hIndex[3]==1){
            kMatrix[30][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==2 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[31][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==2 && b[i].hIndex[1]==1 && b[i].hIndex[3]==1){
            kMatrix[32][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[33][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[3]==1){
            kMatrix[34][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[35][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[3]==1){
            kMatrix[36][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[37][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[3]==1){
            kMatrix[38][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[39][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[3]==1){
            kMatrix[40][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==3 && b[i].hIndex[2]==1 && b[i].hIndex[3]==1){
            kMatrix[41][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[42][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[3]==1){
            kMatrix[43][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[44][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[3]==1){
            kMatrix[45][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[46][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[3]==1){
            kMatrix[47][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[48][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[3]==1){
            kMatrix[49][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[1]==1){
            kMatrix[50][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==3 && b[i].hIndex[2]==1 && b[i].hIndex[3]==1){
            kMatrix[51][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[52][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[3]==1){
            kMatrix[53][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[54][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[3]==1){
            kMatrix[55][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[1]==1){
            kMatrix[56][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==3 && b[i].hIndex[2]==1 && b[i].hIndex[3]==1){
            kMatrix[57][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[2]==1){
            kMatrix[58][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[3]==1){
            kMatrix[59][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[2]==1){
            kMatrix[60][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==3 && b[i].hIndex[1]==1 && b[i].hIndex[3]==1){
            kMatrix[61][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==3 && b[i].hIndex[0]==1 && b[i].hIndex[1]==1){
            kMatrix[62][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==4 && b[i].hIndex[3]==0 ){
            kMatrix[63][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==4 && b[i].hIndex[2]==0 ){
            kMatrix[64][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==4 && b[i].hIndex[3]==0 ){
            kMatrix[65][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==4 && b[i].hIndex[2]==0 ){
            kMatrix[66][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==4 && b[i].hIndex[1]==0 ){
            kMatrix[67][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==4 && b[i].hIndex[0]==0 ){
            kMatrix[68][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==4 && b[i].hIndex[1]==0 ){
            kMatrix[69][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==4 && b[i].hIndex[0]==0 ){
            kMatrix[70][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==5 && b[i].hIndex[3]==0 ){
            kMatrix[71][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==5 && b[i].hIndex[2]==0 ){
            kMatrix[72][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==5 && b[i].hIndex[3]==0 ){
            kMatrix[73][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==5 && b[i].hIndex[2]==0 ){
            kMatrix[74][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==5 && b[i].hIndex[3]==0 ){
            kMatrix[75][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==5 && b[i].hIndex[2]==0 ){
            kMatrix[76][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==5 && b[i].hIndex[3]==0 ){
            kMatrix[77][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==5 && b[i].hIndex[2]==0 ){
            kMatrix[78][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==5 && b[i].hIndex[1]==0 ){
            kMatrix[79][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==5 && b[i].hIndex[0]==0 ){
            kMatrix[80][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==5 && b[i].hIndex[3]==0 ){
            kMatrix[81][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==5 && b[i].hIndex[2]==0 ){
            kMatrix[82][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==5 && b[i].hIndex[1]==0 ){
            kMatrix[83][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==5 && b[i].hIndex[0]==0 ){
            kMatrix[84][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==5 && b[i].hIndex[1]==0 ){
            kMatrix[85][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==5 && b[i].hIndex[0]==0 ){
            kMatrix[86][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==6 && b[i].hIndex[1]==0 ){
            kMatrix[87][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==6 && b[i].hIndex[0]==0 ){
            kMatrix[88][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==6 && b[i].hIndex[1]==0 ){
            kMatrix[89][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==6 && b[i].hIndex[0]==0 ){
            kMatrix[90][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==6 && b[i].hIndex[1]==0 ){
            kMatrix[91][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==6 && b[i].hIndex[0]==0 ){
            kMatrix[92][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==6 && b[i].hIndex[3]==0 ){
            kMatrix[93][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==6 && b[i].hIndex[2]==0 ){
            kMatrix[94][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==6 && b[i].hIndex[1]==0 ){
            kMatrix[95][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==6 && b[i].hIndex[0]==0 ){
            kMatrix[96][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==6 && b[i].hIndex[3]==0 ){
            kMatrix[97][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==6 && b[i].hIndex[1]==0 ){
            kMatrix[98][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==6 && b[i].hIndex[0]==0 ){
            kMatrix[99][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==6 && b[i].hIndex[3]==0 ){
            kMatrix[100][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==6 && b[i].hIndex[2]==0 ){
            kMatrix[101][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==7  ){
            kMatrix[102][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==7  ){
            kMatrix[103][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==7  ){
            kMatrix[104][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==7  ){
            kMatrix[105][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==7  ){
            kMatrix[106][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==7  ){
            kMatrix[107][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==1 && b[i].gIndex==8  ){
            kMatrix[108][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==8  ){
            kMatrix[109][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==3 && b[i].gIndex==8  ){
            kMatrix[110][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==8  ){
            kMatrix[111][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==5 && b[i].gIndex==8  ){
            kMatrix[112][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==2 && b[i].gIndex==9  ){
            kMatrix[113][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==4 && b[i].gIndex==9  ){
            kMatrix[114][b[i].patternID]=b[i].count;
        }
        if(b[i].integralID==6 && b[i].gIndex==9  ){
            kMatrix[115][b[i].patternID]=b[i].count;
        }
    }

    */
    
    for(int i=0; i<133; ++i){
        std::cout << (i+1) << " ";
        for(int j=0; j<15; ++j){
            std::cout << kMatrix[i][j] << " & ";
        }
        std::cout << '\n';
    }
    
    for(int i=0; i<133; ++i){
        int testSum{};
        testSum += kMatrix[i][0];
        for(int j=1; j<8; ++j){
            testSum += 3*kMatrix[i][j];
        }
        for(int j=8; j<15; ++j){
            testSum += 6*kMatrix[i][j];
        }
        if(testSum !=0){
            std::cout << i << '\n';
        }
    }
    std::cout << '\n';
   
    
    return 0;
}
