//
//  main.cpp
//  five taxa
//
//  Created by Andrew Richards on 12/23/19.
//  Copyright © 2019 Andrew Richards. All rights reserved.
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
#include <cstdlib>

struct patternCount{
    std::string pattern;
    int counts[15];
};

struct patternCount2{
    std::string pattern;
    int counts[32];
};

std::vector<patternCount> consolidateCounts(std::vector<patternCount> input1);


std::vector<patternCount> consolidateCounts(std::vector<patternCount> input1){
    std::vector<patternCount> x{};
    std::vector<int> consolidatedInd{};
    for(int i=0; i<input1.size();++i){
        consolidatedInd.push_back(0);
    }
    
    for(int i=0; i<input1.size();++i){
        if(consolidatedInd[i]==0){
            consolidatedInd[i]=1;
            patternCount s{};
            s.pattern=input1[i].pattern;
            for(int k=0; k<15; ++k){
                s.counts[k]=input1[i].counts[k];
            }
        for(int j=(i+1); j<input1.size();++j){
            if(consolidatedInd[j]==0){
                if(input1[i].pattern.compare(input1[j].pattern)==0){
                    consolidatedInd[j]=1;
                               for(int k=0; k<15; ++k){
                        s.counts[k] += input1[j].counts[k];
                    }
                }
            }
        }
              x.push_back(s);
        }
    }
    
    
    return x;
}

int main(int argc, const char * argv[]) {
    // insert code here...
   
    std::vector<std::string> sitePatterns{
      "XXXX", "XXXY", "XXYX", "XYXX", "YXXX", "XXYY", "XYXY", "XYYX",
        "XXYZ", "XYXZ", "XYZX", "YXXZ", "YXZX", "YZXX", "XYZW"
    };
    
   std::vector<patternCount> sitePatternProbs{};

    std::vector<std::vector<int>> coefMat{};
    std::vector<int> zeroVec{0,0,0,0,0,0};
    for(int i=0; i<32; ++i){
        coefMat.push_back(zeroVec);
    }
    for(int i=0; i<2; ++i){
        for(int j=0; j<2; ++j){
            for(int k=0; k<2; ++k){
                for(int l=0; l<2; ++l){
                    for(int m=0; m<2; ++m){
                        if(i==0){
                           coefMat[m+2*l+4*k+8*j+16*i][0]=3;
                        }
                        else{
                             coefMat[m+2*l+4*k+8*j+16*i][0]=-1;
                        }
                        if(j==0){
                           coefMat[m+2*l+4*k+8*j+16*i][1]=3;
                        }
                        else{
                             coefMat[m+2*l+4*k+8*j+16*i][1]=-1;
                        }
                        if(k==0){
                           coefMat[m+2*l+4*k+8*j+16*i][2]=3;
                        }
                        else{
                             coefMat[m+2*l+4*k+8*j+16*i][2]=-1;
                        }
                        if(l==0){
                           coefMat[m+2*l+4*k+8*j+16*i][3]=3;
                        }
                        else{
                             coefMat[m+2*l+4*k+8*j+16*i][3]=-1;
                        }
                        if(m==0){
                           coefMat[m+2*l+4*k+8*j+16*i][4]=3;
                        }
                        else{
                             coefMat[m+2*l+4*k+8*j+16*i][4]=-1;
                        }
                    }
                }
            }
        }
    }
    

    
    std::vector<std::vector<std::vector<int>>> coefArray{};
    for(int i=0; i<15; ++i){
        coefArray.push_back(coefMat);
    }
    
    for(int nodeA=0; nodeA<4; ++nodeA){
        for(int nodeB=0; nodeB<4; ++nodeB){
            for(int nodeC=0; nodeC<4; ++nodeC){
                for(int nodeD=0; nodeD<4; ++nodeD){
                        int index1{};
                        if(nodeA == nodeB && nodeA == nodeC && nodeA == nodeD){
                            index1=0;
                        }
                        if(nodeA == nodeB && nodeA == nodeC && nodeA != nodeD ){
                            index1=1;
                        }
                        if(nodeA == nodeB && nodeA == nodeD  && nodeA != nodeC){
                            index1=2;
                        }
                        if(nodeA != nodeB && nodeA == nodeD  && nodeA == nodeC){
                            index1=3;
                        }
                        if(nodeA != nodeB && nodeB == nodeC  && nodeB == nodeD){
                            index1=4;
                        }
                        if(nodeA==nodeB && nodeC==nodeD && nodeA != nodeC){
                            index1=5;
                        }
                        if(nodeA==nodeC && nodeB==nodeD && nodeA != nodeB ){
                            index1=6;
                        }
                        if(nodeA==nodeD && nodeC==nodeB && nodeA != nodeC){
                            index1=7;
                        }
                        if(nodeA==nodeB && nodeA != nodeC && nodeA != nodeD && nodeC != nodeD){
                            index1=8;
                        }
                        if(nodeA==nodeC && nodeA != nodeB && nodeA != nodeD && nodeB != nodeD){
                            index1=9;
                        }
                        if(nodeA==nodeD && nodeA != nodeB && nodeA != nodeC && nodeB != nodeC){
                            index1=10;
                        }
                        if(nodeB==nodeC && nodeA != nodeB && nodeA != nodeD && nodeB != nodeD){
                            index1=11;
                        }
                        if(nodeB==nodeD && nodeA != nodeC && nodeA != nodeB && nodeB != nodeC){
                            index1=12;
                        }
                        if(nodeC==nodeD && nodeA != nodeB && nodeA != nodeC && nodeB != nodeC){
                            index1=13;
                        }
                        if(nodeA != nodeB && nodeA != nodeC && nodeA != nodeD && nodeB != nodeC &&
                           nodeB != nodeD && nodeC != nodeD ){
                            index1=14;
                        }
                      
                       for(int node0=0; node0<4; ++node0){
                            for(int node1=0; node1<4; ++node1){
                                                            int rowIndex{};
                                        if(node0 != nodeA){
                                            rowIndex += 16;
                                        }

                                        if(node0 != nodeB){
                                            rowIndex += 8;
                                        }

                                        if(node1 != nodeC){
                                            rowIndex += 4;
                                        }

                                        if(node1 != nodeD){
                                            rowIndex += 2;
                                        }

                                        if(node0 != node1){
                                           rowIndex += 1;
                                        }
                                ++coefArray[index1][rowIndex][5];
                                
                                    }
                                }
                    
                    
                    
                            }
                        }
                    }
                }

    
    
   
    std::vector<patternCount> output1{consolidateCounts(sitePatternProbs)};
    
    std::cout<<sitePatternProbs.size() << " " << output1.size() <<'\n';
    
    
    for(int i=0; i<15; ++i){
        int firstInd{1};
        std::cout << sitePatterns[i] << " ";
        for(int j=0; j<output1.size(); ++j){
            if(output1[j].counts[i] !=0){
                if(firstInd==0){
                    std::cout << "+";
                }
                firstInd=0;
                std::cout <<output1[j].counts[i] << "*" << output1[j].pattern;
            }
        }
        std::cout << '\n';
    }
    
/*
    
    for(int i=0; i<32; ++i){
        for(int j=0; j<6; ++j){
            std::cout << coefArray[0][i][j] << " ";
        }
        std::cout << '\n';
    }

*/
    std::vector<patternCount2> finalCoef{};
    for(int i=0; i<15; ++i){
        patternCount2 r{};
        r.pattern=sitePatterns[i];
        finalCoef.push_back(r);
    }
    
    for(int i=0; i<15; ++i){
        for(int j=0; j<32; ++j){
            finalCoef[i].counts[j]=0;
        }
        for(int k=0; k<32; ++k){
                finalCoef[i].counts[0] += coefArray[i][k][5];
                finalCoef[i].counts[1] += coefArray[i][k][0]*coefArray[i][k][5];
                finalCoef[i].counts[2] += coefArray[i][k][1]*coefArray[i][k][5];
                finalCoef[i].counts[3] += coefArray[i][k][2]*coefArray[i][k][5];
                finalCoef[i].counts[4] += coefArray[i][k][3]*coefArray[i][k][5];
                finalCoef[i].counts[5] += coefArray[i][k][4]*coefArray[i][k][5];
                finalCoef[i].counts[6] += coefArray[i][k][0]*coefArray[i][k][1]*coefArray[i][k][5];
                finalCoef[i].counts[7] += coefArray[i][k][0]*coefArray[i][k][2]*coefArray[i][k][5];
                finalCoef[i].counts[8] += coefArray[i][k][0]*coefArray[i][k][3]*coefArray[i][k][5];
                finalCoef[i].counts[9] += coefArray[i][k][0]*coefArray[i][k][4]*coefArray[i][k][5];
                finalCoef[i].counts[10] += coefArray[i][k][1]*coefArray[i][k][2]*coefArray[i][k][5];
                finalCoef[i].counts[11] += coefArray[i][k][1]*coefArray[i][k][3]*coefArray[i][k][5];
                finalCoef[i].counts[12] += coefArray[i][k][1]*coefArray[i][k][4]*coefArray[i][k][5];
                finalCoef[i].counts[13] += coefArray[i][k][2]*coefArray[i][k][3]*coefArray[i][k][5];
                finalCoef[i].counts[14] += coefArray[i][k][2]*coefArray[i][k][4]*coefArray[i][k][5];
                finalCoef[i].counts[15] += coefArray[i][k][3]*coefArray[i][k][4]*coefArray[i][k][5];
                        finalCoef[i].counts[16] += coefArray[i][k][0]*coefArray[i][k][1]*coefArray[i][k][2]*coefArray[i][k][5];
                        finalCoef[i].counts[17] += coefArray[i][k][0]*coefArray[i][k][1]*coefArray[i][k][3]*coefArray[i][k][5];
                        finalCoef[i].counts[18] += coefArray[i][k][0]*coefArray[i][k][1]*coefArray[i][k][4]*coefArray[i][k][5];
                        finalCoef[i].counts[19] += coefArray[i][k][0]*coefArray[i][k][2]*coefArray[i][k][3]*coefArray[i][k][5];
                        finalCoef[i].counts[20] += coefArray[i][k][0]*coefArray[i][k][2]*coefArray[i][k][4]*coefArray[i][k][5];
                        finalCoef[i].counts[21] += coefArray[i][k][0]*coefArray[i][k][3]*coefArray[i][k][4]*coefArray[i][k][5];
                        finalCoef[i].counts[22] += coefArray[i][k][1]*coefArray[i][k][2]*coefArray[i][k][3]*coefArray[i][k][5];
                        finalCoef[i].counts[23] += coefArray[i][k][1]*coefArray[i][k][2]*coefArray[i][k][4]*coefArray[i][k][5];
                        finalCoef[i].counts[24] += coefArray[i][k][1]*coefArray[i][k][3]*coefArray[i][k][4]*coefArray[i][k][5];
                        finalCoef[i].counts[25] += coefArray[i][k][2]*coefArray[i][k][3]*coefArray[i][k][4]*coefArray[i][k][5];
        finalCoef[i].counts[26] += coefArray[i][k][0]*coefArray[i][k][1]*coefArray[i][k][2]*coefArray[i][k][3]*coefArray[i][k][5];
        finalCoef[i].counts[27] += coefArray[i][k][0]*coefArray[i][k][1]*coefArray[i][k][2]*coefArray[i][k][4]*coefArray[i][k][5];
        finalCoef[i].counts[28] += coefArray[i][k][0]*coefArray[i][k][1]*coefArray[i][k][3]*coefArray[i][k][4]*coefArray[i][k][5];
        finalCoef[i].counts[29] += coefArray[i][k][0]*coefArray[i][k][2]*coefArray[i][k][3]*coefArray[i][k][4]*coefArray[i][k][5];
        finalCoef[i].counts[30] += coefArray[i][k][1]*coefArray[i][k][2]*coefArray[i][k][3]*coefArray[i][k][4]*coefArray[i][k][5];
                    finalCoef[i].counts[31] += coefArray[i][k][0]*coefArray[i][k][1]*coefArray[i][k][2]*
                          coefArray[i][k][3]*coefArray[i][k][4]*coefArray[i][k][5];
        }
    }
    
    for(int i=0; i<15; ++i){
        std::cout << finalCoef[i].pattern << " ";
        for(int j=0; j<32; ++j){
            std::cout << finalCoef[i].counts[j] << " ";
        }
        std::cout << '\n';
    }
    
    std::vector<int> coefTest{};
    for(int i=0; i<32; ++i){
        int a{};
        for(int j=0; j<15; ++j){
            if(finalCoef[j].counts[i] !=0){
                a=1;
            }
        }
        coefTest.push_back(a);
    }
    
   /* for(int i=0; i<32; ++i){
        std::cout << coefTest[i] << " ";
    }
    std::cout << '\n';
    
    for(int i=0; i<32; ++i){
        if(coefTest[i]==1){
            std::cout << i << " ";
        }
    }
        std::cout << '\n'; */
    
    return 0;
}
