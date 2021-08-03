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
                                        std::string s{};
                                        if(node0==nodeA){
                                            s.append("r1*");
                                        }
                                        else{
                                            s.append("s1*");
                                        }
                                        if(node0==nodeB){
                                            s.append("r2*");
                                        }
                                        else{
                                            s.append("s2*");
                                        }
                                        if(node1==nodeC){
                                            s.append("r3*");
                                        }
                                        else{
                                            s.append("s3*");
                                        }
                                        if(node1==nodeD){
                                            s.append("r4*");
                                        }
                                        else{
                                            s.append("s4*");
                                        }
                                        if(node0==node1){
                                            s.append("r5");
                                        }
                                        else{
                                            s.append("s5");
                                        }
                                        patternCount r{};
                                        r.pattern=s;
                                        r.counts[index1]=1;
                                         sitePatternProbs.push_back(r);
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
    
    return 0;
}
