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

typedef std::vector< std::vector<double> > doubleMat2;
typedef std::vector< std::vector<int> > intMat2;
typedef std::vector< doubleMat2 > doubleMat3;

long double asymmInt1 (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
long double asymmInt2sep (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
long double asymmInt2tog (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
long double asymmInt3sep (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
long double asymmInt3tog (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
long double asymmInt4sep (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
long double asymmInt4tog (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
long double asymmInt5sep (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
long double asymmInt5tog12 (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
long double asymmInt5tog23 (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
long double asymmInt5togAll (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
long double symInt1 (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
long double symInt2sep (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
long double symInt2tog (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
long double symInt3sep (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
long double symInt3tog (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
long double symInt4sep (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
long double symInt4tog12 (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
long double symInt4tog23 (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
long double symInt4togAll (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3);
std::vector<long double> probSym(double tau1, double tau2, double tau3, double gammaA, double gammaB, double gammaC, double gammaD, doubleMat2 ABmat, doubleMat2 CDmat, doubleMat2 rootMat, double theta);
std::vector<long double> probAsymm(double tau1, double tau2, double tau3, double gammaA, double gammaB, double gammaC, double gammaD, doubleMat2 CDmat, doubleMat2 BCDmat, doubleMat2 rootMat, double theta);

long double symInt1 (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{};  // both below root
    y = exp(tau1+tau2+tau3);
    y *= (exp(-psi1a*(1.0+g1))-exp(-psi1b*(1.0+g1)));
    y *= (exp(-psi2a*(1.0+g2))-exp(-psi2b*(1.0+g2)));
    y *= (exp(-psi3a*(1.0+g3))-exp(-psi3b*(1.0+g3)));
    y /= (1.0+g1);
    y /= (1.0+g2);
    y /= (1.0+g3);
    return y;
}

long double symInt2sep (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{}; // t2<tau3<t1
    y = exp(tau1+tau2+2.0*tau3);
    y *= (exp(-psi1a*(2.0+g1))-exp(-psi1b*(2.0+g1)));
    y *= (exp(-psi2a*(1.0+g2))-exp(-psi2b*(1.0+g2)));
    y *= (exp(-psi3a*(1.0+g3))-exp(-psi3b*(1.0+g3)));
    y /= (2.0+g1);
    y /= (1.0+g2);
    y /= (1.0+g3);
    return y;
}

long double symInt2tog (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{}; // t2<tau3<t1
    y = exp(tau1+tau2+2.0*tau3);
    y *= (exp(-psi2a*(1.0+g2))-exp(-psi2b*(1.0+g2)));
    y *= ((exp(-psi3a*(3.0+g1+g3))-exp(-psi3a*(2.0+g1)-psi3b*(1.0+g3)))/(1.0+g3)
          -(exp(-psi3a*(3.0+g1+g3))-exp(-psi3b*(3.0+g1+g3)))/(3.0+g1+g3));
    y /= (2.0+g1);
    y /= (1.0+g2);
    return y;
}

long double symInt3sep (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{}; // t1<tau3<t2
    y = exp(tau1+tau2+2.0*tau3);
    y *= (exp(-psi1a*(1.0+g1))-exp(-psi1b*(1.0+g1)));
    y *= (exp(-psi2a*(2.0+g2))-exp(-psi2b*(2.0+g2)));
    y *= (exp(-psi3a*(1.0+g3))-exp(-psi3b*(1.0+g3)));
    y /= (1.0+g1);
    y /= (2.0+g2);
    y /= (1.0+g3);
    return y;
}

long double symInt3tog (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{}; // t1<tau3<t2
    y = exp(tau1+tau2+2.0*tau3);
    y *= (exp(-psi1a*(1.0+g1))-exp(-psi1b*(1.0+g1)));
    y *= ((exp(-psi3a*(3.0+g2+g3))-exp(-psi3a*(2.0+g2)-psi3b*(1.0+g3)))/(1.0+g3)
          -(exp(-psi3a*(3.0+g2+g3))-exp(-psi3b*(3.0+g2+g3)))/(3.0+g2+g3));
    y /= (1.0+g1);
    y /= (2.0+g2);
    return y;
}

long double symInt4sep (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{}; //tau3<t1
    y = exp(tau1+tau2+4.0*tau3);
    y *= (exp(-psi1a*(3.0+g1))-exp(-psi1b*(3.0+g1)));
    y *= (exp(-psi2a*(2.0+g2))-exp(-psi2b*(2.0+g2)));
    y *= (exp(-psi3a*(1.0+g3))-exp(-psi3b*(1.0+g3)));
    y /= (3.0+g1);
    y /= (2.0+g2);
    y /= (1.0+g3);
    return y;
}

long double symInt4tog12 (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{}; // tau3<t1
    y = exp(tau1+tau2+4.0*tau3);
    y *= (exp(-psi3a*(1.0+g3))-exp(-psi3b*(1.0+g3)));
    y *= ((exp(-psi2a*(5.0+g1+g2))-exp(-psi2a*(3.0+g1)-psi2b*(2.0+g2)))/(2.0+g2)
          -(exp(-psi2a*(5.0+g1+g2))-exp(-psi2b*(5.0+g1+g2)))/(5.0+g1+g2));
    y /= (3.0+g1);
    y /= (1.0+g3);
    return y;
}

long double symInt4tog23 (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{}; // tau3<t1
    y = exp(tau1+tau2+4.0*tau3);
    y *= (exp(-psi1a*(3.0+g1))-exp(-psi1b*(3.0+g1)));
    y *= ((exp(-psi3a*(3.0+g2+g3))-exp(-psi3a*(2.0+g2)-psi3b*(1.0+g3)))/(1.0+g3)
          -(exp(-psi3a*(3.0+g2+g3))-exp(-psi3b*(3.0+g2+g3)))/(3.0+g2+g3));
    y /= (3.0+g1);
    y /= (2.0+g2);
    return y;
}

long double symInt4togAll (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{}; // tau3<t1
    y = exp(tau1+tau2+4.0*tau3);
    y *= (
        (exp(-psi3a*(6.0+g1+g2+g3))-exp(-psi3a*(5.0+g1+g2)-psi3b*(1.0+g3)))/((1.0+g3)*(2.0+g2))
        -(exp(-psi3a*(6.0+g1+g2+g3))-exp(-psi3a*(3.0+g1)-psi3b*(3.0+g2+g3)))/((3.0+g2+g3)*(2.0+g2))
        -(exp(-psi3a*(6.0+g1+g2+g3))-exp(-psi3a*(5.0+g1+g2)-psi3b*(1.0+g3)))/((5.0+g1+g2)*(1.0+g3))
        +(exp(-psi3a*(6.0+g1+g2+g3))-exp(-psi3b*(6.0+g1+g2+g3)))/((5.0+g1+g2)*(6.0+g1+g2+g3))
          );
    y /= (3.0+g1);
    return y;
}

long double asymmInt1 (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{};  //t1<tau2<t2<tau3
    y = exp(tau1+tau2+tau3);
    y *= (exp(-psi1a*(1.0+g1))-exp(-psi1b*(1.0+g1)));
    y *= (exp(-psi2a*(1.0+g2))-exp(-psi2b*(1.0+g2)));
    y *= (exp(-psi3a*(1.0+g3))-exp(-psi3b*(1.0+g3)));
    y /= (1.0+g1);
    y /= (1.0+g2);
    y /= (1.0+g3);
    return y;
}

long double asymmInt2sep (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{}; //t1<tau2<tau3<t2
    y = exp(tau1+tau2+2.0*tau3);
    y *= (exp(-psi1a*(1.0+g1))-exp(-psi1b*(1.0+g1)));
    y *= (exp(-psi2a*(2.0+g2))-exp(-psi2b*(2.0+g2)));
    y *= (exp(-psi3a*(1.0+g3))-exp(-psi3b*(1.0+g3)));
    y /= (1.0+g1);
    y /= (2.0+g2);
    y /= (1.0+g3);
    return y;
}

long double asymmInt2tog (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{}; //t1<tau2<tau3<t2
    y = exp(tau1+tau2+2.0*tau3);
    y *= (exp(-psi1a*(1.0+g1))-exp(-psi1b*(1.0+g1)));
    y *= ((exp(-psi3a*(3.0+g2+g3))-exp(-psi3a*(2.0+g2)-psi3b*(1.0+g3)))/(1.0+g3)
          -(exp(-psi3a*(3.0+g2+g3))-exp(-psi3b*(3.0+g2+g3)))/(3.0+g2+g3));
    y /= (1.0+g1);
    y /= (2.0+g2);
    return y;
}

long double asymmInt3sep (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{}; // tau2<t1<t2<tau3
    y = exp(tau1+2.0*tau2+tau3);
    y *= (exp(-psi1a*(2.0+g1))-exp(-psi1b*(2.0+g1)));
    y *= (exp(-psi2a*(1.0+g2))-exp(-psi2b*(1.0+g2)));
    y *= (exp(-psi3a*(1.0+g3))-exp(-psi3b*(1.0+g3)));
    y /= (2.0+g1);
    y /= (1.0+g2);
    y /= (1.0+g3);
    return y;
}

long double asymmInt3tog (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{}; // tau2<t1<t2<tau3
    y = exp(tau1+2.0*tau2+tau3);
    y *= (exp(-psi3a*(1.0+g3))-exp(-psi3b*(1.0+g3)));
    y *= ((exp(-psi2a*(3.0+g1+g2))-exp(-psi2a*(2.0+g1)-psi2b*(1.0+g2)))/(1.0+g2)
          -(exp(-psi2a*(3.0+g1+g2))-exp(-psi2b*(3.0+g1+g2)))/(3.0+g1+g2));
    y /= (2.0+g1);
    y /= (1.0+g3);
    return y;
}

long double asymmInt4sep (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{}; // tau2<t1<tau3<t2
    y = exp(tau1+2.0*tau2+2.0*tau3);
    y *= (exp(-psi1a*(2.0+g1))-exp(-psi1b*(2.0+g1)));
    y *= (exp(-psi2a*(2.0+g2))-exp(-psi2b*(2.0+g2)));
    y *= (exp(-psi3a*(1.0+g3))-exp(-psi3b*(1.0+g3)));
    y /= (2.0+g1);
    y /= (2.0+g2);
    y /= (1.0+g3);
    return y;
}

long double asymmInt4tog (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{}; // tau2<t1<tau3<t2
    y = exp(tau1+2.0*tau2+2.0*tau3);
    y *= (exp(-psi1a*(2.0+g1))-exp(-psi1b*(2.0+g1)));
    y *= ((exp(-psi3a*(3.0+g2+g3))-exp(-psi3a*(2.0+g2)-psi3b*(1.0+g3)))/(1.0+g3)
          -(exp(-psi3a*(3.0+g2+g3))-exp(-psi3b*(3.0+g2+g3)))/(3.0+g2+g3));
    y /= (2.0+g1);
    y /= (2.0+g2);
    return y;
}

long double asymmInt5sep (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{}; // tau3<t1
    y = exp(tau1+2.0*tau2+3.0*tau3);
    y *= (exp(-psi1a*(3.0+g1))-exp(-psi1b*(3.0+g1)));
    y *= (exp(-psi2a*(2.0+g2))-exp(-psi2b*(2.0+g2)));
    y *= (exp(-psi3a*(1.0+g3))-exp(-psi3b*(1.0+g3)));
    y /= (3.0+g1);
    y /= (2.0+g2);
    y /= (1.0+g3);
    return y;
}

long double asymmInt5tog12 (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{}; // tau3<t1
    y = exp(tau1+2.0*tau2+3.0*tau3);
    y *= (exp(-psi3a*(1.0+g3))-exp(-psi3b*(1.0+g3)));
    y *= ((exp(-psi2a*(5.0+g1+g2))-exp(-psi2a*(3.0+g1)-psi2b*(2.0+g2)))/(2.0+g2)
          -(exp(-psi2a*(5.0+g1+g2))-exp(-psi2b*(5.0+g1+g2)))/(5.0+g1+g2));
    y /= (3.0+g1);
    y /= (1.0+g3);
    return y;
}

long double asymmInt5tog23 (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{}; // tau3<t1
    y = exp(tau1+2.0*tau2+3.0*tau3);
    y *= (exp(-psi1a*(3.0+g1))-exp(-psi1b*(3.0+g1)));
    y *= ((exp(-psi3a*(3.0+g2+g3))-exp(-psi3a*(2.0+g2)-psi3b*(1.0+g3)))/(1.0+g3)
          -(exp(-psi3a*(3.0+g2+g3))-exp(-psi3b*(3.0+g2+g3)))/(3.0+g2+g3));
    y /= (3.0+g1);
    y /= (2.0+g2);
    return y;
}

long double asymmInt5togAll (double tau1, double tau2, double tau3, double psi1a, double psi1b, double psi2a, double psi2b, double psi3a, double psi3b, long double g1, long double g2, long double g3){
    long double y{}; // tau3<t1
    y = exp(tau1+2.0*tau2+3.0*tau3);
    y *= (
          (exp(-psi3a*(6.0+g1+g2+g3))-exp(-psi3a*(5.0+g1+g2)-psi3b*(1.0+g3)))/((1.0+g3)*(2.0+g2))
          -(exp(-psi3a*(6.0+g1+g2+g3))-exp(-psi3a*(3.0+g1)-psi3b*(3.0+g2+g3)))/((3.0+g2+g3)*(2.0+g2))
          -(exp(-psi3a*(6.0+g1+g2+g3))-exp(-psi3a*(5.0+g1+g2)-psi3b*(1.0+g3)))/((5.0+g1+g2)*(1.0+g3))
          +(exp(-psi3a*(6.0+g1+g2+g3))-exp(-psi3b*(6.0+g1+g2+g3)))/((5.0+g1+g2)*(6.0+g1+g2+g3))
          );
    y /= (3.0+g1);
    return y;
}

std::vector<long double> probSym(double tau1, double tau2, double tau3, double gammaA, double gammaB, double gammaC, double gammaD, doubleMat2 ABmat, doubleMat2 CDmat, doubleMat2 rootMat, double theta){
    std::vector<long double> probs{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double thetaNew {4.0*theta/6.0};
    intMat2 coefMat{
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
    intMat2 permMat{
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
    intMat2 distMat{
        {0,0,0,0,0},
        {1,1,0,0,0},
        {0,0,1,1,0},
        {1,0,1,0,1},
        {1,0,0,1,1},
        {0,1,1,0,1},
        {0,1,0,1,1},
        {1,1,1,1,0},
        {1,1,1,0,1},
        {1,1,0,1,1},
        {1,0,1,1,1},
        {0,1,1,1,1},
        {1,1,1,1,1}
    };
      
    double g0{};
    double g1{};
    double g2{};
    double g3{};
    double d1{};
    double d2{};
    double d3{};
    double d4{};
    double d5{};
    double psi1a{};
    double psi1b{};
    double psi2a{};
    double psi2b{};
    double psi3a{};
    double psi3b{};
    doubleMat2 CDmatFull{CDmat};
    doubleMat2 ABmatFull{ABmat};
    for(int i=0; i<rootMat.size(); ++i){
        CDmatFull.push_back(rootMat[i]);
        ABmatFull.push_back(rootMat[i]);
    }
    
    
    int sizeAB{static_cast<int>(ABmat.size())};
    int sizeCD{static_cast<int>(CDmat.size())};
    
    std::vector< std::vector<long double> > intMat{};
    std::vector<long double> intVec{};
    
    // ((A,B),(C,D))

    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
         psi3a = rootMat[t3Loc][0];
         psi3b = rootMat[t3Loc][1];
        for(int t1Loc=0; (t1Loc-sizeCD)<=t3Loc; ++t1Loc){
            psi1a = CDmatFull[t1Loc][0];
            psi1b = CDmatFull[t1Loc][1];
            for(int t2Loc=0; (t2Loc-sizeAB)<=t3Loc; ++t2Loc){
                psi2a = ABmatFull[t2Loc][0];
                psi2b = ABmatFull[t2Loc][1];
                d3 = tau1*(gammaC-CDmat[0][2]);
                d4 = tau1*(gammaD-CDmat[0][2]);
                d1 = tau2*(gammaA-ABmat[0][2]);
                d2 = tau2*(gammaB-ABmat[0][2]);
                d5 = 0;
                if(t1Loc !=0){
                    for(int i=1; i<=t1Loc; ++i){
                        d3 += CDmatFull[i][0]*(CDmatFull[i-1][2]-CDmatFull[i][2]);
                        d4 += CDmatFull[i][0]*(CDmatFull[i-1][2]-CDmatFull[i][2]);
                    }
                }
                if(t2Loc !=0){
                    for(int i=1; i<=t2Loc; ++i){
                    d1 += ABmatFull[i][0]*(ABmatFull[i-1][2]-ABmatFull[i][2]);
                    d2 += ABmatFull[i][0]*(ABmatFull[i-1][2]-ABmatFull[i][2]);
                    }
                }
                if((t1Loc-sizeCD) < t3Loc){
                    for(int i=(t1Loc+1); i<=(t3Loc+sizeCD); ++i){
                        d5 += CDmatFull[i][0]*(CDmatFull[i-1][2]-CDmatFull[i][2]);
                    }
                }
                if((t2Loc-sizeAB) < t3Loc){
                    for(int i=(t2Loc+1); i<=(t3Loc+sizeAB); ++i){
                        d5 += ABmatFull[i][0]*(ABmatFull[i-1][2]-ABmatFull[i][2]);
                    }
                }
                std::cout << '\n' << "((A,B),(C,D))" << " " << t1Loc << " " << t2Loc << " " << t3Loc << " ";
                for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                    g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                    g1 = thetaNew*(CDmatFull[t1Loc][2]*distMat[distIndex][2] + CDmatFull[t1Loc][2]*distMat[distIndex][3] - CDmatFull[t1Loc][2]*distMat[distIndex][4]);
                    g2 = thetaNew*(ABmatFull[t2Loc][2]*distMat[distIndex][0] + ABmatFull[t2Loc][2]*distMat[distIndex][1] - ABmatFull[t2Loc][2]*distMat[distIndex][4]);
                    g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][4];
                    long double intMult{};
                    if(t1Loc<sizeCD && t2Loc<sizeAB){
                        intMult = symInt1(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc<sizeCD && t2Loc>=sizeAB && (t2Loc-sizeAB)<t3Loc){
                        intMult = symInt3sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc<sizeCD && t2Loc>=sizeAB && (t2Loc-sizeAB)==t3Loc){
                        intMult = symInt3tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t2Loc<sizeAB && t1Loc>=sizeCD && (t1Loc-sizeCD)<t3Loc){
                        intMult = symInt2sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t2Loc<sizeAB && t1Loc>=sizeCD && (t1Loc-sizeCD)==t3Loc){
                        intMult = symInt2tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t2Loc>=sizeAB && t1Loc>=sizeCD && (t1Loc-sizeCD)<(t2Loc-sizeAB) && (t2Loc-sizeAB)<t3Loc){
                        intMult = symInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t2Loc>=sizeAB && t1Loc>=sizeCD && (t1Loc-sizeCD)>(t2Loc-sizeAB) && (t1Loc-sizeCD)<t3Loc){
                        intMult = symInt4sep(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                    }
                    if(t2Loc>=sizeAB && t1Loc>=sizeCD && (t1Loc-sizeCD)==(t2Loc-sizeAB) && (t2Loc-sizeAB)<t3Loc){
                        intMult = symInt4tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        intMult += symInt4tog12(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                    }
                    if(t2Loc>=sizeAB && t1Loc>=sizeCD && (t1Loc-sizeCD)<(t2Loc-sizeAB) && (t2Loc-sizeAB)==t3Loc){
                        intMult = symInt4tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t2Loc>=sizeAB && t1Loc>=sizeCD && (t1Loc-sizeCD)>(t2Loc-sizeAB) && (t1Loc-sizeCD)==t3Loc){
                        intMult = symInt4tog23(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                    }
                    if(t2Loc>=sizeAB && t1Loc>=sizeCD && (t1Loc-sizeCD)==(t2Loc-sizeAB) && (t2Loc-sizeAB)==t3Loc){
                        intMult = symInt4togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        intMult += symInt4togAll(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                    }
                    intVec.push_back(intMult);
                    std::cout << intMult << " ";
                    for(int pattern=0; pattern<15; ++pattern){
                        probs[pattern] += intMult*g0*coefMat[pattern][distIndex]/256.0;
                    }
                }
                intMat.push_back(intVec);
                intVec.clear();
            }
        }
    }
    
    
    // ((A,C),(B,D))

    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
         psi3a = rootMat[t3Loc][0];
         psi3b = rootMat[t3Loc][1];
        for(int t1Loc=0; t1Loc<=t3Loc; ++t1Loc){
            psi1a = rootMat[t1Loc][0];
            psi1b = rootMat[t1Loc][1];
            for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
                psi2a = rootMat[t2Loc][0];
                psi2b = rootMat[t2Loc][1];
                d2 = tau1*(gammaC-CDmat[0][2]);
                d4 = tau1*(gammaD-CDmat[0][2]);
                d1 = tau2*(gammaA-ABmat[0][2]);
                d3 = tau2*(gammaB-ABmat[0][2]);
                d2 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                d4 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                d1 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                d3 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                d5 = 0;
                if(t1Loc !=0){
                    for(int i=1; i<=t1Loc; ++i){
                    d3 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    d4 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(sizeAB>1){
                    for(int i=1; i<sizeAB; ++i){
                        d1 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                        d3 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                    }
                }
                if(sizeCD>1){
                    for(int i=1; i<sizeCD; ++i){
                        d2 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                        d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                    }
                }
                if(t2Loc !=0){
                    for(int i=1; i<=t2Loc; ++i){
                    d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t1Loc< t3Loc){
                    for(int i=(t1Loc+1); i<=t3Loc; ++i){
                        d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t2Loc < t3Loc){
                    for(int i=(t2Loc+1); i<=t3Loc; ++i){
                        d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                std::cout << '\n' << "((A,C),(B,D))"  << " " << t1Loc << " " << t2Loc << " " << t3Loc << " ";
                for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                    g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                    g1 = thetaNew*(rootMat[t1Loc][2]*distMat[distIndex][2] + rootMat[t1Loc][2]*distMat[distIndex][3] - rootMat[t1Loc][2]*distMat[distIndex][4]);
                    g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][0] + rootMat[t2Loc][2]*distMat[distIndex][1] - rootMat[t2Loc][2]*distMat[distIndex][4]);
                    g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][4];
                    long double intMult{};
                    if( t1Loc<t2Loc && t2Loc<t3Loc){
                        intMult = symInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if( t1Loc>t2Loc && t1Loc<t3Loc){
                        intMult = symInt4sep(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                    }
                    if( t1Loc==t2Loc && t1Loc<t3Loc){
                        intMult = symInt4tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        intMult += symInt4tog12(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                    }
                    if( t1Loc<t2Loc && t2Loc==t3Loc){
                        intMult = symInt4tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if( t1Loc>t2Loc && t1Loc==t3Loc){
                        intMult = symInt4tog23(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                    }
                    if( t1Loc==t2Loc && t1Loc==t3Loc){
                        intMult = symInt4togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        intMult += symInt4togAll(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                    }
                    intVec.push_back(intMult);
                    std::cout << intMult << " ";
                    for(int pattern=0; pattern<15; ++pattern){
                        int index1 {permMat[pattern][1]-1};
                        probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                    }
                }
                intMat.push_back(intVec);
                intVec.clear();
            }
        }
    }

    //((A,D),(B,C))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
         psi3a = rootMat[t3Loc][0];
         psi3b = rootMat[t3Loc][1];
        for(int t1Loc=0; t1Loc<=t3Loc; ++t1Loc){
            psi1a = rootMat[t1Loc][0];
            psi1b = rootMat[t1Loc][1];
            for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
                psi2a = rootMat[t2Loc][0];
                psi2b = rootMat[t2Loc][1];
                d4 = tau1*(gammaC-CDmat[0][2]);
                d2 = tau1*(gammaD-CDmat[0][2]);
                d1 = tau2*(gammaA-ABmat[0][2]);
                d3 = tau2*(gammaB-ABmat[0][2]);
                d2 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                d4 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                d1 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                d3 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                d5 = 0;
                if(t1Loc !=0){
                    for(int i=1; i<=t1Loc; ++i){
                    d3 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    d4 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(sizeAB>1){
                    for(int i=1; i<sizeAB; ++i){
                        d1 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                        d3 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                    }
                }
                if(sizeCD>1){
                    for(int i=1; i<sizeCD; ++i){
                        d2 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                        d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                    }
                }
                if(t2Loc !=0){
                    for(int i=1; i<=t2Loc; ++i){
                    d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t1Loc< t3Loc){
                    for(int i=(t1Loc+1); i<=t3Loc; ++i){
                        d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t2Loc < t3Loc){
                    for(int i=(t2Loc+1); i<=t3Loc; ++i){
                        d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                std::cout << '\n' << "((A,D),(B,D))"  << " " << t1Loc << " " << t2Loc << " " << t3Loc << " ";
                for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                    g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                    g1 = thetaNew*(rootMat[t1Loc][2]*distMat[distIndex][2] + rootMat[t1Loc][2]*distMat[distIndex][3] - rootMat[t1Loc][2]*distMat[distIndex][4]);
                    g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][0] + rootMat[t2Loc][2]*distMat[distIndex][1] - rootMat[t2Loc][2]*distMat[distIndex][4]);
                    g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][4];
                    long double intMult{};
                    if( t1Loc<t2Loc && t2Loc<t3Loc){
                        intMult = symInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if( t1Loc>t2Loc && t1Loc<t3Loc){
                        intMult = symInt4sep(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                    }
                    if( t1Loc==t2Loc && t1Loc<t3Loc){
                        intMult = symInt4tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        intMult += symInt4tog12(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                    }
                    if( t1Loc<t2Loc && t2Loc==t3Loc){
                        intMult = symInt4tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if( t1Loc>t2Loc && t1Loc==t3Loc){
                        intMult = symInt4tog23(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                    }
                    if(t1Loc==t2Loc && t1Loc==t3Loc){
                        intMult = symInt4togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        intMult += symInt4togAll(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                    }
                     intVec.push_back(intMult);
                    std::cout << intMult << " ";
                    for(int pattern=0; pattern<15; ++pattern){
                        int index1 {permMat[pattern][2]-1};
                        probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                    }
                }
                intMat.push_back(intVec);
                intVec.clear();
            }
        }
    }
    
    // (A,(B,(C,D)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
         psi3a = rootMat[t3Loc][0];
         psi3b = rootMat[t3Loc][1];
        for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
            psi2a = rootMat[t2Loc][0];
            psi2b = rootMat[t2Loc][1];
            for(int t1Loc=0; (t1Loc-sizeCD)<=t2Loc; ++t1Loc){
                psi1a = CDmatFull[t1Loc][0];
                psi1b = CDmatFull[t1Loc][1];
                d3 = tau1*(gammaC-CDmat[0][2]);
                d4 = tau1*(gammaD-CDmat[0][2]);
                d1 = tau2*(gammaA-ABmat[0][2]);
                d2 = tau2*(gammaB-ABmat[0][2]);
                d1 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                d2 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                d5 = 0;
                if(sizeAB>1){
                    for(int i=1; i<sizeAB; ++i){
                        d1 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                        d2 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                    }
                }
                if(t3Loc !=0){
                    for(int i=1; i<=t3Loc; ++i){
                    d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t2Loc !=0){
                    for(int i=1; i<=t2Loc; ++i){
                    d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t2Loc<t3Loc){
                    for(int i=(t2Loc+1); i<=t3Loc; ++i){
                        d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t1Loc !=0){
                    for(int i=1; i<=t1Loc; ++i){
                    d3 += CDmatFull[i][0]*(CDmatFull[i-1][2]-CDmatFull[i][2]);
                    d4 += CDmatFull[i][0]*(CDmatFull[i-1][2]-CDmatFull[i][2]);
                    }
                }
                if((t1Loc-sizeCD)<t2Loc){
                    for(int i=(t1Loc+1); i<=(t2Loc+sizeCD); ++i){
                        d5 += CDmatFull[i][0]*(CDmatFull[i-1][2]-CDmatFull[i][2]);
                    }
                }
                std::cout << '\n' << "(A,(B,(C,D)))" << " " << t1Loc << " " << t2Loc << " " << t3Loc << " ";
                for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                    g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                    g1 = thetaNew*(CDmatFull[t1Loc][2]*distMat[distIndex][2] + CDmatFull[t1Loc][2]*distMat[distIndex][3] - CDmatFull[t1Loc][2]*distMat[distIndex][4]);
                    g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                    g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                    long double intMult{};
                    if(t1Loc<sizeCD &&  t2Loc<t3Loc){
                        intMult = symInt3sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc<sizeCD &&  t2Loc==t3Loc){
                        intMult = symInt3tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc>=sizeCD && (t1Loc-sizeCD)<t2Loc && t2Loc<t3Loc){
                        intMult = symInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc>=sizeCD && (t1Loc-sizeCD)==t2Loc && t2Loc<t3Loc){
                        intMult = symInt4tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc>=sizeCD && (t1Loc-sizeCD)<t2Loc && t2Loc==t3Loc){
                        intMult = symInt4tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc>=sizeCD && (t1Loc-sizeCD)==t2Loc && t2Loc==t3Loc){
                        intMult = symInt4togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    intVec.push_back(intMult);
                    std::cout << intMult << " ";
                    for(int pattern=0; pattern<15; ++pattern){
                        probs[pattern] += intMult*g0*coefMat[pattern][distIndex]/256.0;
                    }
                }
                intMat.push_back(intVec);
                intVec.clear();
            }
        }
    }
    
    // (A,(C,(B,D)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
         psi3a = rootMat[t3Loc][0];
         psi3b = rootMat[t3Loc][1];
        for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
            psi2a = rootMat[t2Loc][0];
            psi2b = rootMat[t2Loc][1];
            for(int t1Loc=0; t1Loc<=t2Loc; ++t1Loc){
                psi1a = rootMat[t1Loc][0];
                psi1b = rootMat[t1Loc][1];
                d2 = tau1*(gammaC-CDmat[0][2]);
                d4 = tau1*(gammaD-CDmat[0][2]);
                d1 = tau2*(gammaA-ABmat[0][2]);
                d3 = tau2*(gammaB-ABmat[0][2]);
                d1 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                d3 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                d2 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                d4 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                d5 = 0;
                if(sizeAB>1){
                    for(int i=1; i<sizeAB; ++i){
                        d1 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                        d3 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                    }
                }
                if(sizeCD>1){
                    for(int i=1; i<sizeCD; ++i){
                        d2 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                        d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                    }
                }
                if(t3Loc !=0){
                    for(int i=1; i<=t3Loc; ++i){
                    d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t2Loc !=0){
                    for(int i=1; i<=t2Loc; ++i){
                        d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t2Loc<t3Loc){
                    for(int i=(t2Loc+1); i<=t3Loc; ++i){
                        d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t1Loc !=0){
                    for(int i=1; i<=t1Loc; ++i){
                    d3 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    d4 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t1Loc<t2Loc){
                    for(int i=(t1Loc+1); i<=t2Loc; ++i){
                        d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                std::cout << '\n' << "(A,(C,(B,D)))" << " " << t1Loc << " " << t2Loc << " " << t3Loc << " ";
                for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                    g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                    g1 = thetaNew*(rootMat[t1Loc][2]*distMat[distIndex][2] + rootMat[t1Loc][2]*distMat[distIndex][3] - rootMat[t1Loc][2]*distMat[distIndex][4]);
                    g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                    g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                    long double intMult{};
                    if(t1Loc<t2Loc && t2Loc<t3Loc){
                        intMult = symInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc==t2Loc && t2Loc<t3Loc){
                        intMult = symInt4tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc<t2Loc && t2Loc==t3Loc){
                        intMult = symInt4tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc==t2Loc && t2Loc==t3Loc){
                        intMult = symInt4togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    intVec.push_back(intMult);
                    std::cout << intMult << " ";
                    for(int pattern=0; pattern<15; ++pattern){
                        int index1 {permMat[pattern][1]-1};
                        probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                    }
                }
                intMat.push_back(intVec);
                intVec.clear();
            }
        }
    }
    
    // (A,(D,(B,C)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
         psi3a = rootMat[t3Loc][0];
         psi3b = rootMat[t3Loc][1];
        for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
            psi2a = rootMat[t2Loc][0];
            psi2b = rootMat[t2Loc][1];
            for(int t1Loc=0; t1Loc<=t2Loc; ++t1Loc){
                psi1a = rootMat[t1Loc][0];
                psi1b = rootMat[t1Loc][1];
                d4 = tau1*(gammaC-CDmat[0][2]);
                d2 = tau1*(gammaD-CDmat[0][2]);
                d1 = tau2*(gammaA-ABmat[0][2]);
                d3 = tau2*(gammaB-ABmat[0][2]);
                d1 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                d3 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                d2 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                d4 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                d5 = 0;
                if(sizeAB>1){
                    for(int i=1; i<sizeAB; ++i){
                        d1 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                        d3 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                    }
                }
                if(sizeCD>1){
                    for(int i=1; i<sizeCD; ++i){
                        d2 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                        d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                    }
                }
                if(t3Loc !=0){
                    for(int i=1; i<=t3Loc; ++i){
                    d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t2Loc !=0){
                    for(int i=1; i<=t2Loc; ++i){
                    d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t2Loc<t3Loc){
                    for(int i=(t2Loc+1); i<=t3Loc; ++i){
                        d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t1Loc !=0){
                    for(int i=1; i<=t1Loc; ++i){
                    d3 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    d4 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t1Loc<t2Loc){
                    for(int i=(t1Loc+1); i<=t2Loc; ++i){
                        d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                std::cout << '\n' << "(A,(D,(B,C)))" << " " << t1Loc << " " << t2Loc << " " << t3Loc << " ";
                for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                    g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                    g1 = thetaNew*(rootMat[t1Loc][2]*distMat[distIndex][2] + rootMat[t1Loc][2]*distMat[distIndex][3] - rootMat[t1Loc][2]*distMat[distIndex][4]);
                    g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                    g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                    long double intMult{};
                    if(t1Loc<t2Loc && t2Loc<t3Loc){
                        intMult = symInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc==t2Loc && t2Loc<t3Loc){
                        intMult = symInt4tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc<t2Loc && t2Loc==t3Loc){
                        intMult = symInt4tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc==t2Loc && t2Loc==t3Loc){
                        intMult = symInt4togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                     intVec.push_back(intMult);
                    std::cout << intMult << " ";
                    for(int pattern=0; pattern<15; ++pattern){
                        int index1 {permMat[pattern][2]-1};
                        probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                    }
                }
                intMat.push_back(intVec);
                intVec.clear();
            }
        }
    }
    
    //(B,(A,(C,D)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
          psi3a = rootMat[t3Loc][0];
          psi3b = rootMat[t3Loc][1];
         for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
             psi2a = rootMat[t2Loc][0];
             psi2b = rootMat[t2Loc][1];
             for(int t1Loc=0; (t1Loc-sizeCD)<=t2Loc; ++t1Loc){
                 psi1a = CDmatFull[t1Loc][0];
                 psi1b = CDmatFull[t1Loc][1];
                 d3 = tau1*(gammaC-CDmat[0][2]);
                 d4 = tau1*(gammaD-CDmat[0][2]);
                 d2 = tau2*(gammaA-ABmat[0][2]);
                 d1 = tau2*(gammaB-ABmat[0][2]);
                 d1 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                 d2 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                 d5 = 0;
                 if(sizeAB>1){
                     for(int i=1; i<sizeAB; ++i){
                         d1 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                         d2 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                     }
                 }
                 if(t3Loc !=0){
                     for(int i=1; i<=t3Loc; ++i){
                     d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t2Loc !=0){
                     for(int i=1; i<=t2Loc; ++i){
                     d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t2Loc<t3Loc){
                     for(int i=(t2Loc+1); i<=t3Loc; ++i){
                         d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t1Loc !=0){
                     for(int i=1; i<=t1Loc; ++i){
                     d3 += CDmatFull[i][0]*(CDmatFull[i-1][2]-CDmatFull[i][2]);
                     d4 += CDmatFull[i][0]*(CDmatFull[i-1][2]-CDmatFull[i][2]);
                     }
                 }
                 if((t1Loc-sizeCD)<t2Loc){
                     for(int i=(t1Loc+1); i<=(t2Loc+sizeCD); ++i){
                         d5 += CDmatFull[i][0]*(CDmatFull[i-1][2]-CDmatFull[i][2]);
                     }
                 }
                 std::cout << '\n' << "(B,(A,(C,D)))" << " " << t1Loc << " " << t2Loc << " " << t3Loc << " ";
                 for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                     g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                     g1 = thetaNew*(CDmatFull[t1Loc][2]*distMat[distIndex][2] + CDmatFull[t1Loc][2]*distMat[distIndex][3] - CDmatFull[t1Loc][2]*distMat[distIndex][4]);
                     g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                     g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                     long double intMult{};
                     if(t1Loc<sizeCD &&  t2Loc<t3Loc){
                         intMult = symInt3sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc<sizeCD &&  t2Loc==t3Loc){
                         intMult = symInt3tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeCD && (t1Loc-sizeCD)<t2Loc && t2Loc<t3Loc){
                         intMult = symInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeCD && (t1Loc-sizeCD)==t2Loc && t2Loc<t3Loc){
                         intMult = symInt4tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeCD && (t1Loc-sizeCD)<t2Loc && t2Loc==t3Loc){
                         intMult = symInt4tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeCD && (t1Loc-sizeCD)==t2Loc && t2Loc==t3Loc){
                         intMult = symInt4togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     intVec.push_back(intMult);
                     std::cout << intMult << " ";
                     for(int pattern=0; pattern<15; ++pattern){
                         int index1 {permMat[pattern][3]-1};
                         probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                     }
                 }
                 intMat.push_back(intVec);
                 intVec.clear();
             }
         }
     }
    
    //(B,(C,(A,D)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
         psi3a = rootMat[t3Loc][0];
         psi3b = rootMat[t3Loc][1];
        for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
            psi2a = rootMat[t2Loc][0];
            psi2b = rootMat[t2Loc][1];
            for(int t1Loc=0; t1Loc<=t2Loc; ++t1Loc){
                psi1a = rootMat[t1Loc][0];
                psi1b = rootMat[t1Loc][1];
                d2 = tau1*(gammaC-CDmat[0][2]);
                d4 = tau1*(gammaD-CDmat[0][2]);
                d3 = tau2*(gammaA-ABmat[0][2]);
                d1 = tau2*(gammaB-ABmat[0][2]);
                d1 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                d3 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                d2 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                d4 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                d5 = 0;
                if(sizeAB>1){
                    for(int i=1; i<sizeAB; ++i){
                        d1 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                        d3 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                    }
                }
                if(sizeCD>1){
                    for(int i=1; i<sizeCD; ++i){
                        d2 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                        d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                    }
                }
                if(t3Loc !=0){
                    for(int i=1; i<=t3Loc; ++i){
                    d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t2Loc !=0){
                    for(int i=1; i<=t2Loc; ++i){
                    d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t2Loc<t3Loc){
                    for(int i=(t2Loc+1); i<=t3Loc; ++i){
                        d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t1Loc !=0){
                    for(int i=1; i<=t1Loc; ++i){
                    d3 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    d4 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t1Loc<t2Loc){
                    for(int i=(t1Loc+1); i<=t2Loc; ++i){
                        d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                std::cout << '\n' << "(B,(C,(A,D)))" << " " << t1Loc << " " << t2Loc << " " << t3Loc << " ";
                for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                    g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                    g1 = thetaNew*(rootMat[t1Loc][2]*distMat[distIndex][2] + rootMat[t1Loc][2]*distMat[distIndex][3] - rootMat[t1Loc][2]*distMat[distIndex][4]);
                    g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                    g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                    long double intMult{};
                    if(t1Loc<t2Loc && t2Loc<t3Loc){
                        intMult = symInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc==t2Loc && t2Loc<t3Loc){
                        intMult = symInt4tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc<t2Loc && t2Loc==t3Loc){
                        intMult = symInt4tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc==t2Loc && t2Loc==t3Loc){
                        intMult = symInt4togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    intVec.push_back(intMult);
                    std::cout << intMult << " ";
                    for(int pattern=0; pattern<15; ++pattern){
                        int index1 {permMat[pattern][4]-1};
                        probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                    }
                }
                intMat.push_back(intVec);
                intVec.clear();
            }
        }
    }
    
    //(B,(D,(A,C)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
         psi3a = rootMat[t3Loc][0];
         psi3b = rootMat[t3Loc][1];
        for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
            psi2a = rootMat[t2Loc][0];
            psi2b = rootMat[t2Loc][1];
            for(int t1Loc=0; t1Loc<=t2Loc; ++t1Loc){
                psi1a = rootMat[t1Loc][0];
                psi1b = rootMat[t1Loc][1];
                d4 = tau1*(gammaC-CDmat[0][2]);
                d2 = tau1*(gammaD-CDmat[0][2]);
                d3 = tau2*(gammaA-ABmat[0][2]);
                d1 = tau2*(gammaB-ABmat[0][2]);
                d1 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                d3 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                d2 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                d4 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                d5 = 0;
                if(sizeAB>1){
                    for(int i=1; i<sizeAB; ++i){
                        d1 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                        d3 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                    }
                }
                if(sizeCD>1){
                    for(int i=1; i<sizeCD; ++i){
                        d2 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                        d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                    }
                }
                if(t3Loc !=0){
                    for(int i=1; i<=t3Loc; ++i){
                    d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t2Loc !=0){
                    for(int i=1; i<=t2Loc; ++i){
                    d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t2Loc<t3Loc){
                    for(int i=(t2Loc+1); i<=t3Loc; ++i){
                        d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t1Loc !=0){
                    for(int i=1; i<=t1Loc; ++i){
                    d3 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    d4 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t1Loc<t2Loc){
                    for(int i=(t1Loc+1); i<=t2Loc; ++i){
                        d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                std::cout << '\n' << "(B,(D,(A,C)))" << " " << t1Loc << " " << t2Loc << " " << t3Loc << " ";
                for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                    g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                    g1 = thetaNew*(rootMat[t1Loc][2]*distMat[distIndex][2] + rootMat[t1Loc][2]*distMat[distIndex][3] - rootMat[t1Loc][2]*distMat[distIndex][4]);
                    g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                    g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                    long double intMult{};
                    if(t1Loc<t2Loc && t2Loc<t3Loc){
                        intMult = symInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc==t2Loc && t2Loc<t3Loc){
                        intMult = symInt4tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc<t2Loc && t2Loc==t3Loc){
                        intMult = symInt4tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc==t2Loc && t2Loc==t3Loc){
                        intMult = symInt4togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                     intVec.push_back(intMult);
                    std::cout << intMult << " ";
                    for(int pattern=0; pattern<15; ++pattern){
                        int index1 {permMat[pattern][5]-1};
                        probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                    }
                }
                intMat.push_back(intVec);
                intVec.clear();
            }
        }
    }
    
    //(C,(A,(B,D)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
         psi3a = rootMat[t3Loc][0];
         psi3b = rootMat[t3Loc][1];
        for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
            psi2a = rootMat[t2Loc][0];
            psi2b = rootMat[t2Loc][1];
            for(int t1Loc=0; t1Loc<=t2Loc; ++t1Loc){
                psi1a = rootMat[t1Loc][0];
                psi1b = rootMat[t1Loc][1];
                d1 = tau1*(gammaC-CDmat[0][2]);
                d4 = tau1*(gammaD-CDmat[0][2]);
                d2 = tau2*(gammaA-ABmat[0][2]);
                d3 = tau2*(gammaB-ABmat[0][2]);
                d2 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                d3 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                d1 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                d4 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                d5 = 0;
                if(sizeAB>1){
                    for(int i=1; i<sizeAB; ++i){
                        d2 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                        d3 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                    }
                }
                if(sizeCD>1){
                    for(int i=1; i<sizeCD; ++i){
                        d1 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                        d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                    }
                }
                if(t3Loc !=0){
                    for(int i=1; i<=t3Loc; ++i){
                    d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t2Loc !=0){
                    for(int i=1; i<=t2Loc; ++i){
                    d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t2Loc<t3Loc){
                    for(int i=(t2Loc+1); i<=t3Loc; ++i){
                        d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t1Loc !=0){
                    for(int i=1; i<=t1Loc; ++i){
                    d3 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    d4 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t1Loc<t2Loc){
                    for(int i=(t1Loc+1); i<=t2Loc; ++i){
                        d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                std::cout << '\n' << "(C,(A,(B,D)))" << " " << t1Loc << " " << t2Loc << " " << t3Loc << " ";
                for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                    g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                    g1 = thetaNew*(rootMat[t1Loc][2]*distMat[distIndex][2] + rootMat[t1Loc][2]*distMat[distIndex][3] - rootMat[t1Loc][2]*distMat[distIndex][4]);
                    g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                    g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                    long double intMult{};
                    if(t1Loc<t2Loc && t2Loc<t3Loc){
                        intMult = symInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc==t2Loc && t2Loc<t3Loc){
                        intMult = symInt4tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc<t2Loc && t2Loc==t3Loc){
                        intMult = symInt4tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc==t2Loc && t2Loc==t3Loc){
                        intMult = symInt4togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                     intVec.push_back(intMult);
                    std::cout << intMult << " ";
                    for(int pattern=0; pattern<15; ++pattern){
                        int index1 {permMat[pattern][6]-1};
                        probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                    }
                }
                intMat.push_back(intVec);
                intVec.clear();
            }
        }
    }
    
    //(C,(B,(A,D)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
         psi3a = rootMat[t3Loc][0];
         psi3b = rootMat[t3Loc][1];
        for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
            psi2a = rootMat[t2Loc][0];
            psi2b = rootMat[t2Loc][1];
            for(int t1Loc=0; t1Loc<=t2Loc; ++t1Loc){
                psi1a = rootMat[t1Loc][0];
                psi1b = rootMat[t1Loc][1];
                d1 = tau1*(gammaC-CDmat[0][2]);
                d4 = tau1*(gammaD-CDmat[0][2]);
                d3 = tau2*(gammaA-ABmat[0][2]);
                d2 = tau2*(gammaB-ABmat[0][2]);
                d2 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                d3 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                d1 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                d4 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                d5 = 0;
                if(sizeAB>1){
                    for(int i=1; i<sizeAB; ++i){
                        d2 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                        d3 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                    }
                }
                if(sizeCD>1){
                    for(int i=1; i<sizeCD; ++i){
                        d1 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                        d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                    }
                }
                if(t3Loc !=0){
                    for(int i=1; i<=t3Loc; ++i){
                    d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t2Loc !=0){
                    for(int i=1; i<=t2Loc; ++i){
                    d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t2Loc<t3Loc){
                    for(int i=(t2Loc+1); i<=t3Loc; ++i){
                        d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t1Loc !=0){
                    for(int i=1; i<=t1Loc; ++i){
                    d3 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    d4 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                if(t1Loc<t2Loc){
                    for(int i=(t1Loc+1); i<=t2Loc; ++i){
                        d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                    }
                }
                std::cout << '\n' << "(C,(B,(A,D)))" << " " << t1Loc << " " << t2Loc << " " << t3Loc << " ";
                for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                    g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                    g1 = thetaNew*(rootMat[t1Loc][2]*distMat[distIndex][2] + rootMat[t1Loc][2]*distMat[distIndex][3] - rootMat[t1Loc][2]*distMat[distIndex][4]);
                    g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                    g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                    long double intMult{};
                    if(t1Loc<t2Loc && t2Loc<t3Loc){
                        intMult = symInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc==t2Loc && t2Loc<t3Loc){
                        intMult = symInt4tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc<t2Loc && t2Loc==t3Loc){
                        intMult = symInt4tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                    if(t1Loc==t2Loc && t2Loc==t3Loc){
                        intMult = symInt4togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                    }
                     intVec.push_back(intMult);
                    std::cout << intMult << " ";
                    for(int pattern=0; pattern<15; ++pattern){
                        int index1 {permMat[pattern][7]-1};
                        probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                    }
                }
                intMat.push_back(intVec);
                intVec.clear();
            }
        }
    }
    
    
    //(C,(D,(A,B)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
          psi3a = rootMat[t3Loc][0];
          psi3b = rootMat[t3Loc][1];
         for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
             psi2a = rootMat[t2Loc][0];
             psi2b = rootMat[t2Loc][1];
             for(int t1Loc=0; (t1Loc-sizeAB)<=t2Loc; ++t1Loc){
                 psi1a = ABmatFull[t1Loc][0];
                 psi1b = ABmatFull[t1Loc][1];
                 d1 = tau1*(gammaC-CDmat[0][2]);
                 d2 = tau1*(gammaD-CDmat[0][2]);
                 d3 = tau2*(gammaA-ABmat[0][2]);
                 d4 = tau2*(gammaB-ABmat[0][2]);
                 d1 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                 d2 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                 d5 = 0;
                 if(sizeCD>1){
                     for(int i=1; i<sizeCD; ++i){
                         d1 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                         d2 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                     }
                 }
                 if(t3Loc !=0){
                     for(int i=1; i<=t3Loc; ++i){
                     d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t2Loc !=0){
                     for(int i=1; i<=t2Loc; ++i){
                     d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t2Loc<t3Loc){
                     for(int i=(t2Loc+1); i<=t3Loc; ++i){
                         d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t1Loc !=0){
                     for(int i=1; i<=t1Loc; ++i){
                     d3 += ABmatFull[i][0]*(ABmatFull[i-1][2]-ABmatFull[i][2]);
                     d4 += ABmatFull[i][0]*(ABmatFull[i-1][2]-ABmatFull[i][2]);
                     }
                 }
                 if((t1Loc-sizeAB)<t2Loc){
                     for(int i=(t1Loc+1); i<=(t2Loc+sizeAB); ++i){
                         d5 += ABmatFull[i][0]*(ABmatFull[i-1][2]-ABmatFull[i][2]);
                     }
                 }
                 std::cout << '\n' << "(C,(D,(A,B)))" << " " << t1Loc << " " << t2Loc << " " << t3Loc << " ";
                 for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                     g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                     g1 = thetaNew*(ABmatFull[t1Loc][2]*distMat[distIndex][2] + ABmatFull[t1Loc][2]*distMat[distIndex][3] - ABmatFull[t1Loc][2]*distMat[distIndex][4]);
                     g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                     g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                     long double intMult{};
                     if(t1Loc<sizeAB &&  t2Loc<t3Loc){
                         intMult = symInt3sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc<sizeAB &&  t2Loc==t3Loc){
                         intMult = symInt3tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeAB && (t1Loc-sizeAB)<t2Loc && t2Loc<t3Loc){
                         intMult = symInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeAB && (t1Loc-sizeAB)==t2Loc && t2Loc<t3Loc){
                         intMult = symInt4tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeAB && (t1Loc-sizeAB)<t2Loc && t2Loc==t3Loc){
                         intMult = symInt4tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeAB && (t1Loc-sizeAB)==t2Loc && t2Loc==t3Loc){
                         intMult = symInt4togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                      intVec.push_back(intMult);
                     std::cout << intMult << " ";
                     for(int pattern=0; pattern<15; ++pattern){
                         int index1 {permMat[pattern][8]-1};
                         probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                     }
                 }
                 intMat.push_back(intVec);
                 intVec.clear();
             }
         }
     }
    
    // (D,(A,(B,C)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
            psi3a = rootMat[t3Loc][0];
            psi3b = rootMat[t3Loc][1];
           for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
               psi2a = rootMat[t2Loc][0];
               psi2b = rootMat[t2Loc][1];
               for(int t1Loc=0; t1Loc<=t2Loc; ++t1Loc){
                   psi1a = rootMat[t1Loc][0];
                   psi1b = rootMat[t1Loc][1];
                   d4 = tau1*(gammaC-CDmat[0][2]);
                   d1 = tau1*(gammaD-CDmat[0][2]);
                   d2 = tau2*(gammaA-ABmat[0][2]);
                   d3 = tau2*(gammaB-ABmat[0][2]);
                   d2 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                   d3 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                   d1 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                   d4 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                   d5 = 0;
                   if(sizeAB>1){
                       for(int i=1; i<sizeAB; ++i){
                           d2 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                           d3 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                       }
                   }
                   if(sizeCD>1){
                       for(int i=1; i<sizeCD; ++i){
                           d1 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                           d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                       }
                   }
                   if(t3Loc !=0){
                       for(int i=1; i<=t3Loc; ++i){
                       d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                       }
                   }
                   if(t2Loc !=0){
                       for(int i=1; i<=t2Loc; ++i){
                       d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                       }
                   }
                   if(t2Loc<t3Loc){
                       for(int i=(t2Loc+1); i<=t3Loc; ++i){
                           d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                       }
                   }
                   if(t1Loc !=0){
                       for(int i=1; i<=t1Loc; ++i){
                       d3 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                       d4 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                       }
                   }
                   if(t1Loc<t2Loc){
                       for(int i=(t1Loc+1); i<=t2Loc; ++i){
                           d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                       }
                   }
                   std::cout << '\n' << "(D,(A,(B,C)))" << " " << t1Loc << " " << t2Loc << " " << t3Loc << " ";
                   for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                       g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                       g1 = thetaNew*(rootMat[t1Loc][2]*distMat[distIndex][2] + rootMat[t1Loc][2]*distMat[distIndex][3] - rootMat[t1Loc][2]*distMat[distIndex][4]);
                       g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                       g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                       long double intMult{};
                       if(t1Loc<t2Loc && t2Loc<t3Loc){
                           intMult = symInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                       }
                       if(t1Loc==t2Loc && t2Loc<t3Loc){
                           intMult = symInt4tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                       }
                       if(t1Loc<t2Loc && t2Loc==t3Loc){
                           intMult = symInt4tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                       }
                       if(t1Loc==t2Loc && t2Loc==t3Loc){
                           intMult = symInt4togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                       }
                        intVec.push_back(intMult);
                       std::cout << intMult << " ";
                       for(int pattern=0; pattern<15; ++pattern){
                           int index1 {permMat[pattern][9]-1};
                           probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                       }
                   }
                   intMat.push_back(intVec);
                   intVec.clear();
               }
           }
       }
    
    // (D,(B,(A,C)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
            psi3a = rootMat[t3Loc][0];
            psi3b = rootMat[t3Loc][1];
           for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
               psi2a = rootMat[t2Loc][0];
               psi2b = rootMat[t2Loc][1];
               for(int t1Loc=0; t1Loc<=t2Loc; ++t1Loc){
                   psi1a = rootMat[t1Loc][0];
                   psi1b = rootMat[t1Loc][1];
                   d4 = tau1*(gammaC-CDmat[0][2]);
                   d1 = tau1*(gammaD-CDmat[0][2]);
                   d3 = tau2*(gammaA-ABmat[0][2]);
                   d2 = tau2*(gammaB-ABmat[0][2]);
                   d2 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                   d3 += tau3*(ABmat[sizeAB-1][2]-rootMat[0][2]);
                   d1 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                   d4 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                   d5 = 0;
                   if(sizeAB>1){
                       for(int i=1; i<sizeAB; ++i){
                           d2 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                           d3 += ABmat[i][0]*(ABmat[i-1][2]-ABmat[i][2]);
                       }
                   }
                   if(sizeCD>1){
                       for(int i=1; i<sizeCD; ++i){
                           d1 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                           d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                       }
                   }
                   if(t3Loc !=0){
                       for(int i=1; i<=t3Loc; ++i){
                       d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                       }
                   }
                   if(t2Loc !=0){
                       for(int i=1; i<=t2Loc; ++i){
                       d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                       }
                   }
                   if(t2Loc<t3Loc){
                       for(int i=(t2Loc+1); i<=t3Loc; ++i){
                           d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                       }
                   }
                   if(t1Loc !=0){
                       for(int i=1; i<=t1Loc; ++i){
                       d3 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                       d4 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                       }
                   }
                   if(t1Loc<t2Loc){
                       for(int i=(t1Loc+1); i<=t2Loc; ++i){
                           d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                       }
                   }
                   std::cout << '\n' << "(D,(B,(A,C)))" << " " << t1Loc << " " << t2Loc << " " << t3Loc << " ";
                   for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                       g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                       g1 = thetaNew*(rootMat[t1Loc][2]*distMat[distIndex][2] + rootMat[t1Loc][2]*distMat[distIndex][3] - rootMat[t1Loc][2]*distMat[distIndex][4]);
                       g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                       g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                       long double intMult{};
                       if(t1Loc<t2Loc && t2Loc<t3Loc){
                           intMult = symInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                       }
                       if(t1Loc==t2Loc && t2Loc<t3Loc){
                           intMult = symInt4tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                       }
                       if(t1Loc<t2Loc && t2Loc==t3Loc){
                           intMult = symInt4tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                       }
                       if(t1Loc==t2Loc && t2Loc==t3Loc){
                           intMult = symInt4togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                       }
                       intVec.push_back(intMult);
                       std::cout << intMult << " ";
                       for(int pattern=0; pattern<15; ++pattern){
                           int index1 {permMat[pattern][10]-1};
                           probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                       }
                   }
                   intMat.push_back(intVec);
                   intVec.clear();
               }
           }
       }
    
    // (D,(C,(A,B)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
          psi3a = rootMat[t3Loc][0];
          psi3b = rootMat[t3Loc][1];
         for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
             psi2a = rootMat[t2Loc][0];
             psi2b = rootMat[t2Loc][1];
             for(int t1Loc=0; (t1Loc-sizeAB)<=t2Loc; ++t1Loc){
                 psi1a = ABmatFull[t1Loc][0];
                 psi1b = ABmatFull[t1Loc][1];
                 d2 = tau1*(gammaC-CDmat[0][2]);
                 d1 = tau1*(gammaD-CDmat[0][2]);
                 d3 = tau2*(gammaA-ABmat[0][2]);
                 d4 = tau2*(gammaB-ABmat[0][2]);
                 d1 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                 d2 += tau3*(CDmat[sizeCD-1][2]-rootMat[0][2]);
                 d5 = 0;
                 if(sizeCD>1){
                     for(int i=1; i<sizeCD; ++i){
                         d1 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                         d2 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                     }
                 }
                 if(t3Loc !=0){
                     for(int i=1; i<=t3Loc; ++i){
                     d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t2Loc !=0){
                     for(int i=1; i<=t2Loc; ++i){
                     d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t2Loc<t3Loc){
                     for(int i=(t2Loc+1); i<=t3Loc; ++i){
                         d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t1Loc !=0){
                     for(int i=1; i<=t1Loc; ++i){
                     d3 += ABmatFull[i][0]*(ABmatFull[i-1][2]-ABmatFull[i][2]);
                     d4 += ABmatFull[i][0]*(ABmatFull[i-1][2]-ABmatFull[i][2]);
                     }
                 }
                 if((t1Loc-sizeAB)<t2Loc){
                     for(int i=(t1Loc+1); i<=(t2Loc+sizeAB); ++i){
                         d5 += ABmatFull[i][0]*(ABmatFull[i-1][2]-ABmatFull[i][2]);
                     }
                 }
                 std::cout << '\n' << "(D,(A,(A,B)))" << " " << t1Loc << " " << t2Loc << " " << t3Loc << " ";
                 for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                     g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                     g1 = thetaNew*(ABmatFull[t1Loc][2]*distMat[distIndex][2] + ABmatFull[t1Loc][2]*distMat[distIndex][3] - ABmatFull[t1Loc][2]*distMat[distIndex][4]);
                     g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                     g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                     long double intMult{};
                     if(t1Loc<sizeAB &&  t2Loc<t3Loc){
                         intMult = symInt3sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc<sizeAB &&  t2Loc==t3Loc){
                        intMult = symInt3tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeAB && (t1Loc-sizeAB)<t2Loc && t2Loc<t3Loc){
                         intMult = symInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeAB && (t1Loc-sizeAB)==t2Loc && t2Loc<t3Loc){
                         intMult = symInt4tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeAB && (t1Loc-sizeAB)<t2Loc && t2Loc==t3Loc){
                         intMult = symInt4tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeAB && (t1Loc-sizeAB)==t2Loc && t2Loc==t3Loc){
                         intMult = symInt4togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     intVec.push_back(intMult);
                     std::cout << intMult << " ";
                     for(int pattern=0; pattern<15; ++pattern){
                         int index1 {permMat[pattern][11]-1};
                         probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                     }
                 }
                 intMat.push_back(intVec);
                 intVec.clear();
             }
         }
     }
  /*  for(int i=0; i<intMat.size(); ++i){
        for(int j=0; j<intMat[i].size(); ++j){
            std::cout << intMat[i][j] << " ";
        }
        std::cout << '\n';
    }
    
    std::vector<long double> sumVec{};
    for(int i=0; i<intMat.size(); ++i){
        long double testsum{
            intMat[i][0]-intMat[i][1]-intMat[i][2]-intMat[i][3]-intMat[i][4]-intMat[i][5]-intMat[i][6]
                +intMat[i][7]+2.0*intMat[i][8]+2.0*intMat[i][9]+2.0*intMat[i][10]+2.0*intMat[i][11]-4.0*intMat[i][12]
        };
        sumVec.push_back(testsum);
    }
    for(int i=0; i<intMat.size(); ++i){
        std::cout << sumVec[i] << " ";
    }
    std::cout << '\n';
    long double totalSum{};
    for(int i=0; i<intMat.size(); ++i){
        totalSum += sumVec[i];
     }
    std::cout << totalSum << '\n';
    
     std::vector<long double> sumVec2{};
    for(int i=0; i<intMat[0].size(); ++i){
        long double testSum{};
        for(int j=0; j<intMat.size(); ++j){
            testSum += intMat[j][i];
        }
        sumVec2.push_back(testSum);
    }
    for(int i=0; i<sumVec2.size(); ++i){
        std::cout << sumVec2[i] << " ";
    }
    std::cout << '\n';
    
    long double totalSum2{
        sumVec2[0]-sumVec2[1]-sumVec2[2]-sumVec2[3]-sumVec2[4]-sumVec2[5]-sumVec2[6]
            +sumVec2[7]+2.0*sumVec2[8]+2.0*sumVec2[9]+2.0*sumVec2[10]+2.0*sumVec2[11]-4.0*sumVec2[12]
    };
    std::cout << totalSum2 << '\n'; */
    
    
    return probs;
}


std::vector<long double> probAsymm(double tau1, double tau2, double tau3, double gammaA, double gammaB, double gammaC, double gammaD, doubleMat2 CDmat, doubleMat2 BCDmat, doubleMat2 rootMat, double theta){
    std::vector<long double> probs{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double thetaNew {4.0*theta/6.0};
    intMat2 coefMat{
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
    intMat2 permMat{
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
    intMat2 distMat{
        {0,0,0,0,0},
        {1,1,0,0,0},
        {0,0,1,1,0},
        {1,0,1,0,1},
        {1,0,0,1,1},
        {0,1,1,0,1},
        {0,1,0,1,1},
        {1,1,1,1,0},
        {1,1,1,0,1},
        {1,1,0,1,1},
        {1,0,1,1,1},
        {0,1,1,1,1},
        {1,1,1,1,1}
    };
  
    double g0{};
    double g1{};
    double g2{};
    double g3{};
    double d1{};
    double d2{};
    double d3{};
    double d4{};
    double d5{};
    double psi1a{};
    double psi1b{};
    double psi2a{};
    double psi2b{};
    double psi3a{};
    double psi3b{};
    doubleMat2 CDmatFull{CDmat};
    doubleMat2 BCDmatFull{BCDmat};
    for(int i=0; i<rootMat.size(); ++i){
        BCDmatFull.push_back(rootMat[i]);
    }
    for(int i=0; i<BCDmatFull.size(); ++i){
        CDmatFull.push_back(BCDmatFull[i]);
    }
    int sizeCD {static_cast<int>(CDmat.size())};
    int sizeBCD {static_cast<int>(BCDmat.size())};
    
    // ((A,B),(C,D))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
          psi3a = rootMat[t3Loc][0];
          psi3b = rootMat[t3Loc][1];
         for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
             psi2a = rootMat[t2Loc][0];
             psi2b = rootMat[t2Loc][1];
             for(int t1Loc=0; (t1Loc-sizeCD-sizeBCD)<=t3Loc; ++t1Loc){
                 psi1a = CDmatFull[t1Loc][0];
                 psi1b = CDmatFull[t1Loc][1];
                 d3 = tau1*(gammaC-CDmat[0][2]);
                 d4 = tau1*(gammaD-CDmat[0][2]);
                 d1 = tau3*(gammaA-rootMat[0][2]);
                 d2 = tau2*(gammaB-BCDmat[0][2]);
                 d2 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                 d5 = 0;
                 if(sizeBCD>1){
                     for(int i=1; i<sizeBCD; ++i){
                         d2 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                     }
                 }
                 if(t1Loc !=0){
                     for(int i=1; i<=t1Loc; ++i){
                         d3 += CDmatFull[i][0]*(CDmatFull[i-1][2]-CDmatFull[i][2]);
                         d4 += CDmatFull[i][0]*(CDmatFull[i-1][2]-CDmatFull[i][2]);
                     }
                 }
                 if(t2Loc !=0){
                     for(int i=1; i<=t2Loc; ++i){
                     d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if((t1Loc-sizeCD-sizeBCD) < t3Loc){
                     for(int i=(t1Loc+1); i<=(t3Loc+sizeCD+sizeBCD); ++i){
                         d5 += CDmatFull[i][0]*(CDmatFull[i-1][2]-CDmatFull[i][2]);
                     }
                 }
                 if(t2Loc < t3Loc){
                     for(int i=(t2Loc+1); i<=t3Loc; ++i){
                         d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 std::cout << t1Loc << " " << t2Loc << " " << t3Loc << '\n';
                 for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                     g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                     g1 = thetaNew*(CDmatFull[t1Loc][2]*distMat[distIndex][2] + CDmatFull[t1Loc][2]*distMat[distIndex][3] - CDmatFull[t1Loc][2]*distMat[distIndex][4]);
                     g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][0] + rootMat[t2Loc][2]*distMat[distIndex][1] - rootMat[t2Loc][2]*distMat[distIndex][4]);
                     g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][4];
                     long double intMult{};
                     if(t1Loc<sizeCD && t2Loc<t3Loc){
                         intMult = asymmInt2sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc<sizeCD && t2Loc==t3Loc){
                         intMult = asymmInt2tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeCD && t1Loc<(sizeCD+sizeBCD) && t2Loc<t3Loc){
                         intMult = asymmInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeCD && t1Loc<(sizeCD+sizeBCD) && t2Loc==t3Loc){
                         intMult = asymmInt4tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=(sizeCD+sizeBCD) && (t1Loc-sizeCD-sizeBCD)<t2Loc && t2Loc<t3Loc){
                         intMult = asymmInt5sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=(sizeCD+sizeBCD) && (t1Loc-sizeCD-sizeBCD)>t2Loc && (t1Loc-sizeCD-sizeBCD)<t3Loc){
                         intMult = asymmInt5sep(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                     }
                     if(t1Loc>=(sizeCD+sizeBCD) && (t1Loc-sizeCD-sizeBCD)==t2Loc && t2Loc<t3Loc){
                         intMult = asymmInt5tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                         intMult += asymmInt5tog12(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                     }
                     if(t1Loc>=(sizeCD+sizeBCD) && (t1Loc-sizeCD-sizeBCD)<t2Loc && t2Loc==t3Loc){
                         intMult = asymmInt5tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=(sizeCD+sizeBCD) && (t1Loc-sizeCD-sizeBCD)>t2Loc && (t1Loc-sizeCD-sizeBCD)==t3Loc){
                         intMult = asymmInt5tog23(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                     }
                     if(t1Loc>=(sizeCD+sizeBCD) && (t1Loc-sizeCD-sizeBCD)==t2Loc && t2Loc==t3Loc){
                         intMult = asymmInt5togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                         intMult += asymmInt5togAll(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                     }
                  //   intVec.push_back(intMult);
                     for(int pattern=0; pattern<15; ++pattern){
                         probs[pattern] += intMult*g0*coefMat[pattern][distIndex]/256.0;
                     }
                 }
              //   intMat.push_back(intVec);
               //  intVec.clear();
             }
         }
     }
     
    // ((A,C),(B,D))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
          psi3a = rootMat[t3Loc][0];
          psi3b = rootMat[t3Loc][1];
         for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
             psi2a = rootMat[t2Loc][0];
             psi2b = rootMat[t2Loc][1];
             for(int t1Loc=0; (t1Loc-sizeBCD)<=t3Loc; ++t1Loc){
                 psi1a = BCDmatFull[t1Loc][0];
                 psi1b = BCDmatFull[t1Loc][1];
                 d2 = tau1*(gammaC-CDmat[0][2]);
                 d4 = tau1*(gammaD-CDmat[0][2]);
                 d1 = tau3*(gammaA-rootMat[0][2]);
                 d3 = tau2*(gammaB-BCDmat[0][2]);
                 d5 = 0;
                 d2 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                 d4 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                 if(sizeCD>1){
                     for(int i=1; i<sizeCD; ++i){
                         d2 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                         d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                     }
                 }
                 
                 d2 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                 if(sizeBCD>1){
                     for(int i=1; i<sizeBCD; ++i){
                         d2 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                     }
                 }
                 if(t1Loc !=0){
                     for(int i=1; i<=t1Loc; ++i){
                         d3 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                         d4 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     }
                 }
                 if(t2Loc !=0){
                     for(int i=1; i<=t2Loc; ++i){
                     d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if((t1Loc-sizeBCD) < t3Loc){
                     for(int i=(t1Loc+1); i<=(t3Loc+sizeBCD); ++i){
                         d5 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     }
                 }
                 if(t2Loc < t3Loc){
                     for(int i=(t2Loc+1); i<=t3Loc; ++i){
                         d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                     g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                     g1 = thetaNew*(BCDmatFull[t1Loc][2]*distMat[distIndex][2] + BCDmatFull[t1Loc][2]*distMat[distIndex][3] - BCDmatFull[t1Loc][2]*distMat[distIndex][4]);
                     g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][0] + rootMat[t2Loc][2]*distMat[distIndex][1] - rootMat[t2Loc][2]*distMat[distIndex][4]);
                     g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][4];
                     long double intMult{};
                     if(t1Loc<sizeBCD  && t2Loc<t3Loc){
                         intMult = asymmInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc<sizeBCD  && t2Loc==t3Loc){
                         intMult = asymmInt4tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && (t1Loc-sizeBCD)<t2Loc && t2Loc<t3Loc){
                         intMult = asymmInt5sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && (t1Loc-sizeBCD)>t2Loc && (t1Loc-sizeBCD)<t3Loc){
                         intMult = asymmInt5sep(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                     }
                     if(t1Loc>=sizeBCD && (t1Loc-sizeBCD)==t2Loc && t2Loc<t3Loc){
                         intMult = asymmInt5tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                         intMult += asymmInt5tog12(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                     }
                     if(t1Loc>=sizeBCD && (t1Loc-sizeBCD)<t2Loc && t2Loc==t3Loc){
                         intMult = asymmInt5tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && (t1Loc-sizeBCD)>t2Loc && (t1Loc-sizeBCD)==t3Loc){
                         intMult = asymmInt5tog23(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                     }
                     if(t1Loc>=sizeBCD && (t1Loc-sizeBCD)==t2Loc && t2Loc==t3Loc){
                         intMult = asymmInt5togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                         intMult += asymmInt5togAll(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                     }
                  //   intVec.push_back(intMult);
                     for(int pattern=0; pattern<15; ++pattern){
                        int index1 {permMat[pattern][1]-1};
                        probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                     }
                 }
              //   intMat.push_back(intVec);
               //  intVec.clear();
             }
         }
     }
    
    // ((A,D),(B,C))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
          psi3a = rootMat[t3Loc][0];
          psi3b = rootMat[t3Loc][1];
         for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
             psi2a = rootMat[t2Loc][0];
             psi2b = rootMat[t2Loc][1];
             for(int t1Loc=0; (t1Loc-sizeBCD)<=t3Loc; ++t1Loc){
                 psi1a = BCDmatFull[t1Loc][0];
                 psi1b = BCDmatFull[t1Loc][1];
                 d4 = tau1*(gammaC-CDmat[0][2]);
                 d2 = tau1*(gammaD-CDmat[0][2]);
                 d1 = tau3*(gammaA-rootMat[0][2]);
                 d3 = tau2*(gammaB-BCDmat[0][2]);
                 d5 = 0;
                 d2 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                 d4 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                 if(sizeCD>1){
                     for(int i=1; i<sizeCD; ++i){
                         d2 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                         d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                     }
                 }
                 d2 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                 if(sizeBCD>1){
                     for(int i=1; i<sizeBCD; ++i){
                         d2 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                     }
                 }
                 if(t1Loc !=0){
                     for(int i=1; i<=t1Loc; ++i){
                         d3 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                         d4 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     }
                 }
                 if(t2Loc !=0){
                     for(int i=1; i<=t2Loc; ++i){
                     d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if((t1Loc-sizeBCD) < t3Loc){
                     for(int i=(t1Loc+1); i<=(t3Loc+sizeBCD); ++i){
                         d5 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     }
                 }
                 if(t2Loc < t3Loc){
                     for(int i=(t2Loc+1); i<=t3Loc; ++i){
                         d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                     g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                     g1 = thetaNew*(BCDmatFull[t1Loc][2]*distMat[distIndex][2] + BCDmatFull[t1Loc][2]*distMat[distIndex][3] - BCDmatFull[t1Loc][2]*distMat[distIndex][4]);
                     g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][0] + rootMat[t2Loc][2]*distMat[distIndex][1] - rootMat[t2Loc][2]*distMat[distIndex][4]);
                     g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][4];
                     long double intMult{};
                     if(t1Loc<sizeBCD  && t2Loc<t3Loc){
                         intMult = asymmInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc<sizeBCD  && t2Loc==t3Loc){
                         intMult = asymmInt4tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && (t1Loc-sizeBCD)<t2Loc && t2Loc<t3Loc){
                         intMult = asymmInt5sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && (t1Loc-sizeBCD)>t2Loc && (t1Loc-sizeBCD)<t3Loc){
                         intMult = asymmInt5sep(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                     }
                     if(t1Loc>=sizeBCD && (t1Loc-sizeBCD)==t2Loc && t2Loc<t3Loc){
                         intMult = asymmInt5tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                         intMult += asymmInt5tog12(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                     }
                     if(t1Loc>=sizeBCD && (t1Loc-sizeBCD)<t2Loc && t2Loc==t3Loc){
                         intMult = asymmInt5tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && (t1Loc-sizeBCD)>t2Loc && (t1Loc-sizeBCD)==t3Loc){
                         intMult = asymmInt5tog23(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                     }
                     if(t1Loc>=sizeBCD && (t1Loc-sizeBCD)==t2Loc && t2Loc==t3Loc){
                         intMult = asymmInt5togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                         intMult += asymmInt5togAll(tau1, tau2, tau3, psi2a, psi2b, psi1a, psi1b, psi3a, psi3b, g2, g1, g3);
                     }
                  //   intVec.push_back(intMult);
                     for(int pattern=0; pattern<15; ++pattern){
                        int index1 {permMat[pattern][2]-1};
                        probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                     }
                 }
              //   intMat.push_back(intVec);
               //  intVec.clear();
             }
         }
     }
    
    // (A,(B,(C,D)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
          psi3a = rootMat[t3Loc][0];
          psi3b = rootMat[t3Loc][1];
         for(int t2Loc=0; (t2Loc-sizeBCD)<=t3Loc; ++t2Loc){
             psi2a = BCDmatFull[t2Loc][0];
             psi2b = BCDmatFull[t2Loc][1];
             for(int t1Loc=0; (t1Loc-sizeCD)<=t2Loc; ++t1Loc){
                 psi1a = CDmatFull[t1Loc][0];
                 psi1b = CDmatFull[t1Loc][1];
                 d3 = tau1*(gammaC-CDmat[0][2]);
                 d4 = tau1*(gammaD-CDmat[0][2]);
                 d1 = tau3*(gammaA-rootMat[0][2]);
                 d2 = tau2*(gammaB-BCDmat[0][2]);
                 d5 = 0;
                 if(t1Loc !=0){
                     for(int i=1; i<=t1Loc; ++i){
                         d3 += CDmatFull[i][0]*(CDmatFull[i-1][2]-CDmatFull[i][2]);
                         d4 += CDmatFull[i][0]*(CDmatFull[i-1][2]-CDmatFull[i][2]);
                     }
                 }
                 if(t2Loc !=0){
                     for(int i=1; i<=t2Loc; ++i){
                     d2 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     }
                 }
                 if(t3Loc !=0){
                     for(int i=1; i<=t3Loc; ++i){
                     d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if((t1Loc-sizeCD) < t2Loc){
                     for(int i=(t1Loc+1); i<=(t2Loc+sizeCD); ++i){
                         d5 += CDmatFull[i][0]*(CDmatFull[i-1][2]-CDmatFull[i][2]);
                     }
                 }
                 if((t2Loc-sizeBCD) < t3Loc){
                     for(int i=(t2Loc+1); i<=(t3Loc+sizeBCD); ++i){
                         d1 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     }
                 }
                 for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                     g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                     g1 = thetaNew*(CDmatFull[t1Loc][2]*distMat[distIndex][2] + CDmatFull[t1Loc][2]*distMat[distIndex][3] - CDmatFull[t1Loc][2]*distMat[distIndex][4]);
                     g2 = thetaNew*(BCDmatFull[t2Loc][2]*distMat[distIndex][1] + BCDmatFull[t2Loc][2]*distMat[distIndex][4] - BCDmatFull[t2Loc][2]*distMat[distIndex][0]);
                     g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                     long double intMult{};
                     if(t1Loc<sizeCD  && t2Loc<sizeBCD){
                         intMult = asymmInt1(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc<sizeCD  && t2Loc>=sizeBCD && (t2Loc-sizeBCD)<t3Loc){
                         intMult = asymmInt2sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc<sizeCD  && t2Loc>=sizeBCD && (t2Loc-sizeBCD)==t3Loc){
                         intMult = asymmInt2tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeCD  && t2Loc<sizeBCD && (t1Loc-sizeCD)<t2Loc){
                         intMult = asymmInt3sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeCD  && t2Loc<sizeBCD && (t1Loc-sizeCD)==t2Loc){
                         intMult = asymmInt3tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeCD && t1Loc<(sizeCD+sizeBCD) && t2Loc>=sizeBCD && (t2Loc-sizeBCD)<t3Loc){
                         intMult = asymmInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeCD && t1Loc<(sizeCD+sizeBCD) && t2Loc>=sizeBCD && (t2Loc-sizeBCD)==t3Loc){
                         intMult = asymmInt4tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=(sizeCD+sizeBCD) && t2Loc>=sizeBCD && (t1Loc-sizeCD)<t2Loc && (t2Loc-sizeBCD)<t3Loc){
                         intMult = asymmInt5sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=(sizeCD+sizeBCD) && t2Loc>=sizeBCD && (t1Loc-sizeCD)==t2Loc && (t2Loc-sizeBCD)<t3Loc){
                         intMult = asymmInt5tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=(sizeCD+sizeBCD) && t2Loc>=sizeBCD && (t1Loc-sizeCD)<t2Loc && (t2Loc-sizeBCD)==t3Loc){
                         intMult = asymmInt5tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=(sizeCD+sizeBCD) && t2Loc>=sizeBCD && (t1Loc-sizeCD)==t2Loc && (t2Loc-sizeBCD)==t3Loc){
                         intMult = asymmInt5togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                  //   intVec.push_back(intMult);
                     for(int pattern=0; pattern<15; ++pattern){
                        probs[pattern] += intMult*g0*coefMat[pattern][distIndex]/256.0;
                     }
                 }
              //   intMat.push_back(intVec);
               //  intVec.clear();
             }
         }
     }
    
    // (A,(C,(B,D)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
          psi3a = rootMat[t3Loc][0];
          psi3b = rootMat[t3Loc][1];
         for(int t2Loc=0; (t2Loc-sizeBCD)<=t3Loc; ++t2Loc){
             psi2a = BCDmatFull[t2Loc][0];
             psi2b = BCDmatFull[t2Loc][1];
             for(int t1Loc=0; t1Loc<=t2Loc; ++t1Loc){
                 psi1a = BCDmatFull[t1Loc][0];
                 psi1b = BCDmatFull[t1Loc][1];
                 d2 = tau1*(gammaC-CDmat[0][2]);
                 d4 = tau1*(gammaD-CDmat[0][2]);
                 d1 = tau3*(gammaA-rootMat[0][2]);
                 d3 = tau2*(gammaB-BCDmat[0][2]);
                 d5 = 0;
                 d2 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                 d4 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                 if(sizeCD>1){
                     for(int i=1; i<sizeCD; ++i){
                         d2 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                         d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                     }
                 }
                 if(t1Loc !=0){
                     for(int i=1; i<=t1Loc; ++i){
                     d3 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     d4 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     }
                 }
                 if(t2Loc !=0){
                     for(int i=1; i<=t2Loc; ++i){
                     d2 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     }
                 }
                 if(t3Loc !=0){
                     for(int i=1; i<=t3Loc; ++i){
                     d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t1Loc < t2Loc){
                     for(int i=(t1Loc+1); i<=t2Loc; ++i){
                         d5 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     }
                 }
                 if((t2Loc-sizeBCD) < t3Loc){
                     for(int i=(t2Loc+1); i<=(t3Loc+sizeBCD); ++i){
                         d1 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     }
                 }
                 for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                     g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                     g1 = thetaNew*(BCDmatFull[t1Loc][2]*distMat[distIndex][2] + BCDmatFull[t1Loc][2]*distMat[distIndex][3] - BCDmatFull[t1Loc][2]*distMat[distIndex][4]);
                     g2 = thetaNew*(BCDmatFull[t2Loc][2]*distMat[distIndex][1] + BCDmatFull[t2Loc][2]*distMat[distIndex][4] - BCDmatFull[t2Loc][2]*distMat[distIndex][0]);
                     g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                     long double intMult{};
                     if(t1Loc<sizeBCD  && t2Loc<sizeBCD && t1Loc<t2Loc){
                         intMult = asymmInt3sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc<sizeBCD  && t2Loc<sizeBCD && t1Loc==t2Loc){
                         intMult = asymmInt3tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc<sizeBCD && t2Loc>=sizeBCD && (t2Loc-sizeBCD)<t3Loc){
                         intMult = asymmInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc<sizeBCD && t2Loc>=sizeBCD && (t2Loc-sizeBCD)==t3Loc){
                         intMult = asymmInt4tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && t1Loc<t2Loc && (t2Loc-sizeBCD)<t3Loc){
                         intMult = asymmInt5sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && t1Loc==t2Loc && (t2Loc-sizeBCD)<t3Loc){
                         intMult = asymmInt5tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && t1Loc<t2Loc && (t2Loc-sizeBCD)==t3Loc){
                         intMult = asymmInt5tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && t1Loc==t2Loc && (t2Loc-sizeBCD)==t3Loc){
                         intMult = asymmInt5togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                  //   intVec.push_back(intMult);
                     for(int pattern=0; pattern<15; ++pattern){
                        int index1 {permMat[pattern][1]-1};
                        probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                     }
                 }
              //   intMat.push_back(intVec);
               //  intVec.clear();
             }
         }
     }
    
    // (A,(D,(B,C)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
          psi3a = rootMat[t3Loc][0];
          psi3b = rootMat[t3Loc][1];
         for(int t2Loc=0; (t2Loc-sizeBCD)<=t3Loc; ++t2Loc){
             psi2a = BCDmatFull[t2Loc][0];
             psi2b = BCDmatFull[t2Loc][1];
             for(int t1Loc=0; t1Loc<=t2Loc; ++t1Loc){
                 psi1a = BCDmatFull[t1Loc][0];
                 psi1b = BCDmatFull[t1Loc][1];
                 d4 = tau1*(gammaC-CDmat[0][2]);
                 d2 = tau1*(gammaD-CDmat[0][2]);
                 d1 = tau3*(gammaA-rootMat[0][2]);
                 d3 = tau2*(gammaB-BCDmat[0][2]);
                 d5 = 0;
                 d2 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                 d4 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                 if(sizeCD>1){
                     for(int i=1; i<sizeCD; ++i){
                         d2 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                         d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                     }
                 }
                 if(t1Loc !=0){
                     for(int i=1; i<=t1Loc; ++i){
                     d3 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     d4 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     }
                 }
                 if(t2Loc !=0){
                     for(int i=1; i<=t2Loc; ++i){
                     d2 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     }
                 }
                 if(t3Loc !=0){
                     for(int i=1; i<=t3Loc; ++i){
                     d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t1Loc < t2Loc){
                     for(int i=(t1Loc+1); i<=t2Loc; ++i){
                         d5 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     }
                 }
                 if((t2Loc-sizeBCD) < t3Loc){
                     for(int i=(t2Loc+1); i<=(t3Loc+sizeBCD); ++i){
                         d1 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     }
                 }
                 for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                     g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                     g1 = thetaNew*(BCDmatFull[t1Loc][2]*distMat[distIndex][2] + BCDmatFull[t1Loc][2]*distMat[distIndex][3] - BCDmatFull[t1Loc][2]*distMat[distIndex][4]);
                     g2 = thetaNew*(BCDmatFull[t2Loc][2]*distMat[distIndex][1] + BCDmatFull[t2Loc][2]*distMat[distIndex][4] - BCDmatFull[t2Loc][2]*distMat[distIndex][0]);
                     g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                     long double intMult{};
                     if(t1Loc<sizeBCD  && t2Loc<sizeBCD && t1Loc<t2Loc){
                         intMult = asymmInt3sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc<sizeBCD  && t2Loc<sizeBCD && t1Loc==t2Loc){
                         intMult = asymmInt3tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc<sizeBCD && t2Loc>=sizeBCD && (t2Loc-sizeBCD)<t3Loc){
                         intMult = asymmInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc<sizeBCD && t2Loc>=sizeBCD && (t2Loc-sizeBCD)==t3Loc){
                         intMult = asymmInt4tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && t1Loc<t2Loc && (t2Loc-sizeBCD)<t3Loc){
                         intMult = asymmInt5sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && t1Loc==t2Loc && (t2Loc-sizeBCD)<t3Loc){
                         intMult = asymmInt5tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && t1Loc<t2Loc && (t2Loc-sizeBCD)==t3Loc){
                         intMult = asymmInt5tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && t1Loc==t2Loc && (t2Loc-sizeBCD)==t3Loc){
                         intMult = asymmInt5togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                  //   intVec.push_back(intMult);
                     for(int pattern=0; pattern<15; ++pattern){
                        int index1 {permMat[pattern][2]-1};
                        probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                     }
                 }
              //   intMat.push_back(intVec);
               //  intVec.clear();
             }
         }
     }
    
    // (B,(A,(C,D)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
          psi3a = rootMat[t3Loc][0];
          psi3b = rootMat[t3Loc][1];
         for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
             psi2a = rootMat[t2Loc][0];
             psi2b = rootMat[t2Loc][1];
             for(int t1Loc=0; (t1Loc-sizeCD-sizeBCD)<=t2Loc; ++t1Loc){
                 psi1a = CDmatFull[t1Loc][0];
                 psi1b = CDmatFull[t1Loc][1];
                 d3 = tau1*(gammaC-CDmat[0][2]);
                 d4 = tau1*(gammaD-CDmat[0][2]);
                 d2 = tau3*(gammaA-rootMat[0][2]);
                 d1 = tau2*(gammaB-BCDmat[0][2]);
                 d5 = 0;
                 d1 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                 if(sizeBCD>1){
                     for(int i=1; i<sizeBCD; ++i){
                         d1 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                     }
                 }
                 if(t1Loc !=0){
                     for(int i=1; i<=t1Loc; ++i){
                     d3 += CDmatFull[i][0]*(CDmatFull[i-1][2]-CDmatFull[i][2]);
                     d4 += CDmatFull[i][0]*(CDmatFull[i-1][2]-CDmatFull[i][2]);
                     }
                 }
                 if(t2Loc !=0){
                     for(int i=1; i<=t2Loc; ++i){
                     d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t3Loc !=0){
                     for(int i=1; i<=t3Loc; ++i){
                     d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if((t1Loc-sizeCD-sizeBCD) < t2Loc){
                     for(int i=(t1Loc+1); i<=(t2Loc+sizeCD+sizeBCD); ++i){
                         d5 += CDmatFull[i][0]*(CDmatFull[i-1][2]-CDmatFull[i][2]);
                     }
                 }
                 if(t2Loc < t3Loc){
                     for(int i=(t2Loc+1); i<=t3Loc; ++i){
                         d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                     g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                     g1 = thetaNew*(CDmatFull[t1Loc][2]*distMat[distIndex][2] + CDmatFull[t1Loc][2]*distMat[distIndex][3] - CDmatFull[t1Loc][2]*distMat[distIndex][4]);
                     g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                     g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                     long double intMult{};
                     if(t1Loc<sizeCD && t2Loc<t3Loc){
                         intMult = asymmInt2sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc<sizeCD && t2Loc==t3Loc){
                         intMult = asymmInt2tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeCD && t1Loc<(sizeCD+sizeBCD) && t2Loc<t3Loc){
                         intMult = asymmInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeCD && t1Loc<(sizeCD+sizeBCD) && t2Loc==t3Loc){
                         intMult = asymmInt4tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=(sizeCD+sizeBCD) && (t1Loc-sizeCD-sizeBCD)<t2Loc && t2Loc<t3Loc){
                         intMult = asymmInt5sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=(sizeCD+sizeBCD) && (t1Loc-sizeCD-sizeBCD)==t2Loc && t2Loc<t3Loc){
                         intMult = asymmInt5tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=(sizeCD+sizeBCD) && (t1Loc-sizeCD-sizeBCD)<t2Loc && t2Loc==t3Loc){
                         intMult = asymmInt5tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=(sizeCD+sizeBCD) && (t1Loc-sizeCD-sizeBCD)==t2Loc && t2Loc==t3Loc){
                         intMult = asymmInt5togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                  //   intVec.push_back(intMult);
                     for(int pattern=0; pattern<15; ++pattern){
                        int index1 {permMat[pattern][3]-1};
                        probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                     }
                 }
              //   intMat.push_back(intVec);
               //  intVec.clear();
             }
         }
     }
    
    // (B,(C,(A,D)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
          psi3a = rootMat[t3Loc][0];
          psi3b = rootMat[t3Loc][1];
         for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
             psi2a = rootMat[t2Loc][0];
             psi2b = rootMat[t2Loc][1];
             for(int t1Loc=0; t1Loc<=t2Loc; ++t1Loc){
                 psi1a = rootMat[t1Loc][0];
                 psi1b = rootMat[t1Loc][1];
                 d2 = tau1*(gammaC-CDmat[0][2]);
                 d4 = tau1*(gammaD-CDmat[0][2]);
                 d3 = tau3*(gammaA-rootMat[0][2]);
                 d1 = tau2*(gammaB-BCDmat[0][2]);
                 d5 = 0;
                 d1 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                 d2 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                 d4 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                 if(sizeBCD>1){
                     for(int i=1; i<sizeBCD; ++i){
                         d1 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                         d2 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                         d4 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                     }
                 }
                 d2 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                 d4 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                 if(sizeCD>1){
                     for(int i=1; i<sizeCD; ++i){
                         d2 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                         d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                     }
                 }
                 if(t1Loc !=0){
                     for(int i=1; i<=t1Loc; ++i){
                     d3 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     d4 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t2Loc !=0){
                     for(int i=1; i<=t2Loc; ++i){
                     d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t3Loc !=0){
                     for(int i=1; i<=t3Loc; ++i){
                     d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t1Loc < t2Loc){
                     for(int i=(t1Loc+1); i<=t2Loc; ++i){
                         d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t2Loc < t3Loc){
                     for(int i=(t2Loc+1); i<=t3Loc; ++i){
                         d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                     g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                     g1 = thetaNew*(rootMat[t1Loc][2]*distMat[distIndex][2] + rootMat[t1Loc][2]*distMat[distIndex][3] - rootMat[t1Loc][2]*distMat[distIndex][4]);
                     g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                     g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                     long double intMult{};
                     if( t1Loc<t2Loc && t2Loc<t3Loc){
                         intMult = asymmInt5sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if( t1Loc==t2Loc && t2Loc<t3Loc){
                         intMult = asymmInt5tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if( t1Loc<t2Loc && t2Loc==t3Loc){
                         intMult = asymmInt5tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if( t1Loc==t2Loc && t2Loc==t3Loc){
                         intMult = asymmInt5togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                  //   intVec.push_back(intMult);
                     for(int pattern=0; pattern<15; ++pattern){
                        int index1 {permMat[pattern][4]-1};
                        probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                     }
                 }
              //   intMat.push_back(intVec);
               //  intVec.clear();
             }
         }
     }
    
    // (B,(D,(A,C)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
          psi3a = rootMat[t3Loc][0];
          psi3b = rootMat[t3Loc][1];
         for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
             psi2a = rootMat[t2Loc][0];
             psi2b = rootMat[t2Loc][1];
             for(int t1Loc=0; t1Loc<=t2Loc; ++t1Loc){
                 psi1a = rootMat[t1Loc][0];
                 psi1b = rootMat[t1Loc][1];
                 d4 = tau1*(gammaC-CDmat[0][2]);
                 d2 = tau1*(gammaD-CDmat[0][2]);
                 d3 = tau3*(gammaA-rootMat[0][2]);
                 d1 = tau2*(gammaB-BCDmat[0][2]);
                 d5 = 0;
                 d1 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                 d2 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                 d4 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                 d2 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                 d4 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                 if(sizeBCD>1){
                     for(int i=1; i<sizeBCD; ++i){
                         d1 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                         d2 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                         d4 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                     }
                 }
                 if(sizeCD>1){
                     for(int i=1; i<sizeCD; ++i){
                         d2 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                         d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                     }
                 }
                 if(t1Loc !=0){
                     for(int i=1; i<=t1Loc; ++i){
                     d3 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     d4 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t2Loc !=0){
                     for(int i=1; i<=t2Loc; ++i){
                     d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t3Loc !=0){
                     for(int i=1; i<=t3Loc; ++i){
                     d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t1Loc < t2Loc){
                     for(int i=(t1Loc+1); i<=t2Loc; ++i){
                         d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t2Loc < t3Loc){
                     for(int i=(t2Loc+1); i<=t3Loc; ++i){
                         d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                     g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                     g1 = thetaNew*(rootMat[t1Loc][2]*distMat[distIndex][2] + rootMat[t1Loc][2]*distMat[distIndex][3] - rootMat[t1Loc][2]*distMat[distIndex][4]);
                     g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                     g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                     long double intMult{};
                     if( t1Loc<t2Loc && t2Loc<t3Loc){
                         intMult = asymmInt5sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if( t1Loc==t2Loc && t2Loc<t3Loc){
                         intMult = asymmInt5tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if( t1Loc<t2Loc && t2Loc==t3Loc){
                         intMult = asymmInt5tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if( t1Loc==t2Loc && t2Loc==t3Loc){
                         intMult = asymmInt5togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                  //   intVec.push_back(intMult);
                     for(int pattern=0; pattern<15; ++pattern){
                        int index1 {permMat[pattern][5]-1};
                        probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                     }
                 }
              //   intMat.push_back(intVec);
               //  intVec.clear();
             }
         }
     }
    
    // (C,(A,(B,D)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
          psi3a = rootMat[t3Loc][0];
          psi3b = rootMat[t3Loc][1];
         for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
             psi2a = rootMat[t2Loc][0];
             psi2b = rootMat[t2Loc][1];
             for(int t1Loc=0; (t1Loc-sizeBCD)<=t2Loc; ++t1Loc){
                 psi1a = BCDmatFull[t1Loc][0];
                 psi1b = BCDmatFull[t1Loc][1];
                 d1 = tau1*(gammaC-CDmat[0][2]);
                 d4 = tau1*(gammaD-CDmat[0][2]);
                 d2 = tau3*(gammaA-rootMat[0][2]);
                 d3 = tau2*(gammaB-BCDmat[0][2]);
                 d5 = 0;
                 d1 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                 d1 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                 d4 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                 if(sizeBCD>1){
                     for(int i=1; i<sizeBCD; ++i){
                         d1 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                     }
                 }
                 if(sizeCD>1){
                     for(int i=1; i<sizeCD; ++i){
                         d1 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                         d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                     }
                 }
                 if(t1Loc !=0){
                     for(int i=1; i<=t1Loc; ++i){
                     d3 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     d4 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     }
                 }
                 if(t2Loc !=0){
                     for(int i=1; i<=t2Loc; ++i){
                     d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t3Loc !=0){
                     for(int i=1; i<=t3Loc; ++i){
                     d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if((t1Loc-sizeBCD) < t2Loc){
                     for(int i=(t1Loc+1); i<=(t2Loc+sizeBCD); ++i){
                         d5 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     }
                 }
                 if(t2Loc < t3Loc){
                     for(int i=(t2Loc+1); i<=t3Loc; ++i){
                         d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                     g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                     g1 = thetaNew*(BCDmatFull[t1Loc][2]*distMat[distIndex][2] + BCDmatFull[t1Loc][2]*distMat[distIndex][3] - BCDmatFull[t1Loc][2]*distMat[distIndex][4]);
                     g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                     g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                     long double intMult{};
                     if(t1Loc<sizeBCD && t2Loc<t3Loc){
                         intMult = asymmInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc<sizeBCD && t2Loc==t3Loc){
                         intMult = asymmInt4tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && (t1Loc-sizeBCD)<t2Loc && t2Loc<t3Loc){
                         intMult = asymmInt5sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && (t1Loc-sizeBCD)==t2Loc && t2Loc<t3Loc){
                         intMult = asymmInt5tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && (t1Loc-sizeBCD)<t2Loc && t2Loc==t3Loc){
                         intMult = asymmInt5tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if( t1Loc>=sizeBCD && (t1Loc-sizeBCD)==t2Loc && t2Loc==t3Loc){
                         intMult = asymmInt5togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                  //   intVec.push_back(intMult);
                     for(int pattern=0; pattern<15; ++pattern){
                        int index1 {permMat[pattern][6]-1};
                        probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                     }
                 }
              //   intMat.push_back(intVec);
               //  intVec.clear();
             }
         }
     }
    
    // (C,(B,(A,D)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
             psi3a = rootMat[t3Loc][0];
             psi3b = rootMat[t3Loc][1];
            for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
                psi2a = rootMat[t2Loc][0];
                psi2b = rootMat[t2Loc][1];
                for(int t1Loc=0; t1Loc<=t2Loc; ++t1Loc){
                    psi1a = rootMat[t1Loc][0];
                    psi1b = rootMat[t1Loc][1];
                    d1 = tau1*(gammaC-CDmat[0][2]);
                    d4 = tau1*(gammaD-CDmat[0][2]);
                    d3 = tau3*(gammaA-rootMat[0][2]);
                    d2 = tau2*(gammaB-BCDmat[0][2]);
                    d5 = 0;
                    d1 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                    d2 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                    d4 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                    d1 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                    d4 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                    if(sizeBCD>1){
                        for(int i=1; i<sizeBCD; ++i){
                            d1 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                            d2 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                            d4 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                        }
                    }
                    if(sizeCD>1){
                        for(int i=1; i<sizeCD; ++i){
                            d1 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                            d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                        }
                    }
                    if(t1Loc !=0){
                        for(int i=1; i<=t1Loc; ++i){
                        d3 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        d4 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    if(t2Loc !=0){
                        for(int i=1; i<=t2Loc; ++i){
                        d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    if(t3Loc !=0){
                        for(int i=1; i<=t3Loc; ++i){
                        d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    if(t1Loc < t2Loc){
                        for(int i=(t1Loc+1); i<=t2Loc; ++i){
                            d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    if(t2Loc < t3Loc){
                        for(int i=(t2Loc+1); i<=t3Loc; ++i){
                            d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                        g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                        g1 = thetaNew*(rootMat[t1Loc][2]*distMat[distIndex][2] + rootMat[t1Loc][2]*distMat[distIndex][3] - rootMat[t1Loc][2]*distMat[distIndex][4]);
                        g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                        g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                        long double intMult{};
                        if( t1Loc<t2Loc && t2Loc<t3Loc){
                            intMult = asymmInt5sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        }
                        if( t1Loc==t2Loc && t2Loc<t3Loc){
                            intMult = asymmInt5tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        }
                        if( t1Loc<t2Loc && t2Loc==t3Loc){
                            intMult = asymmInt5tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        }
                        if( t1Loc==t2Loc && t2Loc==t3Loc){
                            intMult = asymmInt5togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        }
                     //   intVec.push_back(intMult);
                        for(int pattern=0; pattern<15; ++pattern){
                           int index1 {permMat[pattern][7]-1};
                           probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                        }
                    }
                 //   intMat.push_back(intVec);
                  //  intVec.clear();
                }
            }
        }
    
    // (C,(D,(A,B)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
             psi3a = rootMat[t3Loc][0];
             psi3b = rootMat[t3Loc][1];
            for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
                psi2a = rootMat[t2Loc][0];
                psi2b = rootMat[t2Loc][1];
                for(int t1Loc=0; t1Loc<=t2Loc; ++t1Loc){
                    psi1a = rootMat[t1Loc][0];
                    psi1b = rootMat[t1Loc][1];
                    d1 = tau1*(gammaC-CDmat[0][2]);
                    d2 = tau1*(gammaD-CDmat[0][2]);
                    d3 = tau3*(gammaA-rootMat[0][2]);
                    d4 = tau2*(gammaB-BCDmat[0][2]);
                    d5 = 0;
                    d1 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                    d2 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                    d4 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                    d1 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                    d2 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                    if(sizeBCD>1){
                        for(int i=1; i<sizeBCD; ++i){
                            d1 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                            d2 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                            d4 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                        }
                    }
                    if(sizeCD>1){
                        for(int i=1; i<sizeCD; ++i){
                            d1 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                            d2 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                        }
                    }
                    if(t1Loc !=0){
                        for(int i=1; i<=t1Loc; ++i){
                        d3 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        d4 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    if(t2Loc !=0){
                        for(int i=1; i<=t2Loc; ++i){
                        d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    if(t3Loc !=0){
                        for(int i=1; i<=t3Loc; ++i){
                        d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    if(t1Loc < t2Loc){
                        for(int i=(t1Loc+1); i<=t2Loc; ++i){
                            d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    if(t2Loc < t3Loc){
                        for(int i=(t2Loc+1); i<=t3Loc; ++i){
                            d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                        g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                        g1 = thetaNew*(rootMat[t1Loc][2]*distMat[distIndex][2] + rootMat[t1Loc][2]*distMat[distIndex][3] - rootMat[t1Loc][2]*distMat[distIndex][4]);
                        g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                        g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                        long double intMult{};
                        if( t1Loc<t2Loc && t2Loc<t3Loc){
                            intMult = asymmInt5sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        }
                        if( t1Loc==t2Loc && t2Loc<t3Loc){
                            intMult = asymmInt5tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        }
                        if( t1Loc<t2Loc && t2Loc==t3Loc){
                            intMult = asymmInt5tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        }
                        if( t1Loc==t2Loc && t2Loc==t3Loc){
                            intMult = asymmInt5togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        }
                     //   intVec.push_back(intMult);
                        for(int pattern=0; pattern<15; ++pattern){
                           int index1 {permMat[pattern][8]-1};
                           probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                        }
                    }
                 //   intMat.push_back(intVec);
                  //  intVec.clear();
                }
            }
        }
    
    // (D,(A,(B,C)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
          psi3a = rootMat[t3Loc][0];
          psi3b = rootMat[t3Loc][1];
         for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
             psi2a = rootMat[t2Loc][0];
             psi2b = rootMat[t2Loc][1];
             for(int t1Loc=0; (t1Loc-sizeBCD)<=t2Loc; ++t1Loc){
                 psi1a = BCDmatFull[t1Loc][0];
                 psi1b = BCDmatFull[t1Loc][1];
                 d4 = tau1*(gammaC-CDmat[0][2]);
                 d1 = tau1*(gammaD-CDmat[0][2]);
                 d2 = tau3*(gammaA-rootMat[0][2]);
                 d3 = tau2*(gammaB-BCDmat[0][2]);
                 d5 = 0;
                 d1 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                 d1 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                 d4 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                 if(sizeBCD>1){
                     for(int i=1; i<sizeBCD; ++i){
                         d1 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                     }
                 }
                 if(sizeCD>1){
                     for(int i=1; i<sizeCD; ++i){
                         d1 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                         d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                     }
                 }
                 if(t1Loc !=0){
                     for(int i=1; i<=t1Loc; ++i){
                     d3 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     d4 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     }
                 }
                 if(t2Loc !=0){
                     for(int i=1; i<=t2Loc; ++i){
                     d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if(t3Loc !=0){
                     for(int i=1; i<=t3Loc; ++i){
                     d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 if((t1Loc-sizeBCD) < t2Loc){
                     for(int i=(t1Loc+1); i<=(t2Loc+sizeBCD); ++i){
                         d5 += BCDmatFull[i][0]*(BCDmatFull[i-1][2]-BCDmatFull[i][2]);
                     }
                 }
                 if(t2Loc < t3Loc){
                     for(int i=(t2Loc+1); i<=t3Loc; ++i){
                         d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                     }
                 }
                 for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                     g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                     g1 = thetaNew*(BCDmatFull[t1Loc][2]*distMat[distIndex][2] + BCDmatFull[t1Loc][2]*distMat[distIndex][3] - BCDmatFull[t1Loc][2]*distMat[distIndex][4]);
                     g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                     g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                     long double intMult{};
                     if(t1Loc<sizeBCD && t2Loc<t3Loc){
                         intMult = asymmInt4sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc<sizeBCD && t2Loc==t3Loc){
                         intMult = asymmInt4tog(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && (t1Loc-sizeBCD)<t2Loc && t2Loc<t3Loc){
                         intMult = asymmInt5sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && (t1Loc-sizeBCD)==t2Loc && t2Loc<t3Loc){
                         intMult = asymmInt5tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if(t1Loc>=sizeBCD && (t1Loc-sizeBCD)<t2Loc && t2Loc==t3Loc){
                         intMult = asymmInt5tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                     if( t1Loc>=sizeBCD && (t1Loc-sizeBCD)==t2Loc && t2Loc==t3Loc){
                         intMult = asymmInt5togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                     }
                  //   intVec.push_back(intMult);
                     for(int pattern=0; pattern<15; ++pattern){
                        int index1 {permMat[pattern][9]-1};
                        probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                     }
                 }
              //   intMat.push_back(intVec);
               //  intVec.clear();
             }
         }
     }
    
    // (D,(B,(A,C)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
             psi3a = rootMat[t3Loc][0];
             psi3b = rootMat[t3Loc][1];
            for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
                psi2a = rootMat[t2Loc][0];
                psi2b = rootMat[t2Loc][1];
                for(int t1Loc=0; t1Loc<=t2Loc; ++t1Loc){
                    psi1a = rootMat[t1Loc][0];
                    psi1b = rootMat[t1Loc][1];
                    d4 = tau1*(gammaC-CDmat[0][2]);
                    d1 = tau1*(gammaD-CDmat[0][2]);
                    d3 = tau3*(gammaA-rootMat[0][2]);
                    d2 = tau2*(gammaB-BCDmat[0][2]);
                    d5 = 0;
                    d1 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                    d2 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                    d4 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                    d1 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                    d4 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                    if(sizeBCD>1){
                        for(int i=1; i<sizeBCD; ++i){
                            d1 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                            d2 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                            d4 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                        }
                    }
                    if(sizeCD>1){
                        for(int i=1; i<sizeCD; ++i){
                            d1 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                            d4 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                        }
                    }
                    if(t1Loc !=0){
                        for(int i=1; i<=t1Loc; ++i){
                        d3 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        d4 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    if(t2Loc !=0){
                        for(int i=1; i<=t2Loc; ++i){
                        d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    if(t3Loc !=0){
                        for(int i=1; i<=t3Loc; ++i){
                        d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    if(t1Loc < t2Loc){
                        for(int i=(t1Loc+1); i<=t2Loc; ++i){
                            d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    if(t2Loc < t3Loc){
                        for(int i=(t2Loc+1); i<=t3Loc; ++i){
                            d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                        g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                        g1 = thetaNew*(rootMat[t1Loc][2]*distMat[distIndex][2] + rootMat[t1Loc][2]*distMat[distIndex][3] - rootMat[t1Loc][2]*distMat[distIndex][4]);
                        g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                        g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                        long double intMult{};
                        if( t1Loc<t2Loc && t2Loc<t3Loc){
                            intMult = asymmInt5sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        }
                        if( t1Loc==t2Loc && t2Loc<t3Loc){
                            intMult = asymmInt5tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        }
                        if( t1Loc<t2Loc && t2Loc==t3Loc){
                            intMult = asymmInt5tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        }
                        if( t1Loc==t2Loc && t2Loc==t3Loc){
                            intMult = asymmInt5togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        }
                     //   intVec.push_back(intMult);
                        for(int pattern=0; pattern<15; ++pattern){
                           int index1 {permMat[pattern][10]-1};
                           probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                        }
                    }
                 //   intMat.push_back(intVec);
                  //  intVec.clear();
                }
            }
        }
    
    // (D,(C,(A,B)))
    for(int t3Loc=0; t3Loc<rootMat.size(); ++t3Loc){
             psi3a = rootMat[t3Loc][0];
             psi3b = rootMat[t3Loc][1];
            for(int t2Loc=0; t2Loc<=t3Loc; ++t2Loc){
                psi2a = rootMat[t2Loc][0];
                psi2b = rootMat[t2Loc][1];
                for(int t1Loc=0; t1Loc<=t2Loc; ++t1Loc){
                    psi1a = rootMat[t1Loc][0];
                    psi1b = rootMat[t1Loc][1];
                    d2 = tau1*(gammaC-CDmat[0][2]);
                    d1 = tau1*(gammaD-CDmat[0][2]);
                    d3 = tau3*(gammaA-rootMat[0][2]);
                    d4 = tau2*(gammaB-BCDmat[0][2]);
                    d5 = 0;
                    d1 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                    d2 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                    d4 += tau3*(BCDmat[sizeBCD-1][2]-rootMat[0][2]);
                    d1 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                    d2 += tau2*(CDmat[sizeCD-1][2]-BCDmat[0][2]);
                    if(sizeBCD>1){
                        for(int i=1; i<sizeBCD; ++i){
                            d1 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                            d2 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                            d4 += BCDmat[i][0]*(BCDmat[i-1][2]-BCDmat[i][2]);
                        }
                    }
                    if(sizeCD>1){
                        for(int i=1; i<sizeCD; ++i){
                            d1 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                            d2 += CDmat[i][0]*(CDmat[i-1][2]-CDmat[i][2]);
                        }
                    }
                    if(t1Loc !=0){
                        for(int i=1; i<=t1Loc; ++i){
                        d3 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        d4 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    if(t2Loc !=0){
                        for(int i=1; i<=t2Loc; ++i){
                        d2 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    if(t3Loc !=0){
                        for(int i=1; i<=t3Loc; ++i){
                        d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    if(t1Loc < t2Loc){
                        for(int i=(t1Loc+1); i<=t2Loc; ++i){
                            d5 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    if(t2Loc < t3Loc){
                        for(int i=(t2Loc+1); i<=t3Loc; ++i){
                            d1 += rootMat[i][0]*(rootMat[i-1][2]-rootMat[i][2]);
                        }
                    }
                    for(int distIndex=0; distIndex<distMat.size(); ++distIndex){
                        g0 = exp(-thetaNew*(d1*distMat[distIndex][0] + d2*distMat[distIndex][1] + d3*distMat[distIndex][2] + d4*distMat[distIndex][3] + d5*distMat[distIndex][4]));
                        g1 = thetaNew*(rootMat[t1Loc][2]*distMat[distIndex][2] + rootMat[t1Loc][2]*distMat[distIndex][3] - rootMat[t1Loc][2]*distMat[distIndex][4]);
                        g2 = thetaNew*(rootMat[t2Loc][2]*distMat[distIndex][1] + rootMat[t2Loc][2]*distMat[distIndex][4] - rootMat[t2Loc][2]*distMat[distIndex][0]);
                        g3 = 2.0*thetaNew*rootMat[t3Loc][2]*distMat[distIndex][0];
                        long double intMult{};
                        if( t1Loc<t2Loc && t2Loc<t3Loc){
                            intMult = asymmInt5sep(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        }
                        if( t1Loc==t2Loc && t2Loc<t3Loc){
                            intMult = asymmInt5tog12(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        }
                        if( t1Loc<t2Loc && t2Loc==t3Loc){
                            intMult = asymmInt5tog23(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        }
                        if( t1Loc==t2Loc && t2Loc==t3Loc){
                            intMult = asymmInt5togAll(tau1, tau2, tau3, psi1a, psi1b, psi2a, psi2b, psi3a, psi3b, g1, g2, g3);
                        }
                     //   intVec.push_back(intMult);
                        for(int pattern=0; pattern<15; ++pattern){
                           int index1 {permMat[pattern][11]-1};
                           probs[pattern] += intMult*g0*coefMat[index1][distIndex]/256.0;
                        }
                    }
                 //   intMat.push_back(intVec);
                  //  intVec.clear();
                }
            }
        }
    
    
        return probs;
    }

int main(int argc, const char * argv[]) {

  /*  double breakPoint{};
    std::stringstream input1{argv[1]};
    input1 >> breakPoint;
    
    double upper{};
    std::stringstream input2{argv[2]};
    input2 >> upper;
    
    double lower{};
    std::stringstream input3{argv
        [3]};
    input3 >> lower;
    
    doubleMat2 mat1{{1.0, 1.0+breakPoint, lower}, {1.0+breakPoint, 3.0, upper}};
    // doubleMat2 mat1{{1.0, 3.0, 2.0}};
   // doubleMat2 mat1{{1.0, 2.0, 2.0}};
    doubleMat2 mat2{{2.0, 3.0, 1.0}};
    doubleMat2 mat3{{3.0, 1000.0, 1.0}}; */
    
    /*
    doubleMat3 gammaAB{
    {
        {0.0, 1.0, 1.0}
    },
    {
        {0.0, 0.5, 1.0},
        {0.5, 1.0, 3.0}
    },
    {
        {0.0, 0.3, 3.0},
        {0.3, 0.6, 2.0},
        {0.6, 1.0, 1.0}
    },
    {
        {0.0, 0.25, 0.5},
        {0.25, 0.5, 2.5},
        {0.5, 0.75, 1.5},
        {0.75, 1.0, 3.5}
    }
    };
    
    doubleMat3 gammaCD{
     {
         {0.0, 2.0, 1.0}
     },
     {
         {0.0, 1.0, 1.0},
         {1.0, 2.0, 3.0}
     },
     {
         {0.0, 0.6, 3.0},
         {0.6, 1.4, 2.0},
         {1.4, 2.0, 1.0}
     },
     {
         {0.0, 0.5, 0.5},
         {0.5, 1.0, 2.5},
         {1.0, 1.5, 1.5},
         {1.5, 2.0, 3.5}
     }
     };
    
    doubleMat3 gammaRoot{
    {
        {0.0, 1000.0, 1.0}
    },
        {
        {0.0, 0.5, 3.0},
        {0.5, 1000.0, 1.0}
    },
        {
        {0.0, 0.5, 3.0},
        {0.5, 1.0, 2.0},
        {1.0, 1000.0, 1.0}
    },
        {
        {0.0, 0.5, 4.0},
        {0.5, 1.0, 3.0},
        {1.0, 1.5, 2.0},
        {1.5, 1000.0, 1.0}
    }
    };
    
    
    int ABmatInd{};
     std::stringstream input2{argv[1]};
     input2 >> ABmatInd;
     
     int CDmatInd{};
     std::stringstream input3{argv[2]};
     input3 >> CDmatInd;
     
     int rootMatInd{};
     std::stringstream input4{argv[3]};
     input4 >> rootMatInd;
     
     doubleMat2 CDmat{gammaCD[CDmatInd]};
      doubleMat2 ABmat{gammaAB[ABmatInd]};
    doubleMat2 rootMat{gammaRoot[rootMatInd]};
    
    for(int i=0; i<CDmat.size(); ++i){
        CDmat[i][0] += 1.0;
        CDmat[i][1] += 1.0;
    }
    
    for(int i=0; i<ABmat.size(); ++i){
        ABmat[i][0] += 2.0;
        ABmat[i][1] += 2.0;
    }
    
    for(int i=0; i<rootMat.size(); ++i){
        rootMat[i][0] += 3.0;
        rootMat[i][1] += 3.0;
    }

    std::vector<long double> vec1{probSym(1.0, 2.0, 3.0, 1.0, 1.0, 1.0, 1.0, ABmat, CDmat, rootMat, 0.006)};
 //  std::vector<long double> vec1{probAsymm(1.0, 2.0, 3.0, 1.0, 1.0, 1.0, 1.0, CDmat, BCDmat, rootMat, 0.006)};
    
    std::vector<long double> vec2{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    vec2[0]=4.0*vec1[0];
    for(int i=1; i<8; ++i){
        vec2[i]=12.0*vec1[i];
    }
    for(int i=8; i<15; ++i){
          vec2[i]=24.0*vec1[i];
      }
 */
    
    doubleMat2 CDmat{
   //     {1.0, 1.5, 1.6},
    //   {1.5, 2.0, 1.7},
        {2.0, 3.0, 1.8},
        {3.0, 4.0, 1.9},
     //   {1.8, 2.0, 2.0},
    //    {2.0, 2.25, 2.1},
    //    {2.25, 2.5, 2.2},
    //    {2.5, 2.75, 2.3},
     //   {2.75, 3.0, 2.4}
    };
    
    doubleMat2 BCDmat{
  //    {2.0, 3.0, 1.8},
   //   {3.0, 4.0, 1.9},
     {3.0, 3.5, 2.0},
     {3.5, 4.0, 2.1},
     //   {3.7, 5.0, 1.0},
     //   {2.0, 3.0, 8.0},
     //   {3.0, 4.0, 9.0}
    };
    
    doubleMat2 rootMat{
 //    {3.0, 3.3, 2.0},
  //    {3.3, 3.7, 2.1},
   //   {3.7, 4.0, 2.2},
        {4.0, 4.1, 2.2},
        {4.1, 4.2, 2.3},
        {4.2, 4.3, 2.4},
        {4.3, 1000.0, 1.0}
    };
    
   std::vector<long double> vec1{probSym(2.0, 3.0, 4.0, 1.1, 1.2, 1.3, 1.5, BCDmat, CDmat, rootMat, 0.003)};
    
  //  std::vector<long double> vec1{probAsymm(1.0, 2.0, 4.0, 1.2, 1.3, 1.4, 1.5, CDmat, BCDmat, rootMat, 0.003)};

  //   std::vector<double> asymmProbVec{probAsymm(1.0, 2.00, 3.0, 1.1, 1.0, 1.0, 1.2, 1.0, 1.1, 0.003)};
    
    std::vector<long double> vec2{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    vec2[0]=4.0*vec1[0];
    for(int i=1; i<8; ++i){
        vec2[i]=12.0*vec1[i];
    }
    for(int i=8; i<15; ++i){
          vec2[i]=24.0*vec1[i];
      }
    
    std::cout << '\n';
    for(int i=0; i<15; ++i){
        std::cout << vec2[i] << " ";
    }
    std::cout << '\n';
    double total{};
    for(int i=0; i<15; ++i){
        total += vec2[i];
      }
    std::cout << total << '\n';
    

 
 /*      std::ofstream out1("expected", std::ios::app);
    for(int i=0; i<15; ++i){
          out1 << vec2[i] << " ";
      }
      out1 << '\n';
   */
    
    return 0;
}
