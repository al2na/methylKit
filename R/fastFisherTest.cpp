#include <iostream>
#include <cmath>
#include <cstdlib>

void initLogFacs(double* logFacs, int n);
double logHypergeometricProb(double* logFacs, int a, int b, int c, int d);
double fastFishersExactTest(int a, int b, int c, int d);


extern "C" {
    void fastFisherTest(int* a, int* b, int* c, int* d, int* na, double* pval){
        for(int i=0; i < *na; ++i){
            pval[i] = fastFishersExactTest(a[i], b[i], c[i], d[i]);
        }
    }
}

void initLogFacs(double* logFacs, int n) {
    logFacs[0] = 0;
    for(int i=1; i < n+1; ++i) {
        logFacs[i] = logFacs[i-1] + log((double)i); // only n times of log() calls
    }
}

double logHypergeometricProb(double* logFacs, int a, int b, int c, int d) {
    return logFacs[a+b] + logFacs[c+d] + logFacs[a+c] + logFacs[b+d]
    - logFacs[a] - logFacs[b] - logFacs[c] - logFacs[d] - logFacs[a+b+c+d];
}

double fastFishersExactTest(int a, int b, int c, int d){
    int n = a + b + c + d;
    double* logFacs = new double[n+1]; // *** dynamically allocate memory logFacs[0..n] ***
    initLogFacs(logFacs, n); // *** initialize logFacs array ***
    double logpCutoff = logHypergeometricProb(logFacs,a,b,c,d); // *** logFacs added
    double pFraction = 0;
    for(int x=0; x <= n; ++x) {
        if ( a+b-x >= 0 && a+c-x >= 0 && d-a+x >=0 ) {
            double l = logHypergeometricProb(logFacs,x,a+b-x,a+c-x,d-a+x);
            if ( l <= logpCutoff ) pFraction += exp(l - logpCutoff);
        }
    }
    double logpValue = logpCutoff + log(pFraction);
    return logpValue/log(10.);
}

