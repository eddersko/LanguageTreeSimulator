//
//  CognateSet.cpp
//  BirthDeathSimulation
//
//  Created by Edwin Ko on 2/4/22.
//

#include "CognateSet.hpp"
#include "RandomVariable.hpp"
#include <iostream>

/*
 Initialize values for the cognate set using equilibrian frequencies
 nc: number of characters/cognates
 rng: randon number generator
 eq: equilibrium frequencies
 */
CognateSet::CognateSet(int nc, RandomVariable* rng, std::vector<double> eq) {
        
    numCognates = nc;
    cognates.resize(numCognates);
    
    for (int i = 0; i < numCognates; i++) {
        double u = rng->uniformRv();
        double sum = 0.0;
        for (int j = 0; j < (int)eq.size(); j++) {
            sum += eq[j];
            if (u < sum) {
                cognates[i] = j;
                break;
            }
        }
    }
}

// Deep copy of cognate set
CognateSet::CognateSet(const CognateSet& c) {
    numCognates = c.numCognates;
    cognates.resize(numCognates);
    for (int i = 0; i < numCognates; i++)
        cognates[i] = c.cognates[i];
}


// rate at which the cognate changes
double CognateSet::calculateRate(double** q, int numStates) {
    
    double lambda = 0.0;
    
    for (int i = 0; i < numCognates; i++) {
        int s = cognates[i];
        lambda += -q[s][s]; 
    }
    
    return lambda;
    
}


void CognateSet::changeCognate(double** q, int numStates, RandomVariable* rng) {
    
    double rate = calculateRate(q, numStates);
    double u = rng->uniformRv()*rate;
    double sum = 0.0;
    int cognateIdx = 0;
    
    for (int i = 0; i < numCognates; i++) {
        int st = cognates[i];
        sum += -q[st][st];
        
        if (u < sum) {
            cognateIdx = i;
            break;
        }
    }
    
    u = rng->uniformRv();
    sum = 0.0;
    int currState = cognates[cognateIdx];
    
    for (int i = 0; i < numStates; i++) {
        if (i != currState) {
            sum += -q[currState][i] / q[currState][currState];
            if (u < sum) {
                cognates[cognateIdx] = i;
                break;
            }
        }
        
    }
    
}

void CognateSet::print(void) {
    listCognates();
}

void CognateSet::listCognates(void) {
    for (int i = 0; i < (int)cognates.size(); i++)
        std::cout << cognates[i];
    std::cout << std::endl;
}
