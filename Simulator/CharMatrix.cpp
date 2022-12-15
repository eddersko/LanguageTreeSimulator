//
//  CharMatrix.cpp
//  BirthDeathSimulation
//
//  Created by Edwin Ko on 1/26/22.
//

#include "RandomVariable.hpp"
#include "CharMatrix.hpp"
#include "Tree.hpp"
#include "Node.hpp"
#include "Probability.hpp"
#include <cmath>
#include <iostream>
#include "CognateSet.hpp"
#include "Msg.hpp"

bool compareTimes(Node* n1, Node* n2) {
    return (n1->getTime() < n2->getTime());
}

/*
 Simulate data on a given tree input with horizontal transfer.
 
  t: tree
  q: rate matrix
  ns: number of states
  freqs: equilibrium frequencies
  nc: number of charactes/cognates
  rng: random number generator
  alphaRat: alpha and beta parameter of gamma distribution for rate heterogeneity
  alphaRes: alpha parameter of beta distribution for resilience/resistance to borrowing
  betaRes: beta parameter of beta distribution for resilience/resistance to borrowing
  sharingTimes: times of horizontal transfer events
  delta: likelihood of borrowing locally or distantly
 */
 
CharMatrix::CharMatrix(Tree* t, double** q, int ns, std::vector<double> freqs, int nc, double alphaRat, double alphaRes, double betaRes, double sharingRate, double delta) {
        
    double expectedNumChanges = 2.0;
    double subtreeLength = t->rescale();
    double rateFactor = expectedNumChanges / subtreeLength;
    
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    // determine sharing events
    std::vector<Node*> sourceNodes;
    t->addSharingEvents(&rng, sharingRate/subtreeLength, sourceNodes, delta);

    // simulate data
    numStates = ns;
    numChar = nc;
    resilience = new double[nc];
    siteRate = new double[nc];
    //Probability::Gamma::discretization(rateVar, alphaR, alphaR, 4, false);
    
    std::vector<double> stateFreqs;
    for (int i = 0; i < numStates; i++)
        stateFreqs.push_back(freqs[i]);
    std::vector<Node*>& dpseq = t->getDownPassSequence();
    
    // add CognateSets
    for (Node* n : dpseq)
        n->setCognateSet(new CognateSet(numChar, &rng, stateFreqs));
        
        
    for (int n = (int)dpseq.size()-1; n >= 0; n--) {
        Node* p = dpseq[n];
                
        CognateSet* cs = p->getCognateSet();
        
        if (cs == NULL)
            Msg::error("There should be a cognate set that is not NULL!");
        
        if (p == t->getRoot()) {
            
            // no need to do anything, cognate set already constructed from stationary probabilities
            
        } else {
         
            for (int c = 0; c < numChar; c++) {
                
                int currState = (*p->getAncestor()->getCognateSet())[c];
                double len = (p->getTime() - p->getAncestor()->getTime()) * rateFactor;
                
                resilience[c] = Probability::Beta::rv(&rng, alphaRes, betaRes);
                
                if (alphaRat < 50.0)
                    siteRate[c] = Probability::Gamma::rv(&rng, alphaRat, alphaRat);
                else
                    siteRate[c] = 1.0;
                                
                double v = 0.0;
                double cLen = len * siteRate[c];
                
                while (v < cLen) {
                                        
                    double rate = -q[currState][currState];
                            
                    v += -log(rng.uniformRv())/rate;
                    
                    if (v < cLen) {
                        
                        double u = rng.uniformRv();
                        double sum = 0.0;
                        
                        for (int i = 0; i < numStates; i++) {
                            
                            sum += q[currState][i] / rate;
                            //std::cout << u << " " << sum << std::endl;
                            if (u < sum) {
                                currState = i;
                                break;
                            }
                        }
                    }
                }
                (*cs)[c] = currState;
            }
                                    
        }
    }

    // collect all sourceNodes sorted by times
    sort(sourceNodes.begin(), sourceNodes.end(), compareTimes);
    
    //for (int i = 0; i < (int)sourceNodes.size(); i++)
        //std::cout << sourceNodes[i]->getTime() << std::endl;
        
    //t->print();
    
    // simulate in order sharing then resimulate history from destination node
    for (int i = 0; i < (int)sourceNodes.size(); i++) {
        
        Node* dest = sourceNodes[i]->getDest();
                
        
        if (dest == NULL) {
            continue;
            //Msg::error("dest should not be NULL.");
        }
        
        for (int j = 0; j < numChar; j++) {
            double randomNumber = rng.uniformRv();
            if (randomNumber > resilience[j]) {
                //std::cout << "Sharing between " << dest->getIndex() << " and " << sourceNodes[i]->getIndex() << " in site " << std::endl;
                (*dest->getCognateSet())[j] = (*sourceNodes[i]->getCognateSet())[j];
                simulateSubTree(dest, dest, q, &rng, j);
            }
        }
        
    }
    
    // count number of taxa on tree
    
    numTaxa = 0;    
    for (Node* n : t->getDownPassSequence()) {
        if (n->getIsTip() == true)
            numTaxa++;
        
    }
            
    
    // dynamically allocate the matrix
    matrix = new int*[numTaxa];
    matrix[0] = new int[numTaxa*numChar];
    for (int i = 1; i < numTaxa; i++)
        matrix[i] = matrix[i-1]+numChar;
    for (int i = 0; i < numTaxa; i++)
        for (int j = 0; j < numChar; j++)
            matrix[i][j] = 0;
    
    // put tip cognate sets into matrix
    for (Node* n : t->getDownPassSequence()) {
        if (n->getIsTip() == true) {
            int tipIdx = n->getIndex();
            if (tipIdx >= numTaxa)
                Msg::error("Tip index is too large!");
            CognateSet* cs = n->getCognateSet();
            
            for (int i = 0; i < numChar; i++)
                matrix[tipIdx][i] = (*cs)[i];
        }
    }
    
    
    delete resilience;
    delete siteRate;
}

/*
 Simulate data on a given tree input with internal and external borrowing.
 
  t1: desination tree
  t2: source tree
  q: rate matrix
  ns: number of states
  freqs: equilibrium frequencies
  nc: number of charactes/cognates
  rng: random number generator
  alphaRat: alpha and beta parameter of gamma distribution for rate heterogeneity
  alphaRes: alpha parameter of beta distribution for resilience/resistance to borrowing
  betaRes: beta parameter of beta distribution for resilience/resistance to borrowing
  sharingRate: rates of borrowing events
  ratio: ratio of external borrowing to internal borrowing
  delta: likelihood of borrowing locally or distantly
 */
    
CharMatrix::CharMatrix(Tree* destTree, Tree* sourceTree, double** q, int ns, std::vector<double> freqs, int nc, double alphaRat, double alphaRes, double betaRes, double sharingRate, double ratio, double delta) {
            
    double expectedNumChanges = 2.0;
    double subtreeLength = destTree->rescale();
    double rateFactor = expectedNumChanges/subtreeLength;
    
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    double exSharingRate = sharingRate * ratio;
    double inSharingRate = sharingRate - (sharingRate * ratio);
        
    //std::cout << ratio << std::endl;
    //std::cout << sharingRate << std::endl;
    // determine sharing events
    std::vector<Node*> internalSourceNodes;
    std::vector<Node*> externalSourceNodes;
    
    //std::cout << "subtreeLength: " << subtreeLength << ", inSharingRate: " << (inSharingRate / subtreeLength) << ", delta: " << delta << std::endl;
    
    //std::cout << "sourceNodes size (before): " << sourceNodes.size() << std::endl;
    destTree->addSharingEvents(&rng, inSharingRate / subtreeLength, internalSourceNodes, delta); // internal borrowing events
    //std::cout << "sourceNodes size (after internal): " << sourceNodes.size() << std::endl;
    destTree->addSharingEvents(&rng, sourceTree, exSharingRate / subtreeLength, externalSourceNodes); // external borrowing events
    //std::cout << "sourceNodes size (after external): " << sourceNodes.size() << std::endl;

    // simulate data
    numStates = ns;
    numChar = nc;
    resilience = new double[nc];
    siteRate = new double[nc];
    
    std::vector<Node*>& sourceDpseq = sourceTree->getDownPassSequence();
    std::vector<Node*>& destDpseq = destTree->getDownPassSequence();
    
    std::vector<double> stateFreqs;
    for (int i = 0; i < numStates; i++)
        stateFreqs.push_back(freqs[i]);
                
    // add CognateSets
    for (Node* n : sourceDpseq)
        n->setCognateSet(new CognateSet(numChar, &rng, stateFreqs));
        
    for (Node* n : destDpseq)
        n->setCognateSet(new CognateSet(numChar, &rng, stateFreqs)); // MEMORY LEAK
    
    for (int n = (int)sourceDpseq.size()-1; n >= 0; n--) {
        Node* p = sourceDpseq[n];
                
        CognateSet* cs = p->getCognateSet();
        
        if (cs == NULL)
            Msg::error("There should be a cognate set that is not NULL!");
        
        if (p == sourceTree->getRoot()) {
            
            // no need to do anything, cognate set already constructed from stationary probabilities
            
        } else {
         
            for (int c = 0; c < numChar; c++) {
                
                int currState = (*p->getAncestor()->getCognateSet())[c];
                
                
                
                double len = (p->getTime() - p->getAncestor()->getTime()) * rateFactor;
                
                //std::cout << "raw length: " << (p->getTime() - p->getAncestor()->getTime()) << ", adjusted length: " << len << std::endl;
                
                double v = 0.0;
                
                while (v < len) {
                                        
                    double rate = -q[currState][currState];
                            
                    v += -log(rng.uniformRv())/rate;
                    
                    if (v < len) {
                        
                        double u = rng.uniformRv();
                        double sum = 0.0;
                        
                        for (int i = 0; i < numStates; i++) {
                            
                            sum += q[currState][i] / rate;
                            //std::cout << u << " " << sum << std::endl;
                            if (u < sum) {
                                currState = i;
                                break;
                            }
                        }
                    }
                }
                (*cs)[c] = currState;
            }
        }
    }
    
        
    for (int n = (int)destDpseq.size()-1; n >= 0; n--) {
        Node* p = destDpseq[n];
                
        CognateSet* cs = p->getCognateSet();
        
        if (cs == NULL)
            Msg::error("There should be a cognate set that is not NULL!");
        
        if (p == destTree->getRoot()) {
            
            // no need to do anything, cognate set already constructed from stationary probabilities
            
        } else {
         
            for (int c = 0; c < numChar; c++) {
                
                int currState = (*p->getAncestor()->getCognateSet())[c];
                double len = (p->getTime() - p->getAncestor()->getTime()) * rateFactor;
                
                resilience[c] = Probability::Beta::rv(&rng, alphaRes, betaRes);
                
                if (alphaRat < 50.0)
                    siteRate[c] = Probability::Gamma::rv(&rng, alphaRat, alphaRat);
                else
                    siteRate[c] = 1.0;
                                
                double v = 0.0;
                double cLen = len * siteRate[c];
                
                while (v < cLen) {
                                        
                    double rate = -q[currState][currState];
                            
                    v += -log(rng.uniformRv())/rate;
                    
                    if (v < cLen) {
                        
                        double u = rng.uniformRv();
                        double sum = 0.0;
                        
                        for (int i = 0; i < numStates; i++) {
                            
                            sum += q[currState][i] / rate;
                            //std::cout << u << " " << sum << std::endl;
                            if (u < sum) {
                                currState = i;
                                break;
                            }
                        }
                    }
                }
            (*cs)[c] = currState;
            }
                                    
        }
    }

    std::vector<Node*> sourceNodes;
    
    for (Node* n : internalSourceNodes)
        sourceNodes.push_back(n);
    for (Node* n : externalSourceNodes)
        sourceNodes.push_back(n);
        
    // collect all sourceNodes sorted by times
    sort(sourceNodes.begin(), sourceNodes.end(), compareTimes);
    
    //for (int i = 0; i < (int)sourceNodes.size(); i++)
        //std::cout << sourceNodes[i]->getTime() << std::endl;
        
    //t->print();
    
    // simulate in order sharing then resimulate history from destination node
    for (int i = 0; i < (int)sourceNodes.size(); i++) {
        
        Node* dest = sourceNodes[i]->getDest();
        
        if (dest == NULL) {
            continue;
            //Msg::error("dest should not be NULL.");
        }
        
        for (int j = 0; j < numChar; j++) {
            double randomNumber = rng.uniformRv();
            if (randomNumber > resilience[j]) {
                //std::cout << "Sharing between " << dest->getIndex() << " (" << (*dest->getCognateSet())[j] << ") and " << sourceNodes[i]->getIndex() << " (" << (*sourceNodes[i]->getCognateSet())[j] << ") in site " << std::endl;
                                                                                                
                (*dest->getCognateSet())[j] = (*sourceNodes[i]->getCognateSet())[j];
                simulateSubTree(dest, dest, q, &rng, j);
            }
        }
        
    }
    
    
    // count number of taxa on tree
    
    numTaxa = 0;
    for (Node* n : destTree->getDownPassSequence()) {
        if (n->getIsTip() == true)
            numTaxa++;
    }
    
    // dynamically allocate the matrix
    matrix = new int*[numTaxa];
    matrix[0] = new int[numTaxa*numChar];
    for (int i = 1; i < numTaxa; i++)
        matrix[i] = matrix[i-1]+numChar;
    for (int i = 0; i < numTaxa; i++)
        for (int j = 0; j < numChar; j++)
            matrix[i][j] = 0;
    
    // put tip cognate sets into matrix
    for (Node* n : destTree->getDownPassSequence()) {
        if (n->getIsTip() == true) {
            int tipIdx = n->getIndex();
            if (tipIdx >= numTaxa)
                Msg::error("Tip index is too large!");
            CognateSet* cs = n->getCognateSet();
            
            for (int i = 0; i < numChar; i++)
                matrix[tipIdx][i] = (*cs)[i];
        }
    }
    
       
   
    delete resilience;
    delete siteRate;
    
}

/* This CharMatrix simulates both internal and external borrowing. */

CharMatrix::CharMatrix(Tree* t, double** q, int ns, std::vector<double> freqs, int nc, double alphaRat, double alphaRes, double betaRes, double sharingRate, double ratio, double delta) {
    
    
    double expectedNumChanges = 2.0;
    double subtreeLength = t->rescale();
    double rateFactor = expectedNumChanges / subtreeLength;
    double externalSharingRate = sharingRate * ratio;
    double internalSharingRate = sharingRate - (sharingRate * ratio);
    
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    // determine sharing events
    std::vector<Node*> sourceNodes;
    t->addSharingEvents(&rng, internalSharingRate/subtreeLength, sourceNodes, delta);
    t->addExternalSharingEvents(&rng, externalSharingRate/subtreeLength, sourceNodes);
    
    // simulate data
    numStates = ns;
    numChar = nc;
    resilience = new double[nc];
    siteRate = new double[nc];
    //Probability::Gamma::discretization(rateVar, alphaR, alphaR, 4, false);

    std::vector<double> stateFreqs;
    for (int i = 0; i < numStates; i++)
        stateFreqs.push_back(freqs[i]);
    std::vector<Node*>& dpseq = t->getDownPassSequence();

    // add CognateSets
    for (Node* n : dpseq)
        n->setCognateSet(new CognateSet(numChar, &rng, stateFreqs));
        
    
    // simulating evolution along the tree
    for (int n = (int)dpseq.size()-1; n >= 0; n--) {
        Node* p = dpseq[n];
                
        CognateSet* cs = p->getCognateSet();
        
        if (cs == NULL)
            Msg::error("There should be a cognate set that is not NULL!");
        
        if (p == t->getRoot()) {
            
            // no need to do anything, cognate set already constructed from stationary probabilities
            
        } else {
         
            for (int c = 0; c < numChar; c++) {
                
                int currState = (*p->getAncestor()->getCognateSet())[c];
                double len = (p->getTime() - p->getAncestor()->getTime()) * rateFactor;
                
                resilience[c] = Probability::Beta::rv(&rng, alphaRes, betaRes);
                
                if (alphaRat < 50.0)
                    siteRate[c] = Probability::Gamma::rv(&rng, alphaRat, alphaRat);
                else
                    siteRate[c] = 1.0;
                                
                double v = 0.0;
                double cLen = len * siteRate[c];
                
                while (v < cLen) {
                                        
                    double rate = -q[currState][currState];
                            
                    v += -log(rng.uniformRv())/rate;
                    
                    if (v < cLen) {
                        
                        double u = rng.uniformRv();
                        double sum = 0.0;
                        
                        for (int i = 0; i < numStates; i++) {
                            
                            sum += q[currState][i] / rate;
                            //std::cout << u << " " << sum << std::endl;
                            if (u < sum) {
                                currState = i;
                                break;
                            }
                        }
                    }
                }
                (*cs)[c] = currState;
            }
                                    
        }
    }

    
    // this is where the fun begins!
    
    // collect all sourceNodes sorted by times
    sort(sourceNodes.begin(), sourceNodes.end(), compareTimes);

    //for (int i = 0; i < (int)sourceNodes.size(); i++)
        //std::cout << sourceNodes[i]->getTime() << std::endl;
        
    //t->print();

    // simulate in order sharing then resimulate history from destination node
    for (int i = 0; i < (int)sourceNodes.size(); i++) {
        
        Node* dest = sourceNodes[i]->getDest();
        
        if (dest == NULL)
            continue;
            //Msg::error("dest should not be NULL.");
        
                        
        for (int j = 0; j < numChar; j++) {
            double randomNumber = rng.uniformRv();
            if (randomNumber > resilience[j]) {
                //std::cout << "Sharing between " << dest->getIndex() << " and " << sourceNodes[i]->getIndex() << " in site " << std::endl;
                if (dest->getIsExternal() == true) {
                                        
                    // get the node's cognate set, add new site and set value to 1
                    
                    if (dest->getCognateSet() == NULL)
                        Msg::error("Cognate set should not be NULL.");
                                        
                    dest->getCognateSet()->incrementNumSites(1);
                    addSitesToDescendants(dest, 1);
                    //CognateSet* cs = dest->getCognateSet();
                    //dest->setCognateSet(cs);
                    //delete cs;
                    
                    int siteNum = dest->getCognateSet()->getNumCognates() - 1;
                    
                    // get cognate sets of other nodes at existing at that time, add new site, and set value to ?
                    
                    std::vector<Node*> activeNodes = t->nodesAtTime(dest->getTime());
                    
                    // get descendants of activeNodes and increment number of sites by one
                    
                    for (int k = 0; k < activeNodes.size(); k++) {
                        activeNodes[k]->getCognateSet()->incrementNumSites(-1);
                        addSitesToDescendants(activeNodes[k], -1);
                        
                        //CognateSet* cs = activeNodes[k]->getCognateSet();
                        //activeNodes[k]->setCognateSet(cs);
                        //delete cs;
                        
                    }
                                        
                    simulateSubTree(dest, dest, q, &rng, siteNum);
                                        
                } else {
                    (*dest->getCognateSet())[j] = (*sourceNodes[i]->getCognateSet())[j];
                    simulateSubTree(dest, dest, q, &rng, j);
                }
            }
        }
        
    }

    // count number of taxa on tree

    numTaxa = 0;
    for (Node* n : t->getDownPassSequence()) {
        if (n->getIsTip() == true)
            numTaxa++;
        
    }
            

    // dynamically allocate the matrix
    matrix = new int*[numTaxa];
    matrix[0] = new int[numTaxa*numChar];
    for (int i = 1; i < numTaxa; i++)
        matrix[i] = matrix[i-1]+numChar;
    for (int i = 0; i < numTaxa; i++)
        for (int j = 0; j < numChar; j++)
            matrix[i][j] = 0;

    // put tip cognate sets into matrix
    for (Node* n : t->getDownPassSequence()) {
        if (n->getIsTip() == true) {
            int tipIdx = n->getIndex();
            if (tipIdx >= numTaxa)
                Msg::error("Tip index is too large!");
            CognateSet* cs = n->getCognateSet();
            
            for (int i = 0; i < numChar; i++)
                matrix[tipIdx][i] = (*cs)[i];
        }
    }


    delete resilience;
    delete siteRate;

}


/*
 Simulate data on a given tree input with horizontal transfer and temporal biasing.
 
  t: tree
  q: rate matrix
  ns: number of states
  freqs: equilibrium frequencies
  nc: number of charactes/cognates
  rng: random number generator
  alphaRat: alpha and beta parameter of gamma distribution for rate heterogeneity
  alphaRes: alpha parameter of beta distribution for resilience/resistance to borrowing
  betaRes: beta parameter of beta distribution for resilience/resistance to borrowing
  sharingTimes: times of horizontal transfer events
  delta: likelihood of borrowing locally or distantly
  episolon: likelihood of borrowing at deeper or shallower time depths
 */
 
CharMatrix::CharMatrix(Tree* t, double** q, int ns, std::vector<double> freqs, int nc, double alphaRat, double alphaRes, double betaRes, double sharingRate, double delta, bool borrowNearTips) {
        
    double expectedNumChanges = 2.0;
    double subtreeLength = t->rescale();
    double rateFactor = expectedNumChanges / subtreeLength;
    
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    // determine sharing events
    std::vector<Node*> sourceNodes;
    t->addSharingEvents(&rng, sharingRate / subtreeLength, sourceNodes, delta, borrowNearTips);

    // simulate data
    numStates = ns;
    numChar = nc;
    resilience = new double[nc];
    siteRate = new double[nc];
    //Probability::Gamma::discretization(rateVar, alphaR, alphaR, 4, false);
    
    std::vector<double> stateFreqs;
    for (int i = 0; i < numStates; i++)
        stateFreqs.push_back(freqs[i]);
    std::vector<Node*>& dpseq = t->getDownPassSequence();
    
    // add CognateSets
    for (Node* n : dpseq)
        n->setCognateSet(new CognateSet(numChar, &rng, stateFreqs));
        
        
    for (int n = (int)dpseq.size()-1; n >= 0; n--) {
        Node* p = dpseq[n];
                
        CognateSet* cs = p->getCognateSet();
        
        if (cs == NULL)
            Msg::error("There should be a cognate set that is not NULL!");
        
        if (p == t->getRoot()) {
            
            // no need to do anything, cognate set already constructed from stationary probabilities
            
        } else {
         
            for (int c = 0; c < numChar; c++) {
                                                
                int currState = (*p->getAncestor()->getCognateSet())[c];
                double len = p->getTime() - p->getAncestor()->getTime();
                
                resilience[c] = Probability::Beta::rv(&rng, alphaRes, betaRes);
                
                if (alphaRat < 50.0)
                    siteRate[c] = Probability::Gamma::rv(&rng, alphaRat, alphaRat);
                else
                    siteRate[c] = 1.0;
                                
                double v = 0.0;
                double cLen = len * siteRate[c];
                
                while (v < cLen) {
                                        
                    double rate = -q[currState][currState];
                            
                    v += -log(rng.uniformRv())/rate;
                    
                    if (v < cLen) {
                        
                        double u = rng.uniformRv();
                        double sum = 0.0;
                        
                        for (int i = 0; i < numStates; i++) {
                            
                            sum += q[currState][i] / rate;
                            //std::cout << u << " " << sum << std::endl;
                            if (u < sum) {
                                currState = i;
                                break;
                            }
                        }
                    }
                }
                (*cs)[c] = currState;
            }
                                    
        }
    }

    // collect all sourceNodes sorted by times
    sort(sourceNodes.begin(), sourceNodes.end(), compareTimes);
    
    //for (int i = 0; i < (int)sourceNodes.size(); i++)
        //std::cout << sourceNodes[i]->getTime() << std::endl;
        
    //t->print();
    
    // simulate in order sharing then resimulate history from destination node
    for (int i = 0; i < (int)sourceNodes.size(); i++) {
        
        Node* dest = sourceNodes[i]->getDest();
                
        
        if (dest == NULL) {
            continue;
            //Msg::error("dest should not be NULL.");
        }
        
        for (int j = 0; j < numChar; j++) {
            double randomNumber = rng.uniformRv();
            if (randomNumber > resilience[j]) {
                //std::cout << "Sharing between " << dest->getIndex() << " and " << sourceNodes[i]->getIndex() << " in site " << std::endl;
                (*dest->getCognateSet())[j] = (*sourceNodes[i]->getCognateSet())[j];
                simulateSubTree(dest, dest, q, &rng, j);
            }
        }
        
    }
    
    // count number of taxa on tree
    
    numTaxa = 0;
    for (Node* n : t->getDownPassSequence()) {
        if (n->getIsTip() == true)
            numTaxa++;
        
    }
            
    
    // dynamically allocate the matrix
    matrix = new int*[numTaxa];
    matrix[0] = new int[numTaxa*numChar];
    for (int i = 1; i < numTaxa; i++)
        matrix[i] = matrix[i-1]+numChar;
    for (int i = 0; i < numTaxa; i++)
        for (int j = 0; j < numChar; j++)
            matrix[i][j] = 0;
    
    // put tip cognate sets into matrix
    for (Node* n : t->getDownPassSequence()) {
        if (n->getIsTip() == true) {
            int tipIdx = n->getIndex();
            if (tipIdx >= numTaxa)
                Msg::error("Tip index is too large!");
            CognateSet* cs = n->getCognateSet();
            
            for (int i = 0; i < numChar; i++)
                matrix[tipIdx][i] = (*cs)[i];
        }
    }
    
    
    delete resilience;
    delete siteRate;
}

void CharMatrix::addSitesToDescendants(Node* n, int val) {
    
    std::set<Node*>& nDes = n->getDescendants();
    for (Node* p : nDes) {
        p->getCognateSet()->incrementNumSites(val);
    }
}


void CharMatrix::simulateSubTree(Node* n, Node* r, double** q, RandomVariable* rng, int site) {
    
    if (n != NULL) {
        
        // simulate here
        if (n != r) {
            Node* nAncs = n->getAncestor();
            if (nAncs == NULL)
                Msg::error("We have a problem -- nAncs should not be NULL!");
            double brLen = n->getTime()-nAncs->getTime();
            
            int currState = (*nAncs->getCognateSet())[site]; // get current state
            double v = 0.0;
            while (v < brLen) {
                double rate = -q[currState][currState];
                v += -log(rng->uniformRv()) / rate; // increment time
                if (v < brLen) {
                    double u = rng->uniformRv();
                    double sum = 0.0;
                    for (int j = 0; j < numStates; j++) {
                        if (j != currState) {
                            sum += q[currState][j] / rate;
                            if (u < sum) {
                                currState = j;
                                break;
                            }
                        }
                    }
                }
                (*n->getCognateSet())[site] = currState;
                //if (n->getDescendants().size() == 0)
                  //  matrix[n->getIndex()][site] = currState;
            }
        }
        std::set<Node*>& nDes = n->getDescendants();
        for (Node* p : nDes)
            simulateSubTree(p, r, q, rng, site);
        
    }
    
    
}
   
    

/*
 
 Simulate tree and data with horizontal transfer.
 
  t: tree
  q: rate matrix
  ns: number of states
  freqs: equilibrium frequencies
  nc: number of charactes/cognates
  rng: random number generator
  alpha:
  beta:
  sharingTimes: times of horizontal transfer events
 */
CharMatrix::CharMatrix(double** q, int ns, double* freqs, int nc, RandomVariable* rng, double alpha, double beta, double se) {
    
    numStates = ns;
    numChar = nc;
    std::vector<double> eq;
    for (int i = 0; i < numStates; i++)
        eq.push_back(freqs[i]);
    
    CognateSet* root = new CognateSet(numChar, rng, eq);
    std::vector<CognateSet*> activeLangs;
    activeLangs.push_back(root);
    
    double currTime = 0.0;
    double duration = 1.0;
    double lambda = 5.0;
    double mu = 1.0;
    double sharingProb = 0.5;
    
    while (currTime < duration && (int)activeLangs.size() > 0) {
            
        int numActiveLangs = (int)activeLangs.size();
        double speciationRate = numActiveLangs * lambda;
        double extinctionRate = numActiveLangs * mu;
        double changeRate = 0.0;
        for (int i = 0; i < numActiveLangs; i++)
            changeRate += activeLangs[i]->calculateRate(q, numStates);
        double sharingRate = numActiveLangs * se;
        if (numActiveLangs == 1)
            sharingRate = 0.0;
        double rate = speciationRate + extinctionRate + changeRate + sharingRate;
        
        currTime += -log(rng->uniformRv())/rate;
        if (currTime < duration) {
            double u = rng->uniformRv()*rate;
            if (u <= speciationRate) {
                
                // speciation
                int whichCognateSet = (int)(rng->uniformRv()*numActiveLangs);
                CognateSet* p = activeLangs[whichCognateSet];
                CognateSet* q = new CognateSet(*p);
                CognateSet* r = new CognateSet(*p);
                p->addDescendant(q);
                p->addDescendant(r);
                q->setAncestor(p);
                r->setAncestor(p);
                
                activeLangs[whichCognateSet] = q;
                activeLangs.push_back(r);
                
            } else if (u > speciationRate && u <= (speciationRate+extinctionRate)) {
                
                // extinction
                int whichCognateSet = (int)(rng->uniformRv()*numActiveLangs);
                //CognateSet* p = activeLangs[whichCognateSet];
                
                activeLangs[whichCognateSet] = activeLangs[numActiveLangs-1];
                activeLangs.pop_back();
                
            } else if (u > (speciationRate+extinctionRate) && u <= speciationRate + extinctionRate + sharingRate) {
                
                // share
                CognateSet* source = activeLangs[(int)rng->uniformRv()*numActiveLangs];
                CognateSet* destination = NULL;
                
                do {
                    destination = activeLangs[(int)(rng->uniformRv()*numActiveLangs)];
                } while (source == destination);
                
                std::vector<int>& sourceCogs = source->getCognates();
                std::vector<int>& destCogs = destination->getCognates();
                
                for (int i = 0; i < (int)sourceCogs.size(); i++) {
                    double u = rng->uniformRv();
                    if (u < sharingProb)
                        destCogs[i] = sourceCogs[i];
                }
            } else {
                
                double u = rng->uniformRv();
                double sum = 0.0;
                CognateSet* p = NULL;
                
                // choose lineage
                for (int i = 0; i < numActiveLangs; i++) {
                    sum += activeLangs[i]->calculateRate(q, numStates);
                    if (u < sum) {
                        p = activeLangs[i];
                        break;
                    }
                }
                p->changeCognate(q, numStates, rng);
            }
        }
    }
}


/*
    Simulate data on a given tree input.
    t: tree object
    q: rate matrix
    freqs: stationary frequencies
    nc: number of characters/cognates
    rng: random number generator
 */
CharMatrix::CharMatrix(Tree* t, double** q, int ns, double* freqs, int nc, RandomVariable* rng) {
    
    std::vector<Node*>& dpseq = t->getDownPassSequence();

    // initialize instance variables
    numChar = nc;
    numTaxa = 0;
    numStates = ns;
    for (Node* n : dpseq) {
        if (n->getDescendants().size() == 0)
            numTaxa++;
    }
    
    // dynamically allocate matrix
    matrix = new int*[numTaxa];
    matrix[0] = new int[numTaxa*numChar];
    for (int i = 1; i < numTaxa; i++)
        matrix[i] = matrix[i-1]+numChar;
    for (int i = 0; i < numTaxa; i++)
        for (int j = 0; j < numChar; j++)
            matrix[i][j] = 0;
    
    // simulate the data!
    for (int c = 0; c < numChar; c++) {
        for (int n = (int)dpseq.size()-1; n >= 0; n--) { //start at root
            Node* p = dpseq[n];
            Node* pAncs = p->getAncestor();
            if (pAncs == NULL) { // if root then initualize values using equilibrium frequencies
                double u = rng->uniformRv();
                double sum = 0.0;
                for (int i = 0; i < numStates; i++) {
                    sum += freqs[i];
                    if (u < sum) {
                        p->setState(i);
                        break;
                    }
                }
            } else { // otherwise if not root
                int currState = pAncs->getState(); // get current state
                double brLen = p->getBrLen(); // get br length
                double v = 0.0;
                while (v < brLen) {
                    double rate = -q[currState][currState];
                    v += -log(rng->uniformRv()) / rate; // increment time
                    if (v < brLen) {
                        double u = rng->uniformRv();
                        double sum = 0.0;
                        for (int j = 0; j < numStates; j++) {
                            if (j != currState) {
                                sum += q[currState][j] / rate;
                                if (u < sum) {
                                    currState = j;
                                    break;
                                }
                            }
                        }
                    }
                    p->setState(currState);
                    if (p->getDescendants().size() == 0)
                        matrix[p->getIndex()][c] = currState;
                }
            }
        }
    
    }
}


CharMatrix::~CharMatrix(void) {
    delete [] matrix[0];
    delete [] matrix;
}

void CharMatrix::print(void) {
    
    for (int i = 0; i < numTaxa; i++) {
        std::cout << i << " ";
        for (int j = 0; j < numChar; j++)
            std::cout << matrix[i][j];
        std::cout << std::endl;
    }
}


std::string CharMatrix::getString(void) {
    
    std::string str = "";
    
    for (int i = 0; i < numTaxa; i++) {
        str += std::to_string(i) + " ";
        for (int j = 0; j < numChar; j++)
            str += std::to_string(matrix[i][j]);
        str += "\n";
    }
    
    return str;
}
