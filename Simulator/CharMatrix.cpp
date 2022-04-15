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
  alpha:
  beta:
  sharingTimes: times of horizontal transfer events
 */
 
CharMatrix::CharMatrix(Tree* t, double** q, int ns, double* freqs, int nc, double alphaR, double alphaS, double betaS, double se, double delta) {
        
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    // determine sharing events
    std::vector<Node*> sourceNodes;
    t->addSharingEvents(&rng, se, sourceNodes, delta);

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
                
                resilience[c] = Probability::Beta::rv(&rng, alphaS, betaS);                
                
                if (alphaR < 50.0)
                    siteRate[c] = Probability::Gamma::rv(&rng, alphaR, alphaR);
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
     
     std::vector<double> eq;
     for (int i = 0; i < numStates; i++)
         eq.push_back(freqs[i]);
     
     // initialize instance variables
     double currTime = 0.0;
     double duration = 1.0;
     double sharingProb = 0.5;
         
     std::vector<Node*> activeNodes;
     
     while (currTime < duration) {
                 
         double rate = 0.0;
         
         if (currTime == 0.0) {
             
             activeNodes.push_back(t->getDownPassSequence()[-1]);
             CognateSet* root = new CognateSet(numChar, rng, eq);
             activeNodes.back()->setCognateSet(root);
         
         } else {
             
             activeNodes = t->nodesAtTime(currTime);
             
         }
         
         int numActiveNodes = (int)activeNodes.size();
         for (int i = 0; i < numActiveNodes; i++)
             rate += activeNodes[i]->getCognateSet()->calculateRate(q, numStates);
         std::cout << "rate: " << rate << std::endl;
         currTime += -log(rng->uniformRv())/rate;
         
         std::cout << currTime << std::endl;
         
         if (currTime < duration) {
             
             double u = rng->uniformRv()*rate;
             std::cout << "u: " << u << ", rate: " << rate <<std::endl;
             if (u <= rate) {
                                 
                 double u = rng->uniformRv();
                 double sum = 0.0;
                 CognateSet* p = NULL;
                                 
                 for (int i = 0; i < numActiveNodes; i++) {
                     
                     Node* pAncs = activeNodes[i]->getAncestor();
                     
                     if (pAncs == NULL) {
                         continue;
                     } else {
                         sum += activeNodes[i]->getCognateSet()->calculateRate(q, numStates);
                         if (u < sum) {
                             p = activeNodes[i]->getCognateSet();
                             break;
                         }
                     }
                 }
                 p->changeCognate(q, numStates, rng);
                 p->print();
             }
         }
         
         for (int i = 0; i < numActiveNodes; i++) {
             Node* p = activeNodes[i];
             if (p->getNumDescendents() != 0) {
                 for (Node* n : p->getDescendants())
                     n->setCognateSet(p->getCognateSet());
             }
         }
    
    // simulate data with sharing events
    
    std::vector<Node*>& dpseq = t->getDownPassSequence();

    for (Node* n : dpseq) {
        std::cout << "woof" << std::endl;
        if (n->getSource() == NULL && n->getDest() == NULL) {
            continue;
        } else if (n->getSource() != NULL) {
            // if it has a source, then it is a destination
            Node* source = n->getSource();
            Node* destination = n;
                                
            std::vector<int>& sourceCogs = source->getCognateSet()->getCognates();
            std::vector<int>& destCogs = destination->getCognateSet()->getCognates();
            
            for (int i = 0; i < (int)sourceCogs.size(); i++) {
                double u = rng->uniformRv();
                if (u < sharingProb)
                    destCogs[i] = sourceCogs[i];
            }
            
            // simulate data on destination subtree
            
            
        }
    }
    
    
    



 
CharMatrix::CharMatrix(Tree* t, double** q, int ns, double* freqs, int nc, RandomVariable* rng, double alpha, double beta, double se) {
    
    std::vector<Node*>& dpseq = t->getDownPassSequence();
    
    // initialize instance variables
    numChar = nc;
    numTaxa = 0;
    numStates = ns;
    for (Node* n : dpseq) {
        if (n->getDescendants().size() == 0)
            numTaxa++;
    }
    std::vector<double> eq;
    for (int i = 0; i < numStates; i++)
        eq.push_back(freqs[i]);
    
    // initialize root
    CognateSet* root = new CognateSet(numChar, rng, eq);
    Node* rootNode = dpseq[(int)dpseq.size()-1];
    std::vector<CognateSet*> activeLangs;
    std::vector<Node*> activeNodes;
    activeLangs.push_back(root);
    activeNodes.push_back(rootNode);
    double currTime = 0.0;
    double duration = 1.0;
    double sharingProb = 0.5;
    double changeRate = 0.0;
    double sharingRate = 0.0 * se;
    double rate = changeRate + sharingRate;
    std::vector<int> remove;
    
    while (currTime < duration && (int)(activeNodes.size()) > 0) {

        // iterate through active nodes
        for (int i = 0; i < (int)activeNodes.size(); i++) {
            
            // check if currTime exceeds active node, if so then add descendents as active node
            if (currTime > activeNodes[i]->getTotalBrLen()) {
               
                // if the active node is a tip or extinct lineage, then add to list of active nodes to remove
                if (activeNodes[i]->getNumDescendents() == 0) {
                    remove.push_back(i);
                    continue;
                }
                                
                std::set<Node*>& descendents = activeNodes[i]->getDescendants(); // get descendents
                CognateSet* currCognate = activeLangs[i];
                remove.push_back(i); // add current node to list of active nodes to remove
                
                // iterate through immediate descendents
                for (Node* n : descendents) {
                    CognateSet* q = new CognateSet(*currCognate); // duplicate cognate set
                    activeNodes.push_back(n); // add descendent to active nodes
                    activeLangs.push_back(q); // add cognate set to active langs
                    std::cout << "n: " << n->getIndex() << ", currTime: " << currTime << std::endl;
                }
            }
        }
        
        // active nodes and langs in remove list gets removed
        for (int i : remove) {
            activeNodes[i] = activeNodes[activeNodes.size()-1];
            activeLangs[i] = activeLangs[activeLangs.size()-1];
            activeNodes.pop_back();
            activeLangs.pop_back();
        }
        
        remove.clear(); // clear out list of active nodes to remove
        
        // update instance variables based on current active nodes
        for (int i = 0; i < (int)activeLangs.size(); i++)
            changeRate += activeLangs[i]->calculateRate(q, numStates);
        sharingRate = activeLangs.size() * se; // update sharing rate
        if ((int)activeLangs.size() == 1)
            sharingRate = 0.0;
        rate = sharingRate + changeRate; // update overall rate
        currTime += -log(rng->uniformRv())/rate;
        double u = rng->uniformRv()*rate;
        
        if (u <= sharingRate) {
            CognateSet* source = activeLangs[(int)(rng->uniformRv()*activeLangs.size())];
            CognateSet* destination = NULL;
            
            do {
                destination = activeLangs[(int)(rng->uniformRv()*activeLangs.size())];
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
            for (int i = 0; i < (int)activeLangs.size(); i++) {
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
