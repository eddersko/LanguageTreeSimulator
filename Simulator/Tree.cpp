//
//  Tree.cpp
//  BirthDeathSimulation
//
//  Created by Edwin Ko on 1/14/22.
//

#include "Tree.hpp"
#include "RandomVariable.hpp"
#include "Node.hpp"
#include <cmath>
#include <iostream>
#include "Msg.hpp"
#include "CognateSet.hpp"
#include <sstream>
#include <climits>


// Copy constructor

Tree::Tree(const Tree& t) {
            
    this->time = t.time;
    this->numTaxa = t.numTaxa;
    
    for (int i = 0; i < (int)t.nodes.size(); i++)
        addNode();

    this->root = nodes[t.root->getOffset()];
    
    for (int i = 0; i < (int)t.nodes.size(); i++) {
        
        Node* n = this->nodes[i];
        Node* tN = t.nodes[i];
        
        n->setName(tN->getName());
        n->setIndex(tN->getIndex());
        n->setTime(tN->getTime());
        n->setBrLen(tN->getBrLen());
        n->setIsTip(tN->getIsTip());
        n->setNewick(tN->getNewick());
        n->setRescale(tN->getRescale());
        n->setIsDistance(tN->getIsDistance());
        n->setState(tN->getState());
        
        std::set<Node*>& tDescendants = tN->getDescendants();
        
        for (Node* p : tDescendants) {
            n->addDescendant(nodes[p->getOffset()]);
        }
                
        if (tN->getAncestor() != NULL)
            n->setAncestor(nodes[tN->getAncestor()->getOffset()]);
        
        //else
          //  n->setAncestor(NULL);
        //Msg::error("tN->getAncestor() is NULL");
        
        if (tN->getDest() == NULL)
            n->setDest(NULL);
        else
            n->setDest(nodes[tN->getDest()->getOffset()]);
        
        if (tN->getSource() == NULL)
            n->setSource(NULL);
        else
            n->setSource(nodes[tN->getSource()->getOffset()]);
        
        // don't copy cognate set of copied tree
        n->setCognateSet(NULL);
                        
        
    }
                
    for (Node* n : t.downPassSequence) {
        downPassSequence.push_back(nodes[n->getOffset()]);
    }
    
    
}

/* Horizontal transfer events wih a simulated tree under a birth-death process
 
   lambda: birth rate
   mu: death rate
   t: duration/time
   sharingRate: rate of horizontal transfer
   alpha:
   beta:
   q: rate matrix
   numStates: number of states */
Tree::Tree(double lambda, double mu, double t, double sharingRate, double alpha, double beta, double** q, int numStates) {
    
    RandomVariable& rv = RandomVariable::randomVariableInstance();
    std::vector<Node*> activeLineages;
    std::vector<double> eq(2, 0.5);
    root = addNode();
    
    CognateSet* cs = new CognateSet(100, &rv, eq);
    root->setCognateSet(cs);
}

// Building tree object based on an existing tree topology
Tree::Tree(std::string newickString) {
    
    std::vector<std::string> tokens = parseNewickString(newickString);
    Node* p = NULL;
    bool readingBrlen = false;
    int idx = 0;
    
    for (int i = 0; i < tokens.size(); i++) {
        
        std::string token = tokens[i];
        
        if (token == "(") {
            
            Node* newNode = addNode();
            if (p == NULL) {
                root = newNode;
            } else {
                newNode->setAncestor(p);
                p->addDescendant(newNode);
            }
            p = newNode;
            
        } else if (token == ")" || token == ",") {
            
            if (p->getAncestor() == NULL)
                Msg::error("No ancestor found!");
            p = p->getAncestor();
            
        } else if (token == ":") {
            
            readingBrlen = true;
            
        } else if (token == ";") {
            
            if (p != root)
                Msg::error("Expecting p at root!");
            
        } else {
            if (readingBrlen == false) {
                
                Node* newNode = addNode();
                newNode->setAncestor(p);
                p->addDescendant(newNode);
                newNode->setName(token);
                newNode->setIndex(idx++);
                p = newNode;
                
            } else {
                double v = atof(token.c_str());
                p->setBrLen(v);
            }
            readingBrlen = false;
        }
    }
    initializeDownPassSequence();
    
    for (Node* n : downPassSequence) {
        if (n->getNumDescendents() > 0)
            n->setIndex(idx++);
    }
    //print();
}

/* Simulating tree under birth-death process
    lambda: birth rate
    mu: death rate
    t: duration/time  */
Tree::Tree(double lambda, double mu, double t) {
    
    RandomVariable& rv = RandomVariable::randomVariableInstance();
    root = addNode(-1);
    time = t;
    Node* p = addNode(-1);
    Node* r = addNode(-1);
    root->setTime(0.0);
    root->addDescendant(p);
    root->addDescendant(r);
    p->setAncestor(root);
    r->setAncestor(root);
    std::vector<Node*> activeNodes;
    activeNodes.push_back(p);
    activeNodes.push_back(r);
    
    
    double currentTime = 0.0;
    while (currentTime < time && activeNodes.size() > 0) {
        
        double rate = activeNodes.size() * (lambda + mu);
        currentTime += -log(rv.uniformRv()) / rate;
        
        if (currentTime < time) {
            
            int whichNode = (int)(rv.uniformRv()*activeNodes.size());
            p = activeNodes[whichNode];
            p->setTime(currentTime);
            if (rv.uniformRv() < lambda / (lambda+mu)) { // speciation
                
                Node* q = addNode(-1);
                r = addNode(-1);
                q->setAncestor(p);
                r->setAncestor(p);
                p->addDescendant(q);
                p->addDescendant(r);
                
                activeNodes[whichNode] = q; // replace p with tip node q
                activeNodes.push_back(r); // push r as current tip
                
            } else { // extinction
                
                // replace current node p with last node in array then pop last node from array
                activeNodes[whichNode] = activeNodes[activeNodes.size()-1];
                activeNodes.pop_back();
                
                
            }
        }
    }
    
    // set index and time for tip nodes
    for (int i = 0; i < activeNodes.size(); i++) {
        activeNodes[i]->setTime(t);
        activeNodes[i]->setIndex(i);
        activeNodes[i]->setIsTip(true);
    }
    
    numTaxa = (int)activeNodes.size();
    
    initializeDownPassSequence();
    
    // initalize branch lengths
    for (Node* n : downPassSequence) {
        Node* nAncs = n->getAncestor();
        if (nAncs != NULL)
            n->setBrLen(n->getTime() - nAncs->getTime());        
    }
    
    // set index for all internal nodes
    int idx = (int)activeNodes.size();    
    for (Node* p : downPassSequence) {
        if (p->getIndex()==-1)
            p->setIndex(idx++);
    }
            
    rescale();
    
}


Tree::~Tree(void) {
    for (int i = 0; i < (int)nodes.size(); i++)
        delete nodes[i];
}

Node* Tree::addNode(void) {
    Node* newNode = new Node;
    newNode->setOffset((int)nodes.size()-1);
    nodes.push_back(newNode);
    return newNode;
}

Node* Tree::addNode(int idx) {
    Node* newNode = addNode();
    newNode->setOffset((int)nodes.size()-1);
    newNode->setIndex(idx);
    return newNode;
}

void Tree::initializeDownPassSequence(void) {
    downPassSequence.clear();
    passDown(root);
}

// tree traversal
void Tree::passDown(Node* p) {
    if (p != NULL) {
        std::set<Node*>& pDescendants = p->getDescendants();
        for (Node* n : pDescendants)
            passDown(n);
        downPassSequence.push_back(p);
    }
}


void Tree::print(void) {
    listNode(root, 3);
}

void Tree::listNode(Node* p, int indent) {
    if (p != NULL) {
        std::set<Node*>& pDescendants = p->getDescendants();
        for (int i = 0; i < indent; i++)
            std::cout << " ";
        
        if (p->getAncestor() != NULL)
            std::cout << p->getIndex() << " ( a:" << p->getAncestor()->getIndex() << " ";
        else
            std::cout << p->getIndex() << " ( " << "a:NULL" << " ";
        for (Node* n : pDescendants)
            std::cout << n->getIndex() << " ";
        std::cout << ") ";
        if (p->getDest() != NULL)
            std::cout << "dest: " << p->getDest()->getIndex() << " ";
        if (p->getSource() != NULL)
            std::cout << "source: " << p->getSource()->getIndex() << " ";
        std::cout << p->getTime() << ", " << p->getBrLen() << " ";
        std::cout << p->getNewick() << " ";
        if (p == root)
            std::cout << "Root ";
        std::cout << std::endl;
        for (Node* n : pDescendants)
            listNode(n, indent+3);
    }
}


std::vector<std::string> Tree::parseNewickString(std::string ns) {
    
    std::vector<std::string> tokens;
    std::string word = "";
    
    for (int i = 0; i < ns.length(); i++) {
        char c = ns[i];
        
        if (c == '(' || c == ')' || c == ':' || c == ',' || c == ';') {
            
            std::string tok = std::string(1,c);
            if (word != "")
                tokens.push_back(word);
            word = "";
            tokens.push_back(tok);
            
        } else {
            
            std::string tok = std::string(1,c);
            word += tok;
            
        }
        
        
    }
    
#   if 0
    
    for (int i = 0; i < tokens.size(); i++)
        std::cout << i << "--" << tokens[i] << std::endl;
    
#   endif
    
    return tokens;

}

std::vector<Node*> Tree::nodesAtTime(double time) {

    
    std::vector<Node*> ndes;
    
    for (int i = (int)downPassSequence.size()-1; i >= 0; i--) {
        
        Node* p = downPassSequence[i];
        if (p->getAncestor() == NULL)
            p->setTime(0.0);
        else {
            //p->setTime(p->getAncestor()->getTime()+p->getBrLen());
            //std::cout << "p->getTime(): " << p->getTime() << ", p->getAncestor()->getTime(): " << p->getAncestor()->getTime() << std::endl;
            if ((p->getTime() > time) && (p->getAncestor()->getTime() < time)) {
                //std::cout << "time: " << time <<  ", p->getTime(): " << p->getTime() << ", p->getAncestor()->getTime(): " << p->getAncestor()->getTime() << std::endl;
                ndes.push_back(p);
            }
        }
        
    }
    //std::cout << "nodes at time: " << time << std::endl;
    //for (Node* n : ndes)
    //    std::cout << n->getIndex() << " ";
    //std::cout << std::endl;
    //print();
    
    return ndes;
}


void Tree::addSharingEvents(RandomVariable* rng, double rate, std::vector<Node*>& sourceNodes, double delta) {
        
    // insert nodes into tree at source points
            
    for (Node* n : downPassSequence) {
                        
        Node* nAncs = n->getAncestor();
        
        if (nAncs != NULL) {
            double v = nAncs->getTime();
            double len = n->getBrLen()+v;
            
            while (v < len) {
                
                nAncs = n->getAncestor();
                
                v += -log(rng->uniformRv())/rate;
                
                if (v < len) {
                    Node* p = addNode();
                    p->setTime(v);
                    //std::cout << p->getTime() << std::endl;
                    sourceNodes.push_back(p);
                                        
                    nAncs->addDescendant(p);
                    p->addDescendant(n);
                    nAncs->removeDescendant(n);
                    n->setAncestor(p);
                    p->setAncestor(nAncs);
                                            
                    Node* d = NULL;
                    p->setDest(d);
                }
            }
        }
    }
    
    initializeDownPassSequence();
    
    //std::cout << "source node size: " << sourceNodes.size() << std::endl;
    // add destination nodes to tree
    
    for (Node* n : sourceNodes) {
        reindex();
        double nTime = n->getTime();
        std::vector<Node*> activeNodes = nodesAtTime(nTime);
        std::vector<Node*>::iterator it = find(activeNodes.begin(), activeNodes.end(), n);
        if (it != activeNodes.end())
            Msg::error("Found n " + std::to_string(n->getIndex()) +  " in list! (Not good!)");
        //std::cout << "size of activeNodes: " << activeNodes.size() << std::endl;
        //if (activeNodes.size() == 1)
          //  std::cout << "n's index is " << n->getIndex() << " and time is: " << nTime << std::endl;
                
        if (activeNodes.size() > 0) {
        
            Node* d = NULL;
                                    
            // get the distances from n to all the nodes in activeNodes
            
            std::vector<double> destProbs;
            
            for (int i = 0; i < activeNodes.size(); i++) {
                double dist = findDistance(n, activeNodes[i]); // call function
                dist -= activeNodes[i]->getTime() - nTime;
                if (activeNodes[i]->getTime() - nTime < 0.0)
                    Msg::error("Trying to subtract a negative amount!");
                //std::cout << "findDistance(" << n->getIndex() << ", " << activeNodes[i]->getIndex() << "): " << dist << std::endl;
                destProbs.push_back(dist);
            }
            
            
            
            // calculate the (unscaled) probability to each potential destination
            
            double sum = 0.0;
            
            for (int i = 0; i < destProbs.size(); i++) {
                                
                double x = 1.0 / pow(destProbs[i], delta);
                sum += x;
                destProbs[i] = x;
                
            }
            
            // scale those probabilities
            
            double scaleFactor = 1.0 / sum;
            
            for (int i = 0; i < destProbs.size(); i++)
                destProbs[i] *= scaleFactor;
                            
            // then choose one potential destination in proportion to their probabilities
            
            sum = 0.0;
            double u = rng->uniformRv();
            for (int i = 0; i < destProbs.size(); i++) {
                sum += destProbs[i];
                if (u < sum) {
                    d = activeNodes[i];
                    break;
                }
            }
            
            if (d == NULL)
                Msg::error("Why the hell is d NULL?!");
            
            //std::cout << "d's index: " << d->getIndex() << std::endl;
            
            // hold distances in a vector, make a vector of doubles
            

            //do {
              //   d = activeNodes[(int)(rng->uniformRv()*activeNodes.size())];
            //} while (d == n);
            
            Node* dAncs = d->getAncestor();
            
            if (dAncs == NULL)
                Msg::error("Why the hell is dAncs NULL?!");
            
            Node* dest = addNode();
            dest->setTime(nTime);
            
            d->setAncestor(dest);
            dest->setAncestor(dAncs);
            dAncs->removeDescendant(d);
            dAncs->addDescendant(dest);
            dest->addDescendant(d);
            
            dest->setSource(n);
            n->setDest(dest);
            //std::cout << "destination is: " << dest << std::endl;
            
            initializeDownPassSequence();
            //std::cout << "source: " << n->getIndex() << " " << n->getTime() << ", dest: " << d->getIndex() << std::endl;
            
        }
        
    }
                    

    
    // initialize branch lengths from time
    
    for (Node* n : downPassSequence) {
        
        Node* nAncs = n->getAncestor();
        if (nAncs != NULL)
            n->setBrLen(n->getTime()-nAncs->getTime());
        else
            n->setBrLen(0.0);
            
    }
    
    reindex();
    //print();
    
    
    
}

void Tree::reindex() {

    // check indices of tips
    std::set<int> tipIndices;
    int nT = 0;
        
    for (Node* n : downPassSequence) {
        if (n->getIsTip() == true) {
            tipIndices.insert(n->getIndex());
            nT++;
        }
    }

    if (nT != tipIndices.size())
        Msg::error("Tip indices should be the same size as nT!");
    
    for (int i = 0; i < nT; i++) {
        std::set<int>::iterator it = tipIndices.find(i);
        if (it == tipIndices.end())
            Msg::error("Cannot find tip index!");
        
    }
    
    int interiorIdx = nT;
    for (Node* n : downPassSequence) {
        if (n->getNumDescendents() != 0) {
            n->setIndex(interiorIdx++);
        }
    }
}

/*
std::vector<Node*> desc;
desc.push_back(root);
int i = (int)downPassSequence.size()-1;
root->setIndex(i--);


while (desc.size() > 0) {
    Node* n = desc.back();
    desc.pop_back();
    std::set<Node*> descend = n->getDescendants();
    if (descend.size() != 0) {
        for (Node* m : descend) {
            m->setIndex(i--);
            desc.push_back(m);
        }
    }
}

*/


std::string Tree::getNewickString(void) {
        
    markNodesForPrinting();
    std::stringstream ss;
    
    writeTree(root, ss);
    std::string str = ss.str();
    
    int numLeft = 0, numRight = 0;
    for (int i = 0; i < str.length(); i++) {
        if (str[i]=='(')
            numLeft++;
        else if (str[i]==')')
            numRight++;
    }
        
    if (numLeft != numRight)
        Msg::error("numLeft does not equal to numRight!");
    
    if (numLeft != (numTaxa-1))
        Msg::error("numLeft does not equal to numTaxa - 1! " + std::to_string(numLeft) + " " + std::to_string(numTaxa) + " " + str);
        
    return str;
    
}


void Tree::writeTree(Node* p, std::stringstream& ss) {
        
    bool useFullNames = false;
    
    if (p != NULL) {
        
        std::set<Node*>& descendants = p->getDescendants();
        
        int numMarkedDescendants = 0;
        
        for (Node* n : descendants) {
            if (n->getNewick() == true)
                numMarkedDescendants++;
        }
        
        //std::cout << numMarkedDescendants << std::endl;
        
        if (p->getIsTip() == true) {
            
            if (useFullNames == true)
                ss << p->getName() << ":" << getLength(p);
            else
                ss << p->getIndex()  << ":" << getLength(p);
        } else {
            if (p->getNewick() == true && numMarkedDescendants == 2)
                ss << "(";
        }

        int i = 0;
            
        for (Node* n : descendants) {
            writeTree(n, ss);
            if ((i + 1) != (int)descendants.size() && p->getNewick() == true && numMarkedDescendants == 2)
                ss << ",";
            i++;
        }
        
        if (p->getIsTip() == false && p->getNewick() == true && numMarkedDescendants == 2)
                ss << ")" << ":" << getLength(p);
    }
}
    
    
double Tree::getLength(Node* p) {
        
    Node* n = p;

    double v = 0.0;
    
    do {
        
        v += n->getBrLen();
        n = n->getAncestor();
        
    } while(n != NULL && numMarkedDescendants(n) == 1);
    
    return v;
}
    
int Tree::numMarkedDescendants(Node* p) {
            
    std::set<Node*>& descendants = p->getDescendants();
    int numMarkedDescendants = 0;
    
    for (Node* n : descendants) {
        if (n->getNewick() == true)
            numMarkedDescendants++;
    }
    
    return numMarkedDescendants;
}
    
    
void Tree::markNodesForPrinting(void) {
    
    for (Node* n : downPassSequence) {
        n->setNewick(false);
    }
    
    // mark nodes for printing newick string
    
    for (Node* p : downPassSequence) {
        if (p->getIsTip() == true) {
            Node* q = p;
            while (q != NULL) {
                q->setNewick(true);
                q = q->getAncestor();
            }
        }
    }
}

void Tree::rescale(void) {
    
    for (Node* n : downPassSequence) {
        n->setRescale(false);
    }
        
    for (Node* p : downPassSequence) {
        if (p->getIsTip() == true) {
            Node* q = p;
            
            while (q != NULL) {
                q->setRescale(true);
                q = q->getAncestor();
            }
        }
    }
    
    Node* newRoot = NULL;
    double subtreeLength = 0.0;
    
    for (int i = (int)downPassSequence.size()-1; i >= 0; i--) {
        
        Node* p = downPassSequence[i];
        std::set<Node*>& descendants = p->getDescendants();
        int numMarkedDescendants = 0;
        for (Node* n : descendants) {
            if (n->getRescale() == true)
                numMarkedDescendants++;
        }
        
        if (numMarkedDescendants == 2 && newRoot == NULL)
            newRoot = p;
                    
        if (newRoot == NULL && numMarkedDescendants == 1)
            p->setRescale(false);
        
        if (p->getRescale() == true)
            subtreeLength += p->getBrLen();
        
    }
            
    //root = newRoot;
    //root->setAncestor(NULL);
    //double min = root->getTime();
    
    //initializeDownPassSequence();
    //reindex();
    
    //for (Node* p : downPassSequence) {
      //  p->setTime( (p->getTime()-min) /(1.0-min) );
    //}
    
}

// find distance from LCA
double Tree::findDistance(Node* p, Node* q) {
    
    double distance = 0.0;
    
    for (Node* n : downPassSequence) {
        n->setIsDistance(false);
    }

    Node* r = p;
    //distance = r->getTime();
    while (r != NULL) {
        r->setIsDistance(true);
        if (r->getAncestor() != NULL) {
            distance += r->getTime() - r->getAncestor()->getTime();
            //std::cout << "P " << r->getIndex() << " : r->getTime(): " << r->getTime() << ", r->getAncestor()->getTime(): " << r->getAncestor()->getTime() << ", minus: " << r->getTime() - r->getAncestor()->getTime() <<  ", distance: " << distance << std::endl;
            // accumulate time
        }
        r = r->getAncestor();
    }
     
    r = q;
    
    while (r != NULL) {
        if (r->getAncestor() != NULL) {
            if (r->getIsDistance() == true) {
                distance -= r->getTime() - r->getAncestor()->getTime();
                //std::cout << "Q- " << r->getIndex() << " : r->getTime(): " << r->getTime() << ", r->getAncestor()->getTime(): " << r->getAncestor()->getTime() << ", minus: " << r->getTime() - r->getAncestor()->getTime() <<  ", distance: " << distance << std::endl;
            } else {
                distance += r->getTime() - r->getAncestor()->getTime();
                //std::cout << "Q+ " << r->getIndex() << " : r->getTime(): " << r->getTime() << ", r->getAncestor()->getTime(): " << r->getAncestor()->getTime() << ", minus: " << r->getTime() - r->getAncestor()->getTime() <<  ", distance: " << distance << std::endl;
            }
            
        }
        r = r->getAncestor();
    }
    //print();
    return distance;
    
    
}


Node* Tree::getNodeByIdx(int idx) {
    
    for (Node* n : downPassSequence) {
        if (n->getIndex() == idx)
            return n;
    }
    return NULL;
}
