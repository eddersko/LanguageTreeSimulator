//
//  Tree.hpp
//  BirthDeathSimulation
//
//  Created by Edwin Ko on 1/14/22.
//

#ifndef Tree_hpp
#define Tree_hpp
#include <vector>
#include <string>
#include <sstream>

class RandomVariable;
class Node;

class Tree {
    
public:
    Tree(void) = delete;
    Tree(const Tree& t);
    Tree(std::string newickString);
    Tree(double lambda, double mu, double t);
    Tree(double lambda, double mu, double t, double sharingRate, double alpha, double beta, double** q, int numStates);
    ~Tree(void);
    void initializeDownPassSequence(void);
    void print(void);
    std::vector<Node*>& getDownPassSequence(void) {return downPassSequence;}
    std::vector<Node*> nodesAtTime(double time);
    void addSharingEvents(RandomVariable* rng, double rate, std::vector<Node*>& sourceNodes, double delta);
    void addSharingEvents(RandomVariable* rng, Tree* t, double rate, std::vector<Node*>& sourceNodes);
    Node* getRoot(void) { return root; }
    int getNumTaxa(void) {return numTaxa;}
    void setNumTaxa(int n) { numTaxa = n;}
    std::string getNewickString(void);
    double findDistance(Node* p, Node* q);
    Node* getNodeByIdx(int idx);
    double getMaxTime(void) {return time;}
    int getNodeLength(void) {return (int)nodes.size();}
    
private:
    std::vector<Node*> nodes;
    Node* root;
    double time;
    Node* addNode(void);
    Node* addNode(int idx);
    void passDown(Node* p);
    std::vector<Node*> downPassSequence;    
    void listNode(Node* p, int indent);
    std::vector<std::string> parseNewickString(std::string ns);
    void reindex(void);
    int numTaxa;
    void writeTree(Node* n, std::stringstream&);
    void markNodesForPrinting(void);
    void rescale(void);
    double getLength(Node* p);
    int numMarkedDescendants(Node* p);
};

#endif /* Tree_hpp */
