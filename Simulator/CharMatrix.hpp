//
//  CharMatrix.hpp
//  BirthDeathSimulation
//
//  Created by Edwin Ko on 1/26/22.
//

#ifndef CharMatrix_hpp
#define CharMatrix_hpp
class Tree;
class RandomVariable;
class Node;

class CharMatrix {
    
public:
    CharMatrix(void) = delete;
    ~CharMatrix(void);
    CharMatrix(Tree* t, double** q, int ns, double* freqs, int nc, RandomVariable* rng);
    CharMatrix(double** q, int ns, double* freqs, int nc, RandomVariable* rng, double alpha, double beta, double sharingTimes);
    CharMatrix(Tree* t, double** q, int ns, std::vector<double> freqs, int nc, double alphaRat, double alphaRes, double betaRes, double sharingRate, double delta);
    CharMatrix(Tree* t1, Tree* t2, double** q, int ns, std::vector<double> freqs, int nc, double alphaRat, double alphaRes, double betaRes, double sharingRate, double ratio, double delta);
    void print(void);
    std::string getString(void);
private:
    int numTaxa;
    int numChar;
    int numStates;
    int** matrix;
    double* resilience;
    std::vector<double> rateVar;
    double* siteRate;
    void simulateSubTree(Node* n, Node* r, double** q, RandomVariable* rng, int site);
};

#endif /* CharMatrix_hpp */
