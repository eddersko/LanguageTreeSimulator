//
//  CognateSet.hpp
//  BirthDeathSimulation
//
//  Created by Edwin Ko on 2/4/22.
//

#ifndef CognateSet_hpp
#define CognateSet_hpp
#include <vector>
#include <set>

class RandomVariable;
class CognateSet {
    
public:
    CognateSet(void) = delete;
    CognateSet(int nc, RandomVariable* rng, std::vector<double> eq);
    CognateSet(const CognateSet& c);
    void incrementNumSites(int val);
    int& operator[](size_t i) { return cognates[i] ; }
    const int& operator[](size_t i) const { return cognates[i] ; }
    double calculateRate(double** q, int numStates);
    CognateSet* getAncestor(void) {return ancestor;}
    void setAncestor(CognateSet* anc) { ancestor = anc; }
    void addDescendant(CognateSet* p) { descendants.insert(p); }
    void removeDescendant(CognateSet* p) {descendants.erase(p); }
    void removeAllDescendants(void) { descendants.clear(); }
    std::set<CognateSet*>& getDescendants(void) { return descendants; }
    void changeCognate(double** q, int numStates, RandomVariable* rng);
    std::vector<int>& getCognates(void) {return cognates;}
    void print(void);
    int getNumCognates(void) { return numCognates; }
private:
    std::vector<int> cognates;
    int numCognates;
    CognateSet* ancestor;
    std::set<CognateSet*> descendants;
    void listCognates(void);
};

#endif /* CognateSet_hpp */
