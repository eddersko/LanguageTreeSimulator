//
//  main.cpp
//  Simulator
//
//  Created by Edwin Ko on 3/3/22.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <map>
#include "Tree.hpp"
#include "CharMatrix.hpp"
#include "RandomVariable.hpp"
#include "Node.hpp"
double** makeRateMatrix(void);

int main(int argc, const char * argv[]) {
    
    double expectedNumberOfTips =  20.0;
    double turnoverRate = 0.1; // mu over lambda
    double diversificationRate = log(expectedNumberOfTips) - log(2);
    double mu = diversificationRate * turnoverRate / (1.0 - turnoverRate);
    double lambda = mu + diversificationRate;
                
    double** q = makeRateMatrix();
    double freqs[2] = {0.5, 0.5};
    
    // parameters
        
    double alphaR = 1.0;
    double alphaS = 1.0;
    double betaS = 1.0;
    //double sharingRate = 0.0;
    //double delta = 0.0;
    int numCognates = 1000;
    
    
                     
    // numcognates" + std::to_string(numCognates) + "-tr" + std::to_string(turnoverRate) + "-alphaR" + std::to_string(alphaR) + "-alphaS" + std::to_string(alphaS) + "-betaS" + std::to_string(betaS) + "-sr-delta.nex", std::ofstream::out);
    
    int numTrees = 1000;
    std::vector<Tree*> trees;
    //Tree* t = new Tree(lambda, mu, 1.0);
    //trees.push_back(t);
    
    /*double sum = 0.0;
    
    while ((int)trees.size() < numTrees) {
        Tree* t = new Tree(lambda, mu, 1.0);
        for (Node* n : t->getDownPassSequence()) {
            if (n->getIsTip() == true)
                sum++;
        }
        trees.push_back(t);
    }
            
    std::cout << "number of taxa: " << sum / 1000 << std::endl;*/
    
    
    while ((int)trees.size() < numTrees) {
        Tree* t = new Tree(lambda, mu, 1.0);
        if (t->getNumTaxa() == (int)expectedNumberOfTips)
            trees.push_back(t);
        else
            delete t;
    }
     
    
    //std::map<Tree*, CharMatrix*> simulatedMatrices;
    
    for (double sharingRate = 1.00; sharingRate <= 1.00; sharingRate+=0.05) {
        
        
        std::ofstream fw("/Users/edwinko/Downloads/treedist_numcognates" + std::to_string(numCognates) + "-tr" + std::to_string(turnoverRate) + "-alphaR" + std::to_string(alphaR) + "-alphaS" + std::to_string(alphaS) + "-betaS" + std::to_string(betaS) + "-sr" + std::to_string(sharingRate) + "-delta.nex", std::ios::out);
        std::ofstream gw("/Users/edwinko/Downloads/brlen_numcognates" + std::to_string(numCognates) + "-tr" + std::to_string(turnoverRate) + "-alphaR" + std::to_string(alphaR) + "-alphaS" + std::to_string(alphaS) + "-betaS" + std::to_string(betaS) + "-sr" + std::to_string(sharingRate) + "-delta.nex", std::ios::out);
        
        for (double delta = -10.0; delta <=  10.0; delta+=0.5) {
            
        int count = 1;
        
        for (Tree* tree : trees) {
            Tree t(*tree);
            
            fw << std::fixed << std::setprecision(2) << std::endl;
            fw << "#NEXUS\n\n";
            fw << "[!replicate " << std::to_string(count) << "]\n";
            fw << "begin data;\n";
            fw << "dimensions nchar=" << std::to_string(numCognates) << " ntax=" << std::to_string(t.getNumTaxa()) << ";\n";
            fw << "format missing=? gap=- datatype=standard symbols=\"01\";\n";
            fw << "Matrix\n\n";
            
            gw << std::fixed << std::setprecision(2) << std::endl;
            gw << "#NEXUS\n\n";
            gw << "[!replicate " << std::to_string(count) << "]\n";
            gw << "begin data;\n";
            gw << "dimensions nchar=" << std::to_string(numCognates) << " ntax=" << std::to_string(t.getNumTaxa()) << ";\n";
            gw << "format missing=? gap=- datatype=standard symbols=\"01\";\n";
            gw << "Matrix\n\n";
            
            CharMatrix* cm = new CharMatrix(&t, q, 2, freqs, numCognates, alphaR, alphaS, betaS, sharingRate, delta);
            
            //simulatedMatrices.insert(std::make_pair(tree, cm));
                            
            //std::cout << tree->getNewickString() << std::endl;
            //tree->print();
            //break;
            
            fw << cm->getString();
            fw << ";\n";
            fw << "end;\n\n";
            fw << "begin trees;\n";
            fw << "tree true_tree=" << t.getNewickString() << ";\n";
            fw << "end;\n\n";
            fw << "begin paup;\n";
            fw << "set autoclose=yes warnreset=no increase=auto;\n";
            fw << "deroot;\n";
            fw << "savetrees file=/Users/edwinko/Downloads/temp.nex replace=yes format=altnex;\n";
            fw << "set criterion=parsimony;\n";
            fw << "hsearch;\n";
            fw << "savetrees file=/Users/edwinko/Downloads/temp.nex format=altnex append=yes;\n";
            fw << "gettrees file=/Users/edwinko/Downloads/temp.nex allBlocks=yes;\n";
            if (count == 1)
                                                                
                fw << std::fixed << std::setprecision(2) << "treedist reftree=1 file=/Users/edwinko/Downloads/treedist_result_nc" << numCognates << "_tr" << turnoverRate << "_alphaR" << alphaR << "_alphaS" << alphaS << "_betaS" << betaS << "_sr" << sharingRate << "_delta" << delta << ".txt replace=yes;\n";
            else
                fw << std::fixed << std::setprecision(2) << "treedist reftree=1 file=/Users/edwinko/Downloads/treedist_result_nc" << numCognates << "_tr" << turnoverRate << "_alphaR" << alphaR << "_alphaS" << alphaS << "_betaS" << betaS << "_sr" << sharingRate << "_delta" << delta << ".txt append=yes;\n";
            fw << "end;\n";
                        
            gw << cm->getString();
            gw << ";\n";
            gw << "end;\n\n";
            gw << "begin trees;\n";
            gw << "tree true_tree=" << t.getNewickString() << ";\n";
            gw << "end;\n\n";
            gw << "begin paup;\n";
            gw << "set autoclose=yes warnreset=no increase=auto;\n";
            gw << "set criterion=likelihood;\n";
            gw << "lset clock=yes genFreq=equal;\n";
            gw << "lscores 1;\n";
                                    
            if (count == 1)
                gw << std::fixed << std::setprecision(2) << "savetrees file=/Users/edwinko/Downloads/brlen_result_nc" << numCognates << "_tr" << turnoverRate << "_alphaR" << alphaR << "_alphaS" << alphaS << "_betaS" << betaS << "_sr" << sharingRate << "_delta" << delta << ".txt replace=yes format=altnex brLens=yes ;\n";
            else
                gw << std::fixed << std::setprecision(2) << "savetrees file=/Users/edwinko/Downloads/brlen_result_nc" << numCognates << "_tr" << turnoverRate << "_alphaR" << alphaR << "_alphaS" << alphaS << "_betaS" << betaS << "_sr" << sharingRate << "_delta" << delta << ".txt append=yes format=altnex brLens=yes ;\n";
            gw << "end;\n";
            
            count++;
            delete cm;
            
        }
            std::cout << "sr: " << sharingRate << ", delta: " << delta << std::endl;
        }
        
        fw.close();
        gw.close();
    }



    
    return 0;
}


double** makeRateMatrix (void) {
    
    double** m = new double*[2];
    m[0] = new double[4];
    for (int i = 1; i < 2; i++)
        m[i] = m[i-1]+2;
    
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++) {
            if (i == j)
                m[i][j] = -1.0;
            else
                m[i][j] = 1.0;
        }
    
    return m;
}
