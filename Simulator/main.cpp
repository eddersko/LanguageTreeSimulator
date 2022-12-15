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
double** makeDolloRateMatrix(std::vector<double> bf);
//std::vector<double> makeDolloBaseFreq(double lambda, double mu);

int main(int argc, const char * argv[]) {
    
    // parameters
        
    double alphaR = 1.0;
    double alphaS = 1.0;
    double betaS = 1.0;
    int numCognates = 1000;
    int numTrees = 1;
    
    // numcognates" + std::to_string(numCognates) + "-tr" + std::to_string(turnoverRate) + "-alphaR" + std::to_string(alphaR) + "-alphaS" + std::to_string(alphaS) + "-betaS" + std::to_string(betaS) + "-sr-delta.nex", std::ofstream::out);
    
    double expectedNumberOfTips = 20.0;
    
    for (double sharingRate = 0.0; sharingRate <= 5.0; sharingRate += 1.0) {
    
    for (double turnoverRate = 0.0; turnoverRate <= 0.8; turnoverRate += 0.2) {
    
    //double turnoverRate = 0.8; // mu over lambda
    double diversificationRate = log(expectedNumberOfTips) - log(2);
    double mu = diversificationRate * turnoverRate / (1.0 - turnoverRate);
    double lambda = mu + diversificationRate;    
    double ratio = 1.0;
        
    std::vector<double> freqs = {0.9, 0.1};
                    
    //std::vector<double> freqs = makeDolloBaseFreq(lambda, mu);
    // std::vector<double>
    //double arrayOfFreqs [5][2] = {{0.5, 0.5}, {0.6, 0.4}, {0.7, 0.3}, {0.8, 0.2}, {0.9, 0.1}};
    std::string modelName = "sd";
    
    double** q = makeRateMatrix();
    
    std::vector<Tree*> trees;
    Tree* exTree = NULL;
        
    while ((int)trees.size() < numTrees) {
        Tree* t = new Tree(lambda, mu, 1.0);        
        if (t->getNumTaxa() == (int)expectedNumberOfTips)
            trees.push_back(t);
        else
            delete t;
    }
            
    //while (exTree == NULL) {
    //    Tree* t = new Tree(lambda, mu, 1.0);
        //    if (t->getNumTaxa() == (int)expectedNumberOfTips)
        //         exTree = t;
        // else
        //  delete t;
        //}
                
        //for (int i = 0; i < 5;  i++) {
    
        
            //std::cout << std::to_string(freqs[0]) << ", " << std::to_string(freqs[0]) << std::endl;
            
    //std::map<Tree*, CharMatrix*> simulatedMatrices;
                    //
            std::ofstream fw("/Users/edwinko/Dropbox/Berkeley/Computational Phylolinguistics/borrowing/data/tr/" + modelName +  "_treedist_numtaxa" + std::to_string(expectedNumberOfTips) + "-numcognates" + std::to_string(numCognates) + "-bf" + std::to_string(freqs[0]) + "_" + std::to_string(freqs[1]) + "-tr" + std::to_string(turnoverRate) + "-alphaR" + std::to_string(alphaR) + "-alphaS" + std::to_string(alphaS) + "-betaS" + std::to_string(betaS) + "-srRatio" + std::to_string(ratio)  + "-sr" + std::to_string(sharingRate) + "-delta.nex", std::ios::out);
            std::ofstream gw("/Users/edwinko/Dropbox/Berkeley/Computational Phylolinguistics/borrowing/data/tr/" + modelName +  "_brlen_numtaxa" + std::to_string(expectedNumberOfTips) + "-numcognates" + std::to_string(numCognates) + "-bf" + std::to_string(freqs[0]) + "_" + std::to_string(freqs[1]) + "-tr" + std::to_string(turnoverRate) + "-alphaR" + std::to_string(alphaR) + "-alphaS" + std::to_string(alphaS) + "-betaS" + std::to_string(betaS) + "-srRatio" + std::to_string(ratio)  + "-sr" + std::to_string(sharingRate) + "-delta.nex", std::ios::out);
                              
                for (double delta = -10.0; delta <= 10.1; delta += 5) {
                    
                int count = 1;
                
                for (Tree* tree : trees) {
                    
                    Tree t1(*tree);
                    
                    fw << std::fixed << std::setprecision(2) << std::endl;
                    fw << "#NEXUS\n\n";
                    fw << "[!replicate " << std::to_string(count) << "]\n";
                    fw << "begin data;\n";
                    fw << "dimensions nchar=" << std::to_string(numCognates) << " ntax=" << std::to_string(t1.getNumTaxa()) << ";\n";
                    fw << "format missing=? gap=- datatype=standard symbols=\"01\";\n";
                    fw << "Matrix\n\n";
                    
                    gw << std::fixed << std::setprecision(2) << std::endl;
                    gw << "#NEXUS\n\n";
                    gw << "[!replicate " << std::to_string(count) << "]\n";
                    gw << "begin data;\n";
                    gw << "dimensions nchar=" << std::to_string(numCognates) << " ntax=" << std::to_string(t1.getNumTaxa()) << ";\n";
                    gw << "format missing=? gap=- datatype=standard symbols=\"01\";\n";
                    gw << "Matrix\n\n";
                      
                    
                    //CharMatrix* cm = new CharMatrix(&t1, q, 2, freqs, numCognates, alphaR, alphaS, betaS, sharingRate, delta);
                    CharMatrix* cm = new CharMatrix(&t1, q, 2, freqs, numCognates, alphaR, alphaS, betaS, sharingRate, ratio, delta);
                    //simulatedMatrices.insert(std::make_pair(tree, cm));
                  
                    fw << cm->getString();
                    fw << ";\n";
                    fw << "end;\n\n";
                    fw << "begin trees;\n";
                    fw << "tree true_tree=" << t1.getNewickString() << ";\n";
                    fw << "end;\n\n";
                    fw << "begin paup;\n";
                    fw << "set autoclose=yes warnreset=no increase=auto;\n";
                    fw << "deroot;\n";
                    fw << "savetrees file=/Users/edwinko/Downloads/temp.nex replace=yes format=altnex;\n";
                    fw << "set criterion=parsimony;\n";
                    if (freqs[0] == 0.9)
                        fw << "hsearch swap=none;\n";
                    else
                        fw << "hsearch swap;\n";
                    //fw << "pset collapse=no;\n";
                    fw << "savetrees file=/Users/edwinko/Downloads/temp.nex format=altnex append=yes;\n";
                    fw << "gettrees file=/Users/edwinko/Downloads/temp.nex allBlocks=yes;\n";

                    if (count == 1)
                        fw << std::fixed << std::setprecision(2) << "treedist reftree=1 file=/Users/edwinko/Downloads/" + modelName +  "_treedist_result_numTaxa" << expectedNumberOfTips << "_nc" << numCognates << "_bf" << freqs[0] << "-" << freqs[1] << "_tr" << turnoverRate << "_alphaR" << alphaR << "_alphaS" << alphaS << "_betaS" << betaS << "_ratio" << ratio << "_sr" << sharingRate << "_delta" << delta << ".txt replace=yes;\n";
                        
                    else
                        fw << std::fixed << std::setprecision(2) << "treedist reftree=1 file=/Users/edwinko/Downloads/" + modelName +  "_treedist_result_numTaxa" << expectedNumberOfTips << "_nc" << numCognates << "_bf" << freqs[0] << "-" << freqs[1] << "_tr" << turnoverRate << "_alphaR" << alphaR << "_alphaS" << alphaS << "_betaS" << betaS << "_ratio" << ratio << "_sr" << sharingRate  << "_delta" << delta << ".txt append=yes;\n";
                    
                    fw << "end;\n";
                                
                    gw << cm->getString();
                    gw << ";\n";
                    gw << "end;\n\n";
                    gw << "begin trees;\n";
                    gw << "tree true_tree=" << t1.getNewickString() << ";\n";
                    gw << "end;\n\n";
                    gw << "begin paup;\n";
                    gw << "set autoclose=yes warnreset=no increase=auto;\n";
                    gw << "set criterion=likelihood;\n";
                    gw << "lset clock=yes genFreq=equal;\n";
                    gw << "lscores 1;\n";
                                            
                    if (count == 1)
                        
                        gw << std::fixed << std::setprecision(2) << "savetrees file=/Users/edwinko/Downloads/" + modelName + "_brlen_result_numTaxa" << expectedNumberOfTips << "_nc" << numCognates << "_bf" << freqs[0] << "-" << freqs[1] << "_tr" << turnoverRate << "_alphaR" << alphaR << "_alphaS" << alphaS << "_betaS" << betaS << "_ratio" << ratio << "_sr" << sharingRate <<  "_delta" << delta << ".txt replace=yes format=altnex brLens=yes ;\n";
                    else
                        gw << std::fixed << std::setprecision(2) << "savetrees file=/Users/edwinko/Downloads/" + modelName + "_brlen_result_numTaxa" << expectedNumberOfTips << "_nc" << numCognates << "_bf" << freqs[0] << "-" << freqs[1] << "_tr" << turnoverRate << "_alphaR" << alphaR << "_alphaS" << alphaS << "_betaS" << betaS << "_ratio" << ratio << "_sr" << sharingRate << "_delta" << delta << ".txt append=yes format=altnex brLens=yes ;\n";
                    gw << "end;\n";
                    
                    count++;
                    delete cm;
                                
                }

                    
                }
                fw.close();
                gw.close();
             }
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



/*

 Corresponding MATLAB code:
 
 vocab=lambda/mu; % ratio of birth / death rate
 L=10*round(vocab);

 x=vocab/L;
 BF=[1 x]; % base frequencies
 BF=BF./sum(BF)

 Dollo with {stationary distribution}
 Dollo with {1/2, 1/2}
 Dollo with root with 0, 0 -> 1 good, but 1-> 0 is really small
 
 PAUP: swap = none (ends search on the starting tree)
 start with default being stepwise, other options are NJ, current, or treenum
 addseq = simple, closest, (play around with them)
 
 TBR (default) more trees with neighors than NII (every tree has few neighbors)
 
 */

/*
std::vector<double> makeDolloBaseFreq(double r10, double r01) {
    //double vocab = lambda / mu;
    //double L = 10 * round(vocab);
    //double x = vocab / L;
    //double sum = 1 + x;
    std::vector<double> freqs(2);
    freqs[0] = r10/(r10 + r01);
    freqs[1] = r01/(r10 + r01);
    return freqs;
}*/

/*
    Corresponding MATLAB code:
 
    Q=[-1 1;1 -1]*diag(BF);
    Q=Q-diag(diag(Q));
    Q=Q-diag(sum(Q'));
    Q=Q./abs(BF*diag(Q));
*/


double** makeDolloRateMatrix(std::vector<double> bf) {
    
    double** m = new double*[2];
    m[0] = new double[4];
    for (int i = 1; i < 2; i++)
        m[i] = m[i-1]+2;
    
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++) {
            if (i == j && i == 0)
                m[i][j] = -1.0 * bf[i+1]; // 0 0
            else if (i == j && i == 1)
                m[i][j] = -1.0 * bf[i-1]; // 1 1
            else if (i != j && i == 0)
                m[i][j] = 1.0 * bf[j]; // 0 1
            else
                m[i][j] = 1.0 * bf[j]; // 1 0
        }
    
    double sum = 0.0;
    
    for (int i = 0; i < 2; i++)
        sum += bf[i] * m[i][i];
    
    sum = abs(sum);
    
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            m[i][j] /= sum;
    
    return m;
    
}
    


