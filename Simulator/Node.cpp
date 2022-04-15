//
//  Node.cpp
//  BirthDeathSimulation
//
//  Created by Edwin Ko on 1/14/22.
//

#include "Node.hpp"
#include "CognateSet.hpp"

Node::Node(void) {
    index = 0;
    ancestor = NULL;
    time = 0.0;
    name = "";
    v = 0.0;
    destination = NULL;
    source = NULL;
    isTip = false;
    isNewick = false;
    isRescale = false;
    isDistance = false;
    myCognate = NULL;
    offset = 0;
}

Node::~Node(void) {
    if (myCognate != NULL)
        delete myCognate;
}
 
/*double Node::getTotalBrLen(void) {
    
    Node* pAncs = ancestor;
    double totalBrLen = v;
    
    while (pAncs != NULL) {
        totalBrLen += pAncs->getBrLen();
        pAncs = pAncs->getAncestor();
    }
    
    return totalBrLen;
    
}
*/
