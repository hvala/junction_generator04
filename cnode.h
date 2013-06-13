
#ifndef CNODE003_H
#define CNODE003_H

#include <new>

#include "junction.h"


using namespace std;


// COMMENTS ON THE CNODE CLASS
//
// Background/Implementation notes:
// CNodes, chromosome nodes, are used to implement a chromosome as a doubly-linked list.
// They are modelled after the DLLs described in Drozdek 3rd. ed. with a couple of slight
// additions. Each CNode contains a pointer to the node in the proximal orientation on the chromosome
// and another in the distal orientation. The nodes that corresponde to centromeres and telomeres have null
// pointers in the proximal and distal orientation respectively. In addition, rather than containing their own
// data, a third pointer points to a junction in the common population "junction pool". By doing this, junctions
// can be used on multiple chromosomes, individuals, and sub-populations (demes) simultaneously. In addition,
// whenever a CNode is instantiated, it increments the number of occurences of the junction it points to, and
// decrements this variable when its destructor is called. Finally, each CNode contains a pointer back to the
// chromosome on which it resides. This allows one to access chromosomes and individuals from the junction itself.

class CNode
{
  public:
    // constructor prototype
    CNode(Junction *j, int i);
    // destructor prototype
    ~CNode();

    // setters
    void setDistP ( CNode * d ) { distP = d; }
    void setProxP ( CNode * p ) { proxP = p; }
    void setJunction (Junction * j);

    // getters
    Junction* getJunction() { return junctionP; }
    double getJPosition() { return (*junctionP).getPosition(); }
    double getJAncestry() {  return (*junctionP).getAncestry(); }
    CNode * getDistP() { return distP; }
    CNode * getProxP() { return proxP; }
    CNode* newAddyAssign();


    // new node pool
    static CNode* s_newAddy;
    void deallocate();

    // operators
    //CNode operator++();     // move to the next node in the distal direction on the chromosome
    //CNode operator--();     // move to the next node in the proximal direction on the chromosome
    //void* operator new(size_t size);
    //void operator delete(void*);
    // methods


  private:
    // variables
    Junction *junctionP;
    CNode *distP;
    CNode *proxP;
    int location;
};

// INITIALIZE STATIC VARIABLE S_NEWADDY
CNode* CNode::s_newAddy = 0;

// CNODE CONSTRUCTOR
CNode::CNode(Junction *j, int i = 0):
    junctionP(j), location(i)
{
    if (junctionP != 0 )
    {
        (*junctionP).incNumOccur();
        (*junctionP).incLocCount(location);
    }
    proxP = 0;
    distP = 0;
}

// CNODE DESTRUCTOR -- makes sure to decrement the number of occurences of the junction
CNode::~CNode()
{
    deallocate();
}

// MAKE THE CNODE NULL -- POINT TO NOTHING -- BEFORE SENDING IT TO THE GARBAGE COLLECTOR
void CNode::deallocate()
{
    if ( junctionP != 0 )
    {
        (*junctionP).deallocCN(location);
    }
    junctionP = 0;
}

// GET NEXT PRE-ALLOCATED CNODE ADDRESS -- OR ALLOCATE NEW ONES AND RETURN THE LAST
CNode* CNode::newAddyAssign()
{
    CNode* oldAddy;

    if ( s_newAddy != 0)
    {
        oldAddy = s_newAddy;
        s_newAddy = (*s_newAddy).getDistP();
    }
    else
    {
        for ( int i = 0 ; i < 1 ; i++ )
        {
            CNode* newAddy = new CNode(0,0);
            (*newAddy).setDistP(s_newAddy);
            (*newAddy).setProxP(0);
            s_newAddy = newAddy;
        }
        oldAddy = s_newAddy;
        s_newAddy = (*s_newAddy).getDistP();
    }

    return oldAddy;
}


// JUNCTION SETTER FUNCTION -- if the junction is reset, we need to make sure the change in occurence is accounted for
void CNode::setJunction (Junction * j)
{
    (*junctionP).decNumOccur();
    (*junctionP).decLocCount(location);
    junctionP = j;
    if ( j != 0 )
    {
        (*junctionP).incNumOccur();
        (*junctionP).incLocCount(location);
    }

}

#endif // CNODE003_H
