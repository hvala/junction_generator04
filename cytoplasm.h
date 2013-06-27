
#ifndef CYTOPLASM_H
#define CYTOPLASM_H

#include "chromosome.h"

using namespace std;

class Cytoplasm : public Chromosome
{
    public:
        Cytoplasm(CNode* cen, CNode* tel, int n, double l, ChrType t, double p );
        Cytoplasm(Cytoplasm& a);
        virtual ~Cytoplasm();

        // getters
        double virtual getParB() { return 0.0; }
        double virtual getLength() { return 0.0; }
        ChrType virtual getType() {return type; }


        // display functions
        void virtual displayChromosome();
        void virtual textChromosome();

        // Operators
        //Cytoplasm& operator = ( Cytoplasm& c );
        virtual Chromosome* duplicateChr();


    private:

};

// AUTOSOME CONSTRUCTOR
Cytoplasm::Cytoplasm(CNode * cen, CNode * tel, int n, double l = 0, ChrType t = M, double p = 1):
    Chromosome(cen, tel, n, 0, t, 1)
{

}

// AUTOSOME DESTRUCTOR
Cytoplasm::~Cytoplasm()
{
    (*telomere).setDistP( CNode::s_newAddy );     // Add the old CNodes to the list of allocated memory spaces for CNode

    CNode::s_newAddy = centromere;                  // Set the first node pointer to the centromere's location

    CNode *chrNode = centromere;

    while ( (*chrNode).getDistP() != 0 )
    {
        CNode *dumpNode = chrNode;
        chrNode = (*dumpNode).getDistP();
        (*dumpNode).deallocate();
    }

    (*chrNode).deallocate();
}

// CHROMOSOME COPY CONSTRUCTOR
Cytoplasm::Cytoplasm(Cytoplasm& a):
    Chromosome(a)
{

}

// DUPLICATE A CYTOPLASM TO BE INHERITED
Chromosome* Cytoplasm::duplicateChr()
{
    CNode * newCen = new CNode( (*centromere).getJunction() );          // make a new CNode that is a copy of the centromere
    CNode * newTel = new CNode( (*telomere).getJunction() );            // copy the telomere
    Chromosome * repChr = new Cytoplasm(newCen, newTel, number, length, type);       // make a new chromosome using the copied centromere and telomere
    CNode * tempNode = (*centromere).getDistP();                        // make pointer that starts on the CNode after the centromere on the template chromosome
    CNode * repNode = (*repChr).getCentromere();                        // make a point that starts at the centromere on the replicate chromosome

    while ( (*tempNode).getDistP() != 0 )
    {
        CNode * nextNode = new CNode( (*tempNode).getJunction() );  // copy the new chrNode to the replicate chromosome
        (*nextNode).setProxP(repNode);                              // link the new node to the preceeding replicate node
        (*repNode).setDistP(nextNode);
        repNode = (*repNode).getDistP();                            // advance the replicate Node, and the template Node
        tempNode = (*tempNode).getDistP();
    }

    (*repNode).setDistP( (*repChr).getTelomere() );                 // when at the telomere, link the last node replicated to the new telomere
    (*newTel).setProxP(repNode);

    return repChr;
}

/*
// OVERLOAD ASSIGNMENT OPERATOR TO CALL THE COPY CONSTRUCTOR
Cytoplasm& Cytoplasm::operator=(Cytoplasm& c)
{
    if ( this != &c)
    {
        del_Nodes();
        centromere = c.getCentromere();
        telomere = c.getTelomere();
        length = c.getLength();
        type = c.getType();
        number = c.getNumber();
    }

    return *this;
}*/

// DISPLAY THE CHROMOSOME'S LIST OF JUNCTIONS IN THE TERMINAL WINDOW
void Cytoplasm::displayChromosome()
{
    CNode *chrNode = centromere;
    Junction *jun;

    cout << "[" << type << "]_";

    while ( (*chrNode).getDistP() != 0 )
    {
        jun = (*chrNode).getJunction();
        cout << "(" << (*jun).getChromosome() << "," << (*jun).getPosition() << "," << (*jun).getAncestry() << "," << (*jun).getNumOccur() <<  ")_" ;
        CNode *nextNode = (*chrNode).getDistP();
        chrNode = nextNode;
    }

    jun = (*chrNode).getJunction();
    cout << "(" << (*jun).getChromosome() << "," << (*jun).getPosition() << "," << (*jun).getAncestry() << "," << (*jun).getNumOccur() << ")" << endl;
}

// PRINT THE CHROMOSOME'S LIST OF JUNCTIONS TO A .TXT FILE
void Cytoplasm::textChromosome()
{
    CNode *chrNode = centromere;
    Junction *jun;

    errorCheck << "[" << type << "]_";

    while ( (*chrNode).getDistP() != 0 )
    {
        jun = (*chrNode).getJunction();
        errorCheck << "(" << jun << "," ; //(*jun).getChromosome() << "," << (*jun).getPosition() << "," << (*jun).getAncestry() << "," << (*jun).getNumOccur() << ")" << endl;
        (*jun).textJunction(); errorCheck << ")" << endl;
        CNode *nextNode = (*chrNode).getDistP();
        chrNode = nextNode;
    }

    jun = (*chrNode).getJunction();
    errorCheck << "(" << jun << "," ; //(*jun).getChromosome() << "," << (*jun).getPosition() << "," << (*jun).getAncestry() << "," << (*jun).getNumOccur() << ")" << endl;
    (*jun).textJunction(); errorCheck << ")" << endl;
    errorCheck << "$$$$$$" << endl;

}

#endif //CYTOPLASM_H

