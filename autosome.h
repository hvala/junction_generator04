
#ifndef AUTOSOME_H
#define AUTOSOME_H

#include "chromosome.h"

using namespace std;

class Autosome : public Chromosome
{
    public:
        Autosome(CNode* cen, CNode* tel, int n, double l, ChrType t, double p );
        Autosome(Autosome& a);
        virtual ~Autosome();

        // getters
        int virtual getNumber() { return number; }
        double virtual getLength() { return length; }
        double virtual getParB() { return 1.0; }

        // chromosome functions
        virtual int positionAnc(double p);


        // display functions
        void virtual displayChromosome();
        void virtual textChromosome();

        // Operators
        //Autosome& operator = ( Autosome& c );
        virtual Chromosome* duplicateChr();


    private:

};

// AUTOSOME CONSTRUCTOR
Autosome::Autosome(CNode * cen, CNode * tel, int n, double l, ChrType t = A, double p = 1):
    Chromosome(cen, tel, n, l, t, p)
{

}

// AUTOSOME DESTRUCTOR
Autosome::~Autosome()
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
Autosome::Autosome( Autosome& a):
    Chromosome(a)
{

}

// DUPLICATE AN AUTOSOME TO BE INHERITED
Chromosome* Autosome::duplicateChr()
{
    CNode * newCen = new CNode( (*centromere).getJunction() );          // make a new CNode that is a copy of the centromere
    CNode * newTel = new CNode( (*telomere).getJunction() );            // copy the telomere
    Chromosome * repChr = new Autosome(newCen, newTel, number, length);       // make a new chromosome using the copied centromere and telomere
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

//  GENOTYPE A POSITION ON THE CHROMOSOME
int Autosome::positionAnc(double p)
{
    int anc = -1;

    CNode* distNode = (*centromere).getDistP();

    while ( (*distNode).getJPosition() <= p && (*distNode).getDistP() != 0 )
    {
        distNode = (*distNode).getDistP();
    }

    CNode* proxNode = (*distNode).getProxP();
    anc = (*proxNode).getJAncestry();

    return anc;
}

/*
// OVERLOAD ASSIGNMENT OPERATOR TO CALL THE COPY CONSTRUCTOR
Autosome& Autosome::operator=(Autosome& c)
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
}
*/

// DISPLAY THE CHROMOSOME'S LIST OF JUNCTIONS IN THE TERMINAL WINDOW
void Autosome::displayChromosome()
{
    CNode *chrNode = centromere;
    Junction *jun;

    cout << "[" << type << ":" << number << "(" << length << "M)]_";

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
void Autosome::textChromosome()
{
    CNode *chrNode = centromere;
    Junction *jun;

    errorCheck << "[" << type << ":" << number << "(" << length << "M)]_";

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

#endif //AUTOSOME_H
