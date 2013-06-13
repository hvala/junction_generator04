
#ifndef CHROMOSOME003_H
#define CHROMOSOME003_H

#include <vector>
#include <algorithm>

#include "junction.h"
#include "cnode.h"

using namespace std;

extern int generation;
extern ofstream errorCheck;
extern vector<Junction*> junctionPool;


// COMMENTS ON THE CHROMOSOME CLASS
//
// Background notes:
// The chromosome in these simuations provides a means and structure for storing and manipulation
// information about ancestry throughout the genome. As such, they are thought of as ordered lists
// of junctions. Using junctions, we can ignore intevening genotypic information, thus making the
// chromsomes much smaller than they might be in other types of forward simulations. They also provide
// a scaffold upon which recombination is to occur, as well as being the main structure that is
// inherited by an individual's offspring
//
// This program will also feature special types of chromosomes: the sex chromosomes as well as cytoplasmic
// inheritance mechanism. Sex chromosomes will contain loci that indicate the sex of individuals, and will
// also not recombine in heterogametic sexes. Cytoplasmic factors will be inherited in single copy, either
// maternally or paternally depending on the biology of the organism being studied
//
// Implementation notes:
// Each chromosome is implemented (roughly speaking) as a circular, doubly-linked list. The chromosome object,
// itself, is not this list but merely the starting point for operating on the list. Three parameters are required
// to instantiate an chromosome object. They are:
//      1) centromere - a pointer to a CNode which indicates the 0 M position on a chromosome(or chromosome arm--
//                note: the program was developed with mouse genomes in mind, of which all chromosomes are
//                      telocentric. This will be adapted in later versions to accomodate other chromosomal
//                      architectures -- metacentric, etc.
//      2) telomere - a pointer to CNode which indicates the distal end of the chromosome
//      3) type - a variable of type ChrType ( an enum defined in jgenum003.h) indicating the type of chromosome
//              i.e. autosomal vs. sex chromosome or cytoplasmic. The type may also imply the length.
//      4) length - the genetic length of the chromosome in Morgans, also the position of the telomere
//
// Of the above, 1, 2, and 3 are required to instantiate a chromosome. 4 is implied by three.
//
// Chromosomes also provide functions to calculate basic summary statistics, such as number of junctions, and
// the average tract length of ancestry tracts on the chromosome. They can also make use of markers used for
// interrogating ancestry at specific positions in the genome. In addition, chromosomes provide functions for
// displaying the information on them.
//
//

enum ChrType { A, W, X, Y, Z, M, C, CP };

class Chromosome
{
    public:
        // Chromosome constructor prototype
        Chromosome(CNode * cen, CNode * tel, int n, double l, ChrType t, double p);
        // Chromosome copy constructor prototype
        Chromosome(Chromosome& c);
        // Chromosome destructor prototype
        virtual ~Chromosome() {};

        // Setters
        void setCentromere(CNode* c) { centromere = c; }
        void setTelomere(CNode* t) { telomere = t; }

        // Getters
        CNode *getCentromere() { return centromere; }
        CNode *getTelomere() { return telomere; }
        virtual double getLength() {};
        ChrType getType() { return type; }
        virtual int getNumber() {};
        virtual double getParB() {};

        // Biological functions
        virtual Chromosome* duplicateChr() {};
        virtual int positionAnc(double p) {};  // return the ancestry at a position p (any position)

        // Summary Statistic Calcuators
        int calc_NumJunctions();
        double calc_ATL();

        // Data management and display
        virtual void displayChromosome() {};
        virtual void textChromosome() {};
        vector<int> junByWindow(double w);
        //void fs_tractLengths();
        void del_Nodes();
        CNode** dup_Nodes();

        // Operators
        //Chromosome& operator = (Chromosome& c);

    protected:
        CNode* centromere;      // pointer to the centromere node
        CNode* telomere;        // pointer to the telomere node
        double length;          // the length of the chromosome
        ChrType type;
        int number;
        double parB;
};

// CHROMOSOME CONSTRUCTOR
Chromosome::Chromosome(CNode * cen, CNode * tel, int n, double l, ChrType t, double p ):
    centromere(cen), telomere(tel), number(n), length(l), type(t), parB(p)
{
    centromere = cen;
    telomere = tel;
}


// CHROMOSOME COPY CONSTRUCTOR
Chromosome::Chromosome( Chromosome& c):
    length(c.length), type(c.type), number(c.number), parB(c.parB)
{
    CNode** newCT = c.dup_Nodes();
    centromere = newCT[0];
    telomere = newCT[1];
}
/*
// OVERLOAD ASSIGNMENT OPERATOR
Chromosome& Chromosome::operator = (Chromosome& c)
{
    if ( this != &c)
    {
        del_Nodes();
        centromere = c.getCentromere();
        telomere = c.getTelomere();
        length = c.getLength();
        type = c.getType();
    }

    return *this;
}*/

//REMOVE ALL THE NODES ON THE CHROMOSOME FROM MEMORY
void Chromosome::del_Nodes()
{
    CNode *chrNode = centromere;

    while ( (*chrNode).getDistP() != 0 )
    {
        CNode *delNode = chrNode;
        chrNode = (*delNode).getDistP();
        delete delNode;
    }

    delete chrNode;
}

//DUPLICATE THE NODES ON A CHROMOSOME
CNode** Chromosome::dup_Nodes()
{
    CNode** dupCT = new CNode* [2];

    dupCT[0] = new CNode( (*centromere).getJunction() );  // make a new CNode that is a copy of the centromere
    dupCT[1] = new CNode( (*telomere).getJunction() );    // copy the telomere
    //Chromosome * repChr = new Chromosome(newCen, newTel, length);       // make a new chromosome using the copied centromere and telomere
    CNode * tempNode = (*dupCT[0]).getDistP();               // make pointer that starts on the CNode after the centromere on the template chromosome
    CNode * repNode = dupCT[0] ;               // make a point that starts at the centromere on the replicate chromosome

    while ( (*tempNode).getDistP() != 0 )
    {
        CNode * nextNode = new CNode( (*tempNode).getJunction() );  // copy the new chrNode to the replicate chromosome
        (*nextNode).setProxP(repNode);                              // link the new node to the preceeding replicate node
        (*repNode).setDistP(nextNode);
        repNode = (*repNode).getDistP();                            // advance the replicate Node, and the template Node
        tempNode = (*tempNode).getDistP();
    }

    (*repNode).setDistP( dupCT[1] );             // when at the telomere, link the last node replicated to the new telomere
    (*dupCT[1]).setProxP(repNode);

    return dupCT;
}

// COUNT THE NUMBER OF JUNCTIONS ON THE CHROMOSOME
int Chromosome::calc_NumJunctions()
{
    int n = 0;
    CNode chrNode = *centromere;

    while ( chrNode.getDistP() != 0 )
    {
        if ( chrNode.getJunction() != 0 )
        {
            n++;
        }
        CNode * nextNode = chrNode.getDistP();
        chrNode = *nextNode;
    }

    return n;
}

// COUNT THE NUMBER OF JUNCTION PER WINDOW ON THE CHROMOSOME
vector<int> Chromosome::junByWindow(double w)
{
    vector<int> windows;
    CNode * node = centromere;
    int i = 1;

    for ( i = 1 ; i * w < length ; i++ )
    {
        int n = 0;
        do
        {
            node = (*node).getDistP();
            if ( (*node).getJPosition() < i * w )
            {
                n++;
            }

        } while ( (*node).getJPosition() <= i * w && (*node).getDistP() != 0 );

        windows.push_back(n);
    }

    return windows;
}

// CALCULATE THE AVERAGE ANCESTRY TRACT LENGTH ON THE CHROMOSOME
double Chromosome::calc_ATL()
{
    int n = calc_NumJunctions();
    double atl = length / n;
    return atl;
}
/*
// DISPLAY THE CHROMOSOME'S LIST OF JUNCTIONS IN THE TERMINAL WINDOW
void Chromosome::displayChromosome()
{
    CNode *chrNode = centromere;
    Junction *jun;

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
void Chromosome::textChromosome()
{
    CNode *chrNode = centromere;
    Junction *jun;

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
*/

/*
// PRINT THE TRACT LENGTHS OF A CHROMOSOME TO A .CSV FILE
double Chromosome::fs_tractLengths()
{
    int n = 0;
    double atl = 0;
    CNode chrNode = *centromere;
    Junction *prevJ = (*centromere).getJunction();
    double prevPosition = (*prevJ).getPosition();

    CNode * nextNode = chrNode.getDistP();
    chrNode = *nextNode;

    while ( chrNode.getDistP() != 0 )
    {
        Junction *jun = chrNode.getJunction();
        if ( chrNode.getJunction() != 0 )
        {
            n++
        }



        CNode * nextNode = chrNode.getDistP();
        chrNode = *nextNode;
        prevJ = (*chrNode).getJunction();
        prevPosition = (*prevJ).getPosition();
    }


    return atl;
}*/


#endif // CHROMOSOME003_H
