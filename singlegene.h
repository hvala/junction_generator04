#ifndef SINGLEGENE_H
#define SINGLEGENE_H

#include "chromosome.h"


class SingleGene : public Gene
{
    public:
        SingleGene(int chr, double pos, int phen, double pv, int s, int n, Gene* i);
        virtual ~SingleGene();

        // Getters
        int getChr() { return chromosome; }
        double getPos() { return position; }
        int getPhen() { return phenotype; }
        double getPhenVal() { return phen_value; }
        virtual Gene* getInteractor() { return 0; }

        // Setters
        virtual void setInteractor() {}

        // Biological Functions
        virtual int genotype() {}
        virtual double phenotype() {}

    private:
        int* numAlleles;
}

// SINGLE LOCUS GENE CONSTRUCTOR
SingleGene::SingleGene(int chr, double pos, int phen, double pv, int s, int n = 1, Gene* i = 0):
    Gene(chr, pos, phen, pv, s, n, i)
{
    int x = 0;
    for ( int j = 0 ; j < s ; j++ )
    {
        x += j;
    }

    f_landscape = new double [x];


}

// SINGLE LOCUS GENE DESTRUCTOR
SingleGene::~SingleGene
{

}

// GET THE GENOTYPE FOR THE LOCUS
int SingleGene::genotype(Chromosome* c, Chromosome* d)
{
    int gOne = (*c).positionAnc(position);
    int gTwo = (*d).positionAnc(position);

}

// DETERMINE THE PHENOTYPIC VALUE OF THE TRAIT DETERMINED BY THIS GENE
double SingleGene::phenotype()
{

}

#endif // SINGLEGENE_H
