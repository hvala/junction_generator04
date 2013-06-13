#ifndef GENOTYPE_H
#define GENOTYPE_H

#include "chromosome.h"

class Genotype
{
    public:
        Genotype(int l, int a);
        ~Genotype();

        // Getters
        int getNumLoci() { return numLoci; }
        int getNumAnc() { return numAncestries; }

        //Setters


        // Genotype functions
        void genotype(Individual* i, Gene g);
        int calc_Index();


    private:
        int numLoci;
        int* numAncestries;
        int* genotype;

}

// GENOTYPE CONSTRUCTOR
Genotype::Genotype(int l, int a):
    numLoci(l), numAncestries(a)
{
    int x = numLoci * 2 + 1;
    genotype = new int[x];
    genotype[x-1] = 0;
}

// GENOTYPE DESTRUCTOR
{
    delete [] genotype;
}

// CALCULATE THE INDEX OF THE GENOTYPE
double Genotype::calc_Index()
{
    int g = 0;
    vector<int> geno;
    while( genotype[g] != 0 )
    {
        int h = genotype[g];
        geno.push_back()
    }
}

#endif // GENOTYPE_H
