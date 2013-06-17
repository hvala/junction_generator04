#ifndef GENE_H
#define GENE_H

#include <vector>
#include <algorithm>

#include "chromosome.h"


class Gene
{
    public:
        Gene(int chr, double pos, vector<double> ae);


        // Getters
        int getChr() { return chromosome; }
        double getPos() { return position; }
        double getAddEffect(int i) { return add_effects[i]; }
        int getDominance(int i) { return dominance[i]; }

        // Setters


        // Biological Functions
        vector<int> genotype( vector<Chromosome*> );

    protected:
        int chromosome;
        double position;
        vector<double> add_effects;
        vector<int> dominance;

};

// GENE CONSTRUCTOR
Gene::Gene(int chr, double pos, vector<double> ae):
    chromosome(chr), position(pos), add_effects(ae)
{

}


// GET THE GENOTYPE FOR A GENE IN AN INDIVIDUAL
vector<int> Gene::genotype(vector<Chromosome*> g)
{
    vector<int> genotype;
    int x = chromosome * 2;

    genotype.push_back( (*g[x]).positionAnc(position) );
    genotype.push_back( (*g[x+1]).positionAnc(position) );

    sort( genotype.begin(), genotype.end() );

    return genotype;
}


#endif // GENE_H
