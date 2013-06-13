#ifndef GENE_H
#define GENE_H

#include <vector>
#include <algorithm>

#include "chromosome.h"


class Gene
{
    public:
        Gene(int chr, double pos, int phen, vector<double> ae, vector<int> d);
        virtual ~Gene();

        // Getters
        int getChr() { return chromosome; }
        double getPos() { return position; }
        //int getPhen() { return phenotype; }
        double getAddEffect(int i) { return add_effect[i]; }
        int getDominance(int i) { return dominance[i]; }

        // Setters


        // Biological Functions
        vector<int> genotype( vector<Chromosome*> );
        double getPhenotype( vector<Chromosome*> );

    protected:
        int chromosome;
        double position;
        //int phenotype;
        vector<double> add_effects;
        vector<int> dominance;
        //Landscape landScape;

};

// GENE CONSTRUCTOR
Gene::Gene(int chr, double pos, int phen, vector<double> ae, vector<int> d):
    chromosome(chr), position(pos), phenotype(phen), add_effects(ae), dominance(d)
{

}


// GENE DESTRUCTOR
Gene::~Gene()
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

// GET THE PHENOTYPE FOR THE INDIVIDUAL BASED ON THE GENOTYPE OF THE GENE
double Gene::getPhenotype(vector<Chromosome*> g)
{
    vector<int> g_type = genotype(g);

}

#endif // GENE_H
