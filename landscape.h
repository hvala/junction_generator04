#ifndef LANDSCAPE_H
#define LANDSCAPE_H

#include <string>
#include <map>
#include <cmath>

#include "gene.h"
#include "interaction_graphs.h"
#include "jungen_acc.cpp"


class Landscape
{
    public:
        Landscape( double a, int na, vector<Gene*> l, vector<IntGraph*> iG );

    //getters


    //setters


    //landscape functions
        Phenotype findPhenotype(int i) { return l_Map[i]; }

    private:
        double exp_value;
        int numAncestries;
        vector<Gene*> loci;
        vector<IntGraph*> epiInt;
        map<string,double> l_Map;

};

// LANDSCAPE CONSTRUCTOR
Landscape::Landscape( double a, int na, vector<Gene*> l, vector<IntGraph*> iG ):
    exp_value(a), numAncestries(na), loci(l), epiInt(iG)
{
    int landSize = pow(sigmaNum(numAncestries), loci.size() );

    for (int i = 0 ; i < landSize ; i++ )
    {
        vector<int> key = change_base(i, numAncestries);
        string newKey = vec_ints_to_string(key);

        // iterate over the key and calculate the phenotype
        vector<int>::iterator iterK, iterL;
        vector<Gene*>::iterator iterG1, iterG2 = loci.begin();
        vector<IntGraph*>::iterator iterIG;

        double g_addEff = 0;
        double g_epiEff = 0;

        // get the additive effects
        for ( iterK = key.begin(), iterL = iterK + 1 ; iterL < key.end() ; iterK+=2 )
        {
            g_addEff += (**iterG1).getAddEffect(*iterK) + (**iterG1).getAddEffect(*iterL);
            ++iterG1;
        }



        for( iterIG = epiInt.begin() ; iterIG < epiInt.end() ; iterIG++ )
        {
            int x, y = 0;
            iterG1 = loci.begin();
            iterG2 = loci.begin();

            while( iterG1 < loci.end() && (*iterG1) != (**iterIG).getLocusA() && (*iterG1) != (**iterIG).getLocusB() )
            {

            }




        }
    }

}


#endif // LANDSCAPE_H
