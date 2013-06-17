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
        double findPhenotype(string i) { return l_Map[i]; }
    // data management and display
        void printLandscape();

    private:
        double exp_value;
        double errorVar;
        double phenMin;
        double phenMax;
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
    cout << "b1 " << landSize << endl;
    for (int i = 0 ; i < landSize ; i++ )
    {
        vector<int> key = change_base(i, numAncestries);
        string newKey = vec_ints_to_string(key);
        double value = exp_value;
        cout << "b2 " << endl;
        // iterate over the key and calculate the phenotype
        vector<int>::iterator iterK, iterL;
        vector<Gene*>::iterator iterG1 = loci.begin();
        vector<Gene*>::iterator iterG2;
        vector<IntGraph*>::iterator iterIG;

        double g_addEff = 0;
        double g_epiEff = 0;

        // get the additive effects
        for ( iterK = key.begin(), iterL = iterK + 1 ; iterL < key.end() ; iterK+=2, iterL += 2 )
        {
            g_addEff += (**iterG1).getAddEffect(*iterK) + (**iterG1).getAddEffect(*iterL);
            cout << "b2.1 " << g_addEff << endl;
            ++iterG1;
        }
        cout << "b3 " << endl;

        for( iterIG = epiInt.begin() ; iterIG < epiInt.end() ; iterIG++ )
        {
            int w = 0;
            int x = 0;
            int y = 0;
            int z = 0;
            iterG1 = loci.begin();
            iterG2 = loci.begin();

            for (iterG1 = loci.begin() ; iterG1 < loci.end() ; iterG1++, w++ )
            {
                if ( (*iterG1) == (**iterIG).getLocusA() || (*iterG1) == (**iterIG).getLocusB() )
                {
                    y = w;  // y is now the index of the locus at one of the interacting loci in the gene list
                    for (iterG2 = loci.begin() ; iterG2 < loci.end() ; iterG2++, x++ )
                    {
                        if ( (*iterG2) == (**iterIG).getLocusA() || (*iterG2) == (**iterIG).getLocusB() )
                        {
                            z = x; // z is now the index of locus at the other interacting locus in the gene list
                        }
                    }
                }
            }
            cout << "b3 " << endl;
            // now that we have the indices we can check to see if the key indicates epistatic selection
            int a, b, c, d;

            a = key[y];
            b = key[y+1];
            c = key[z];
            d = key[z+1];

            g_epiEff += (**iterIG).calc_Selection(a,b,c,d);
            cout << "b4" << endl;
        }

        value = value + g_addEff + g_epiEff;
        l_Map[newKey] = value;
        cout << "b10" << endl;
    }

}

// print out the landscape keys and values
void Landscape::printLandscape()
{
    map<string, double>::iterator iterL;
    for(iterL = l_Map.begin() ; iterL != l_Map.end() ; iterL++ )
    {
        cout << (*iterL).first << "\t" << (*iterL).second << endl;
    }
}

#endif // LANDSCAPE_H
