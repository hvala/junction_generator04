#ifndef METAPOP_H
#define METAPOP_H

#include <vector>
#include <algorithm>

#include "deme.h"

using namespace std;

extern gsl_rng *r;
extern vector<double> migrationRates;

class MetaPop
{
    public:
        MetaPop(int d, vector<int> s, vector<double> a, vector<int> l );
        ~MetaPop();

        // setters
        void setNumDemes(int d) { numDemes = d; }

        // getters
        int getNumDemes() { return numDemes; }
        vector<Deme*> getDemes() { return demes; }
        Deme* getDeme(int i) { return demes[i]; }

        // Hybrid zone functions
        void stepStone();
        void migrate( Deme* d, Deme* e, double m);
        void migrateSource( Deme* d, double m , int s);
        void newGeneration();
        void longRangeDisp(Deme* d);


        // Data and summaries
        void junByWindow(double w, const char* f);


    private:
        int numDemes;
        vector<Deme*> demes;
        vector<int> popSizes;
        vector<double> ancProps;
        vector<int> locations;

};

// META POPULATION CONSTRUCTOR
MetaPop::MetaPop(int d, vector<int> s, vector<double> a, vector<int> l):
    numDemes(d), popSizes(s), ancProps(a), locations(l)
{
    for( int i = 0 ; i < numDemes; i++ )
    {
        Deme* d = new Deme(popSizes[i], ancProps[i], locations[i]);
        demes.push_back(d);
    }
}

// META POPULATION DESTRUCTOR
MetaPop::~MetaPop()
{
    vector<Deme*>::iterator iterD;
    for ( iterD = demes.begin() ; iterD < demes.end() ; iterD++ )
    {
        delete *iterD;
    }
}

// MIGRATE INDIVIDUALS AMONG DEMES IN THE WHOLE HYBRID ZONE AND SOURCE POPULATIONS -- STEPPING STONE MODEL
void MetaPop::stepStone()
{
    vector<Deme*>::iterator iterD;
    vector<Deme*>::iterator iterE;
    vector<double>::iterator iterMig = migrationRates.begin();

    for ( iterD = demes.begin() ; iterD <= demes.end() ; iterD++)
    {
        if ( iterD == demes.begin() )
        {
            migrateSource( (*iterD), (*iterMig), 0);
            iterMig++;
        }
        else if ( iterD == demes.end() )
        {
            iterE = iterD - 1;
            migrateSource( (*iterE), (*iterMig), 1);
            //longRangeDisp(*iterE);
        }
        else
        {
            iterE = iterD-1;
            migrate( *iterE, *iterD, (*iterMig));
            //longRangeDisp(*iterE);
            iterMig++;
        }
    }
}

// EXCHANGE RANDOM INDIVIUDALS BETWEEN TWO DEMES -- MIGRATION
void MetaPop::migrate( Deme* d, Deme* e, double m)
{
    double mu = (*d).getSize() * m / 2;
    int migrants;

    if ( mu <= 1 )
    {
        mu = gsl_ran_poisson(r, mu);
    }

    migrants = mu;

    for ( int i  = 0 ; i < migrants ; i++ )
    {
        int x;
        int y;

        do
        {
            x = gsl_rng_uniform(r) * (*d).getSize();

        } while ( (*((*d).getMembers()[x])).migrated );

        do
        {
            y = gsl_rng_uniform(r) * (*e).getSize();

        } while ( (*((*d).getMembers()[y])).migrated );

        Individual* tempInd = (*d).getMembers()[x];
        (*d).immigrate(x, (*e).getMembers()[y]);
        (*e).immigrate(y, tempInd);
    }
}

// IMMIGRATE AN INDIVIDUAL FROM A SOURCE POPULATION
void MetaPop::migrateSource( Deme* d, double m, int s)
{
    int migrants = (*d).getSize() * m / 2;
    for ( int i  = 0 ; i < migrants ; i++ )
    {
        //Individual*& indOne = (*d).getMembers()[0];
        int x;

        do
        {
            x = gsl_rng_uniform(r) * (*d).getSize();

        } while ( (*((*d).getMembers()[x])).migrated );

        Chromosome* migChrOne = (*chromosomePool[s]).duplicateChr();
        Chromosome* migChrTwo = (*chromosomePool[s]).duplicateChr();

        Individual* migrant = new Individual(migChrOne, migChrTwo);
        delete (*d).getMembers()[x];
        (*d).immigrate(x, migrant);
    }
}

// SAMPLE MIGRANTS FROM RANDOM DEMES THROUGHOUT THE META-POPULATION -- AS PER KIMURA AND WEISS (1964), M-INFINITY TERM
void MetaPop::longRangeDisp(Deme* d)    // d receives the random deme's emigrant
{
    double lrdRate = 0.01;
    int numMigrants = lrdRate * (*d).getSize(); //gsl_ran_poisson(r, lrdRate);

    for ( int i = 0 ; i < numMigrants ; i++ )
    {
        Deme* e;                    // The randomly chosen deme that supplies the immigrant to d
        Individual* randIndImm;     // the emigrant to d
        Individual* randIndEmm;     // the emigrant to e
        int imm;
        int emm;

        do
        {
            e = demes[gsl_rng_uniform(r) * numDemes];
            imm = gsl_rng_uniform(r) * (*e).getSize();
            randIndImm = (*e).getMembers()[imm];
        }
        while ( ((*randIndImm).migrated) );

        do
        {
            emm = gsl_rng_uniform(r) * (*d).getSize();
            randIndEmm = (*d).getMembers()[emm];
        }
        while ( ((*randIndEmm).migrated) );

        Individual* tempInd = (*d).getInd(emm);
        (*d).immigrate(emm, (*e).getInd(imm) );
        (*e).immigrate(imm, tempInd);
    }
}

// MAKE A NEW GENERATION OF THE HYBRID ZONE
void MetaPop::newGeneration()
{
    (*this).stepStone();

    vector<Deme*>::iterator iterD;

    for ( iterD = demes.begin() ; iterD < demes.end() ; iterD++ )
    {
        (*(*iterD)).newGeneration();
    }
}

// PRINT THE NUMBER OF JUNCTIONS PER WINDOW OF A CHROMOSOME TO A FILE FOR EACH INDIVIDUAL IN THE POPULATION
void MetaPop::junByWindow(double winSize, const char* fileName)
{
    Individual* indName = 0;
    int demeNumber = -1;

    ofstream junctionFile;
    junctionFile.open(fileName);

    vector<Deme*>::iterator iterD;
    for( iterD = demes.begin() ; iterD < demes.end() ; iterD++ )
    {
        demeNumber = (**iterD).getLocation();

        vector<Individual*>::iterator iterI;
        vector<Individual*> demeMembers = (**iterD).getMembers();

        for ( iterI = demeMembers.begin() ; iterI < demeMembers.end() ; iterI++ )
        {
            Chromosome* chrOne = (**iterI).getChrOneP();
            Chromosome* chrTwo = (**iterI).getChrTwoP();

            vector<int> junOne = (*chrOne).junByWindow(winSize);
            vector<int> junTwo = (*chrTwo).junByWindow(winSize);

            double oneATL = (*chrOne).calc_ATL();
            int oneNJ = (*chrOne).calc_NumJunctions();
            double twoATL = (*chrTwo).calc_ATL();
            int twoNJ = (*chrTwo).calc_NumJunctions();

            junctionFile << (*iterI) << "-1," << demeNumber << "," << oneNJ << "," << oneATL ;

            vector<int>::iterator iterW;
            for( iterW = junOne.begin() ; iterW < junOne.end() ; iterW++ )
            {
                junctionFile << "," << (*iterW) ;
            }

            junctionFile << endl << (*iterI) << "-2," << demeNumber << "," << twoNJ << "," << twoATL ;

            for( iterW = junTwo.begin() ; iterW < junTwo.end() ; iterW++ )
            {
                junctionFile << "," << (*iterW) ;
            }

            junctionFile << endl;
        }
    }

    junctionFile.close();
}

#endif // METAPOP_H
