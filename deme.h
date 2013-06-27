
#ifndef DEME_H
#define DEME_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <ctime>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "individual.h"
#include "chromosome.h"

using namespace std;

extern gsl_rng *r;
extern int generation;
extern vector<Chromosome*> chromosomePool;
extern bool incompatibilityOn;
extern ofstream errorCheck;
extern vector<Junction*> junctionPool;
extern int nMarkers;


class Deme
{
    public:
        // constructor prototype
        Deme(int s, vector<double> a, int l, vector<Gamete*> gP);
        ~Deme();

        // getters
        int getSize() { return size; }
        vector<double> getAncProp() { return ancProp; }
        vector<Individual*> getMembers() { return members; }
        int getLocation() { return demeLocation;}
        Individual* getInd(int i) { return members[i]; }

        // setters
        void immigrate(int i, Individual* j);

        // population functions
        void newGeneration();
        void generationSummary();
        Individual* randInd();
        void purgeMembers();

        // summary statistic calculators

        // data management
        void textMembers();

    private:
        int size;
        vector<double> ancProp;
        int demeLocation;
        vector<Individual*> members;
        int gen;
};

// DEME CONSTRUCTOR
Deme::Deme(int s, vector<double> a, int l, vector<Gamete*> gP):
    size(s), ancProp(a), demeLocation(l)
{

    unsigned int n[ancProp.size()];

    //double* a = &v[0];
    //gsl_ran_multinomial (const gsl_rng * r, size_t K, unsigned int N, const double p[], unsigned int n[])
    gsl_ran_multinomial(r, ancProp.size(), size, &ancProp[0], n);
    gen = generation;

    for( int i = 0 ; i < ancProp.size() ; i++ )
    {
        for( int j = 0 ; j < n[i] ; j++ )
        {

            //Individual(Gamete* matGamete, Gamete* patGamete, int l = 0, bool m = false ):
            int sexChr = gsl_ran_bernoulli(r, 0.5);
            int G1 = 2 * i + 1;
            int G2 = 2 * i + 1 - sexChr;
            Individual* newInd = new Individual( gP[G1], gP[G2], demeLocation );

            members.push_back(newInd);
        }
    }


/*    ancFreq = new double [nMarkers];
    for( int i = 0 ; i < nMarkers ; i++ )
    {
        ancFreq[i] = 0;
    }
*/
}

//DEME DESTRUCTOR
Deme::~Deme()
{
    vector<Individual*>::iterator iterI;
    for ( iterI = members.begin() ; iterI < members.end() ; iterI++ )
    {
        delete *iterI;
    }
}

// MAKE A NEW GENERATION IN A DEME FROM THE PREVIOUS GENERATION
void Deme::newGeneration()
{
    vector<Individual*> newMembers;
    for (int i = 0 ; i < size ; i++ )
    {
        Individual* parentOne = 0;
        Individual* parentTwo = 0;
        bool parOnePass = false;
        bool parTwoPass = false;

        // select parents at random, then accept or reject them in accord with their reproductive fitness
        do
        {
            parentOne = randInd();
            double rejectOne = gsl_rng_uniform(r);

            if ( rejectOne < (*parentOne).getRFitness() )
            {
                parOnePass = true;
            }
        }   while ( parOnePass == false );

        do
        {
            parentTwo = randInd();
            double rejectTwo = gsl_rng_uniform(r);

            if ( rejectTwo < (*parentTwo).getRFitness() )
            {
                parTwoPass = true;
            }

        }   while ( parTwoPass == false );

        // once parents have been selected, make a gamete from each, and use them to make a new individual
        Gamete* gameteOne = (*parentOne).makeGamete();
        Gamete* gameteTwo = (*parentTwo).makeGamete();

        Individual* offspring = new Individual(gameteOne, gameteTwo, demeLocation);
        //(*gameteOne).~Chromosome();
        //(*gameteTwo).~Chromosome();

        newMembers.push_back(offspring);

    }
    vector<Individual*>::iterator iterM;

    // once the new generation is complete, delete the old
    for ( iterM = members.begin(); iterM < members.end() ; iterM++ )
    {
        delete *iterM;
    }

    members = newMembers;
}

// IMMIGRATE A NEW INDIVIDUAL -- add a new indivual from elsewhere
void Deme::immigrate(int i, Individual* j)
{
    //cout << "a " << endl;
    //Individual* oldMember = members[i];
    members[i] = j;
    (*members[i]).setLocation(demeLocation);
    (*members[i]).migrated = true;              // set migrated to true so it won't migrate again this generation
    //delete oldMember;
}



// CHOOSE AN INDIVIDUAL AT RANDOM FROM THE DEME
Individual* Deme::randInd()
{
    int x = gsl_rng_uniform(r) * size;
    return members[x];
}

// PURGE THE OLD DEME FROM MEMORY
void Deme::purgeMembers()
{
    vector<Individual*>::iterator iter;
    for( iter = members.begin() ; iter < members.end() ; iter++ )
    {
        delete *iter;
    }
}

// PRINT THE MEMORY ADDRESSES OF DEME MEMBERS TO A TEXT FILE
void Deme::textMembers()
{
    vector<Individual*>::iterator iter;
    for( iter = members.begin() ; iter < members.end() ; iter++ )
    {
        errorCheck << (*iter) << endl;
        //(*(*iter)).textChromosomes(3);
    }

    errorCheck << "*****************" << endl;
}



#endif // DEME_H
