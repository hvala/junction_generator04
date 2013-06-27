
#ifndef JUNGEN_ACC.CPP
#define JUNGEN_ACC.CPP

#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>

#include "gamete.h"
#include "landscape.h"
#include "chromosome.h"
#include "autosome.h"
#include "sex_chromosome.h"
#include "cytoplasm.h"

using namespace std;

extern vector<Junction*> junctionPool;

// Take the sum over integers from 1 to n
int sigmaNum(int n)
{
    int sum = 0;
    for ( int i = 1 ; i <= n ; i++ )
    {
        sum += i;
    }
    return sum;
}

// calculate a factorial
int factorial (int n)
{
    int fact = n;
    for( int m = n - 1 ; m > 1 ; m-- )
    {
        fact = fact * m;
    }
    return fact;
}

// partial factorial -- calculate a product over n to k
int part_factorial (int n, int k)
{
    int p_fact = n;
    for ( int m = n - 1 ; m > k ; m-- )
    {
        p_fact = p_fact * m;
    }
    return p_fact;
}

// calculate n choose k
int n_choose_k (int n, int k)
{
    int nck = 0;
    nck = part_factorial(n, n-k) / factorial(k);
    return nck;
}

// express an integer as a string from a different base system -- base < 10
vector<int> change_base(int x, int base, int numDigits)
{
    vector<int> newNumber;

    // find the largest base power to go into the x
    int i = 0;

    do
    {
        ++i;
    } while ( pow(base, i) <= x );

    for( numDigits ; numDigits > i ; numDigits-- )
    {
        newNumber.push_back(0);
    }

    --i;

    // count the number of times the next greatest base power fits into the number
    while ( i >= 0 )
    {
        int digit = 0;
        do
        {
            digit++;
        } while ( pow(base, i) * digit <= x );
        newNumber.push_back(--digit);

        x = x - ( pow(base, i) * digit );

        --i;
    }

    return newNumber;
}

// take a vector of digits representing a number in a base system 'b', and convert it to a decimal number
int convert_to_decimal(vector<int> number, int base)
{
    vector<int>::iterator iterD;
    int decimal = 0;
    int i;
    for(iterD = number.begin(), i = number.size() - 1 ; iterD < number.end() ; ++iterD, --i)
    {
        decimal = decimal + (*iterD) * pow(base, i);
    }
    return decimal;
}

// convert a vector of ints to a string
string vec_ints_to_string(vector<int> i)
{
    stringstream str;
    vector<int>::iterator iterI;
    for( iterI = i.begin() ; iterI < i.end() ; iterI++ )
    {
        str << (*iterI);
    }

    string key = str.str();
    return key;
}

// Treat the effect of epistatic interactions as recessive
double recEdges(int e, double s, int f = 4)
{
    double v;
    if ( e == f ) { v = s;}
    else { v = 0; }
    return v;
}

// Treat the effect of epistatic interactons as dominant
double domEdges(int e, double s, int f = 4)
{
    double v;
    if ( e > 0 ) { v = s; }
    else { v = 0; }
    return v;
}

// Treat the effect of epistatic interactions as additive
double addEdges(int e, double s, int f = 4)
{
    double edges = e;
    double total = f;
    double v = s * ( edges / total );
    return v;
}

//MAKE A GAMETE POOL FOR A POPULATION WITH XY CHROMOSOMES
vector<Gamete*> XYgametePool(int anc, vector<double> lengths, vector<ChrType> types, double p = 0.001)
{
    vector<Gamete*> gametePool;

    for ( int i = 0 ; i < anc ; i++ )
    {
        int chrN = 0;
        vector<Chromosome*> chrPool;
        Chromosome* hetGamChr;

        Junction* cen = new Junction(chrN, 0, i);
        junctionPool.push_back(cen);
        Junction* tel = new Junction(chrN, lengths[0], i);
        junctionPool.push_back(tel);


        CNode *cent1 = new CNode(cen);
        (*cent1).setProxP(0);
        CNode *telo1 = new CNode(tel);
        (*telo1).setDistP(0);

        (*cent1).setDistP(telo1);
        (*telo1).setProxP(cent1);

        hetGamChr = new SexChromosome(cent1, telo1, chrN, lengths[0], Y);
        chrPool.push_back(hetGamChr);


        vector<ChrType>::iterator iterT;
        vector<double>::iterator iterL;
        chrN++;

        for(iterT = types.begin(), iterL = lengths.begin() ; iterT < types.end() , iterL < lengths.end() ; iterT++, iterL++ )
        {
            Junction* cen = new Junction(chrN, 0, i);
            junctionPool.push_back(cen);
            Junction* tel = new Junction(chrN, *iterL, i);
            junctionPool.push_back(tel);

            CNode *cent1 = new CNode(cen);
            (*cent1).setProxP(0);
            CNode *telo1 = new CNode(tel);
            (*telo1).setDistP(0);

            (*cent1).setDistP(telo1);
            (*telo1).setProxP(cent1);

            Chromosome* newChr;
            if ( *iterT == X )
            {
                newChr = new SexChromosome(cent1, telo1, chrN, *iterL, *iterT, p);

            }
            else if ( (*iterT) == M || (*iterT) == C || (*iterT) == CP )
            {
                newChr = new Cytoplasm(cent1, telo1, chrN, 0, *iterT, 1);

            }
            else
            {
                newChr = new Autosome(cent1, telo1, chrN, *iterL, *iterT);

            }



            chrPool.push_back(newChr);

            chrN++;
        }



        vector<Chromosome*> maleGenome;
        maleGenome.push_back(chrPool[0]);



        vector<Chromosome*>::iterator iterC;
        for( iterC = chrPool.begin() + 2 ; iterC < chrPool.end() ; iterC++ )
        {
            maleGenome.push_back(*iterC);
        }

        Gamete* maleGamete = new Gamete(maleGenome);
        gametePool.push_back(maleGamete);

        vector<Chromosome*> femaleGenome;
        for( iterC = chrPool.begin() + 1 ; iterC < chrPool.end() ; iterC++ )
        {
            femaleGenome.push_back(*iterC);
        }

        Gamete* femaleGamete = new Gamete(femaleGenome);
        gametePool.push_back(femaleGamete);

    }

    return gametePool;
}

//MAKE A GAMETE POOL FOR A POPULATION WITH ZW CHROMOSOMES

//MAKE A GAMETE POOL FOR A POPULATION OF DIPLOIDS WITHOUT SEX CHROMOSOMES

#endif //JUNGEN_ACC.CPP
