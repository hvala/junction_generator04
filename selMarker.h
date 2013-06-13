#ifndef SELECTEDMARKER_H
#define SELECTEDMARKER_H

using namespace std;

#include <cmath>
#include <fstream>
#include <iostream>

#include "chromosome.h"


class SelMarker
{
    public:
        SelMarker( int c, double p, double s, int d, int t );

        // setters
        void setSelCoeff(double s) { selCoeff = s; }
        void setDominance(int d) { dominance = d; }
        void setEnvSel(vector<double> e) { envSel = e; }

        // getters
        double getSelCoeff() { return selCoeff; }
        int getDominance() { return dominance; }
        vector<double> getEnvSel() { return envSel; }

        // Selection Application
        void sel_Ind(Individual* i);
        int g_Chr(Chromosome* c);
        double dominantSel (vector<int> a);
        double additiveSel (vector<int> a);
        double recessiveSel (vector<int> a);
        double underDominantSel (vector<int> a);
        double overDominantSel (vector<int> a);

    private:
        int chromosome;
        double position;
        double selCoeff;            // A selection coefficient used to calculate the reduction in fitness, depends on type how its used
        int dominance;              // Dominant = 2, Additive = 1, Recessive = 0, Underdominant = 3, Overdominant = 4
        int selType;                // DEVELOPMENTAL = 0,  EVIRONMENTAL = 1, REPRODUCTIVE = 2
        vector<double> envSel;
};

// SELMARKER CONSTRUCTOR
SelMarker::SelMarker(int c, double p, double s, int d, int t):
    chromosome(c), position(p), selCoeff(s), dominance(d), selType(t)
{

}
// SELECTION APPLICATION FUNCTIONS

// GENOTYPE THE POSITION ON BOTH CHROMOSOMES OF AN INDIVIDUAL AND ALTER THE FITNESS OF THE INDIVIDUAL ACCORIDINGLY
void SelMarker::sel_Ind(Individual* i)
{
    vector<int> ancs;
    Chromosome* chr = (*i).getChrOneP();
    int ancOne = g_Chr(chr);
    ancs.push_back(ancOne);
    chr = (*i).getChrTwoP();
    ancOne = g_Chr(chr);
    ancs.push_back(ancOne);
    double fitRedux = 0;

    switch(dominance)
    {
        case 4:
            fitRedux = overDominantSel(ancs);
            break;
        case 3:
            fitRedux = underDominantSel(ancs);
            break;
        case 2:
            fitRedux = dominantSel(ancs);
            break;
        case 1:
            fitRedux = additiveSel(ancs);
            break;
        case 0:
            fitRedux = recessiveSel(ancs);
            break;
        default:
            fitRedux = recessiveSel(ancs);
            break;
    }

    switch(selType)
    {
        case 0:
            (*i).setDFitness(fitRedux);
            break;
        case 1:
            fitRedux = fitRedux * envSel[(*i).getLocation()];
            (*i).setEFitness(fitRedux);
            break;
        case 2:
            (*i).setRFitness(fitRedux);
            break;
        default:
            (*i).setRFitness(fitRedux);
            break;
    }
}

// SELECTION FROM DOMINANT LOCI TO AN INDIVIDUAL
double SelMarker::dominantSel(vector<int> ancs)
{
    double s = 0;
    int ancSum = 0;
    ancSum += ( ancs[0] + ancs[1] );

    if( ancSum )
    {
        s = selCoeff;
    }

    return s;
}

// SELECTION FROM ADDITIVE LOCI TO AN INDIVIDUAL
double SelMarker::additiveSel(vector<int> ancs)
{
    double s = 0;
    int ancSum = 0;
    ancSum += ( ancs[0] + ancs[1] );

    if( ancSum == 1 )
    {
        s = selCoeff / 2;
    }
    else if ( ancSum == 2 )
    {
        s = selCoeff;
    }

    return s;
}

// SELECTION FROM ADDITIVE LOCI TO AN INDIVIDUAL
double SelMarker::recessiveSel(vector<int> ancs)
{
    double s = 0;
    int ancSum = 0;
    ancSum += ( ancs[0] + ancs[1] );

    if( ancSum == 2 )
    {
        s = selCoeff / 2;
    }

    return s;
}

// SELECTION FROM UNDERDOMINANT LOCI TO AN INDIVIDUAL
double SelMarker::underDominantSel(vector<int> ancs)
{
    double s = 0;
    int ancSum = 0;
    ancSum += ( ancs[0] + ancs[1] );

    if( ancSum == 1 )
    {
        s = selCoeff ;
    }

    return s;
}

// SELECTION FROM OVERDOMINANT LOCI TO AN INDIVIDUAL
double SelMarker::overDominantSel(vector<int> ancs)
{
    double s = 0;
    int ancSum = 0;
    ancSum += ( ancs[0] + ancs[1] );

    if( ancSum == 2 || ancSum == 0 )
    {
        s = selCoeff;
    }

    return s;
}

// GENOTYPE THE MARKER ON A SINGLE CHROMOSOME
int SelMarker::g_Chr(Chromosome* c)
{
    int anc = (*c).positionAnc(position);
    return anc;
}

/*
// APPLY SELECTION TO ALL INDIVIDUALS IN A GIVEN DEME
void SelMarker::sel_Deme(Deme* d)
{
    vector<Individual*> dMembers = (*d).getMembers();
    vector<Individual*>::iterator iterI;
    double ancFreq = 0;
    int i = 0;

    for( iterI = dMembers.begin() ; iterI < dMembers.end() ; iterI++ )
    {
        sel_Ind( *iterI );
    }
}

// DETERMINE THE ANCESTRY FREQUENCIES OF ALL DEMES OF A METAPOPULATION
void SelMarker::sel_Demes(MetaPop* mp)
{
    vector<Deme*>::iterator iterD;

    for ( iterD = (*mp).getDemes().begin() ; iterD < (*mp).getDemes().end() ; iterD++ )
    {
        sel_Deme( (*iterD) );
    }
}
*/



/*
// DETERMINE THE ANCESTRY FREQUENCY IN A METAPOPULATION
void Marker::sel_Meta(MetaPop* mp)
{
    vector<Deme*> mpDemes = (*mp).getDemes();
    vector<Deme*>::iterator iterD;

    int i = 0;
    int totalSize = 0;
    double w_demeFreq = 0;

    for ( iterD = mpDemes.begin() ; iterD < mpDemes.end() ; iterD++ )
    {
        int demeSize =  (*(*iterD)).getSize();

        w_demeFreq += demeFreqs[i] * demeSize;
        totalSize += demeSize;
        i++;
    }

    mpFreq = w_demeFreq / totalSize;
}

// DETERMINE THE VARIANCE IN ANCESTRY FREQUENCY AMONG DEMES IN A METAPOPULATION
void Marker::gVar_Meta()
{
    double avg_freq;
    double sum_squares;
    vector<double>::iterator iter;

    int i = 0;
    for( iter = demeFreqs.begin() ; iter < demeFreqs.end() ; iter++ )
    {
        avg_freq = avg_freq + (*iter);
        i++;
    }

    avg_freq = avg_freq / i;

    i = 0;
    for( iter = demeFreqs.begin() ; iter < demeFreqs.end() ; iter++ )
    {
        sum_squares = sum_squares + pow( ( *iter - avg_freq ), 2 );
        i++;
    }

    if ( i == 1)
    {
        i++;
    }

    mpVarFreq = sum_squares / ( i - 1 );
}
*/
#endif // SELECTEDMARKER_H
