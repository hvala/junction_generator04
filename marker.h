#ifndef MARKER_H
#define MARKER_H

#include <cmath>
#include <fstream>
#include <iostream>

#include "chromosome.h"
#include "individual.h"
#include "deme.h"
#include "metaPop.h"

using namespace std;

extern ofstream errorCheck;
extern ofstream os_clines;
extern int sim;

class Marker
{
    public:
        // CONSTRUCTOR PROTOTYPE
        Marker(int c, double p);

        // Getters
        double getPosition() { return mark_position; }
        int getChromosome() { return mark_chromosome; }
        vector<double> getDemeFreq() { return demeFreqs; }
        double getMPFreq() {return mpFreq; }

        // Setters
        void setPosition ( double p ) { mark_position = p; }
        void setChromosome ( int c ) { mark_chromosome = c; }
        void setDemeFreqs ( vector<double> f ) { demeFreqs = f; }

        // Genotyping functions
        int g_Chr(Chromosome* c);                    // genotype the position on a single chromosome
        vector<int> g_Ind(Individual* i);           // genotype the position on both chromosomes of an individual
        double g_Deme(Deme* d);                     // determine the ancestry frequency in given deme
        void g_Demes(MetaPop* mp);                   // determine the ancestry frequencies in all demes of a metapopulation -- returns it to demeFreq
        void g_Meta(MetaPop* mp);                  // determine the ancestry frequency in a metapopulation
        void gVar_Meta();               // determine the variance in ancestry freq among demes in a metapopulation

        // Data Management functions
        void textMarker();

    private:
        int mark_chromosome;            // the chromosome the marker is on
        double mark_position;           // the physical/genomic position of the marker
        vector<double> demeFreqs;       // a vector of marker ancestry frequencies, one for each deme
        double mpFreq;                  // the marker ancestry frequency in the meta-population
        double mpVarFreq;               // the variance in ancestry frequency among demes

};

// MARKER CONSTRUCTOR
Marker::Marker(int c = -1, double p = -1):
    mark_chromosome(c), mark_position(p)
{

}

// GENOTYPING FUNCTIONS

// GENOTYPE THE MARKER ON A SINGLE CHROMOSOME
int Marker::g_Chr(Chromosome* c)
{
    int anc = (*c).positionAnc(mark_position);
    return anc;
}

// GENOTYPE THE POSITION ON BOTH CHROMOSOMES OF AN INDIVIDUAL
vector<int> Marker::g_Ind(Individual* i)
{
    vector<int> ancs;
    Chromosome* chr = (*i).getChrOneP();
    int ancOne = g_Chr(chr);
    ancs.push_back(ancOne);
    chr = (*i).getChrTwoP();
    ancOne = g_Chr(chr);
    ancs.push_back(ancOne);
    return ancs;
}


// DETERMINE THE ANCESTRY FREQUENCY AT THE MARKER IN A GIVEN DEME
double Marker::g_Deme(Deme* d)
{
    vector<Individual*> dMembers = (*d).getMembers();
    vector<Individual*>::iterator iterI;
    double ancFreq = 0;
    int i = 0;

    for( iterI = dMembers.begin() ; iterI < dMembers.end() ; iterI++ )
    {
        vector<int> ancs = g_Ind( *iterI );
        ancFreq += ( ancs[0] + ancs[1] );
        i++;
    }

    ancFreq = ancFreq / ( 2 * dMembers.size() );

    return ancFreq;
}

// DETERMINE THE ANCESTRY FREQUENCIES OF ALL DEMES OF A METAPOPULATION
void Marker::g_Demes(MetaPop* mp)
{
    vector<Deme*>::iterator iterD;

    for ( iterD = (*mp).getDemes().begin() ; iterD < (*mp).getDemes().end() ; iterD++ )
    {
        double ancFreq = g_Deme( (*iterD) );
        demeFreqs.push_back(ancFreq);
    }
}

// DETERMINE THE ANCESTRY FREQUENCY IN A METAPOPULATION
void Marker::g_Meta(MetaPop* mp)
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

// PRINT THE MARKER, AND ITS INFO TO A TEXT FILE
void Marker::textMarker()
{
    //errorCheck << "Marker: " << mark_chromosome << "_" << mark_position << "\tgeneration" << generation << endl;
    //errorCheck << mpFreq << "\t" << mpVarFreq << endl;

    vector<double>::iterator iter;
    os_clines << sim << ",";
    for( iter = demeFreqs.begin() ; iter < demeFreqs.end() ; iter++ )
    {
        if ( iter == demeFreqs.begin() )
        {
            os_clines << (*iter);
        }
        else
        {
            os_clines << "," << (*iter) ;
        }
    }
    os_clines << endl; //<< "***************" << endl;
}

#endif // MARKER_H
