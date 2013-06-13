
#ifndef EPISTATICMARKER_H
#define EPISTATICMARKER_H


using namespace std;

#include <cmath>
#include <fstream>
#include <iostream>

#include "chromosome.h"



class EpiMarker
{
    public:
        EpiMarker(vector<int> c, vector<double> p, double s, int d, int t);
        ~EpiMarker();

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
        int g_Chr(Chromosome* c, double p);
        double dominantSel (vector<int> a);
        double additiveSel (vector<int> a);
        double recessiveSel (vector<int> a);
        double underDominantSel (vector<int> a);
        double overDominantSel (vector<int> a);

        // Landscape Builiding
        double* build_2DD();
        double* build_2RR();
        double* build_2DR();
        //double* build_custom(int n, vector<double> w);

    private:
        int numLoci;
        vector<int> chromosomes;
        vector<double> positions;
        double selCoeff;            // A selection coefficient used to calculate the reduction in fitness, depends on type how its used
        int dominance;              // Rec-Rec = 0, Dom-Rec = 1, Dom-Dom = 2
        int selType;                // DEVELOPMENTAL = 0,  EVIRONMENTAL = 1, REPRODUCTIVE = 2
        vector<double> envSel;
        double* landscape;
};

// EPIMARKER CONSTRUCTOR
EpiMarker::EpiMarker(vector<int> c, vector<double> p, double s, int d, int t):
    chromosomes(c), positions(p), selCoeff(s), dominance(d), selType(t)
{
    numLoci = chromosomes.size();

    switch(d)
    {
        case 2:
            landscape = build_2DD();
            break;
        case 1:
            landscape = build_2DR();
            break;
        case 0:
            landscape = build_2RR();
            break;
        default:
            landscape = build_2RR();
            break;
    }
}

// EPIMARKER DESTRUCTOR
EpiMarker::~EpiMarker()
{
    delete [] landscape;
}

// SELECTION APPLICATION FUNCTIONS

// GENOTYPE THE POSITION ON BOTH CHROMOSOMES OF AN INDIVIDUAL AND ALTER THE FITNESS OF THE INDIVIDUAL ACCORIDINGLY
void EpiMarker::sel_Ind(Individual* i)
{
    vector<int> ancs;

    for ( int a = 0; a < numLoci ; a++ )
    {
        Chromosome* chr = (*i).getChrOneP();
        int ancOne = g_Chr(chr, positions[a]);

        chr = (*i).getChrTwoP();
        int ancTwo = g_Chr(chr, positions[a]);

        int ancIndex = ancOne + ancTwo;
        ancs.push_back(ancIndex);
    }

    int landscapeIndex = 0;
    vector<int>::iterator iterA;
    for ( iterA = ancs.begin() ; iterA < ancs.end() ; iterA++ )
    {
        int b = 0;
        landscapeIndex += (*iterA) * pow(3,b);
        b++;
    }
    double fitRedux = landscape[landscapeIndex];

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
double EpiMarker::dominantSel(vector<int> ancs)
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
double EpiMarker::additiveSel(vector<int> ancs)
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
double EpiMarker::recessiveSel(vector<int> ancs)
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
double EpiMarker::underDominantSel(vector<int> ancs)
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
double EpiMarker::overDominantSel(vector<int> ancs)
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
int EpiMarker::g_Chr(Chromosome* c, double p)
{
    int anc = (*c).positionAnc(p);
    return anc;
}

// FITNESS LANDSCAPE BUILDING
// BUILD A FITNESS LANDSCAPE FOR A TWO-LOCUS DOMINANT-DOMINANT DMI
double* EpiMarker::build_2DD()
{
    int x = pow(3, 2);
    double* l = new double [x];
    for ( int a = 0 ; a < x ; a++ )
    {
        *(l+a) = 1;
    }

    l[5] -= selCoeff;  l[6] -= selCoeff;    l[8] -= selCoeff;    l[9] -= selCoeff;

}

// BUILD A FITNESS LANDSCAPE FOR A TWO-LOCUS RECESSIVE-RECESSIVE DMI
double* EpiMarker::build_2RR()
{
    int x = pow(3, 2);
    double* l = new double [x];
    for ( int a = 0 ; a < x ; a++ )
    {
        *(l+a) = 1;
    }

    l[9] -= selCoeff;
}

// BUILD A FITNESS LANDSCAPE FOR A TWO-LOCUS DOMINANT-RECESSIVE DMI
double* EpiMarker::build_2DR()
{
    int x = pow(3, 2);
    double* l = new double [x];
    for ( int a = 0 ; a < x ; a++ )
    {
        *(l+a) = 1;
    }

    l[6] -= selCoeff;  l[9] -= selCoeff;
}
/*
// BUILD A CUSTOM FITNESS LANDSCAPE FOR ARBITRARY NUMBERS OF LOCI AND GENETIC ARCHITECTURE
double* build_custom(int n, vector<double> w)
{
    int x = pow(3, n);
    double* l = new double [x];

    for ( int a = 0 ; a < x ; a++ )
    {
        *(l+a) = w[a];
    }
}*/
#endif // EPISTATICMARKER_H
