#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include <sstream>
#include <string>
#include <ctime>
#include <sys/stat.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "junction.h"
#include "cnode.h"
#include "chromosome.h"
#include "individual.h"
#include "deme.h"
#include "metaPop.h"
#include "marker.h"

using namespace std;



class Experiment
{
    public:
        // Experiment constructor prototype
        Experiment();

        // Setters
        void setExpName(string n) { experimentName = n; }
        void setNumSim(int s) { numSim = s; }
        void setNumGens(int g) { numGens = g; }
        void setNumDemes( int d) { numDemes = d; }
        void setNumAnc( int a) { numAncestries = a; }
        void setDemeSizes( vector<int> s) { demeSizes = s; }
        void setAncProps( vector<double> p) { ancProportions = p; }
        void setDemeLocs( vector<int> l) { demeLocations = l; }
        void setMigRates( vector<double> m) { migRates = m; }
        void setIntOpt( int i) { intOpt = i; }
        void setNumChr( int c) { numChromosomes = c; }
        void setChrLens( vector<double> l) { chrLengths = l; }
        void setSummaries( vector<int> s) { summaries = s; }
        void setNumMarkers( int m ) { numMarkers = m; }


        // Getters
        string getExpName() { return experimentName; }
        int getNumSim() { return numSim; }
        int getNumGens() { return numGens; }
        int getNumDemes() { return numDemes; }
        int getNumAnc() { return numAncestries; }
        vector<int> getDemeSizes() { return demeSizes; }
        vector<double> getAncProps() { return ancProportions; }
        vector<int> getDemeLocs() { return demeLocations; }
        vector<double> getMigRates() { return migRates; }
        int getIntOpt() { return intOpt; }
        int setNumChr() { return numChromosomes; }
        vector<double> getChrLens() { return chrLengths; }
        vector<int> getSummaries() { return summaries; }
        int getNumMarkers() { return numMarkers; }

        // Experiment functions
        void experimentInput();         // Run an interface to collect experiment parameters from the user
        void runExperiment();           // Run the experiment
        void cleanJPool();              // Remove all junctions that no longer exist
        void drainJPool();              // Remove all junctions from the junction pool
        void textJPool();               // Print all the junctions to a text file
        vector<Gamete*> makeGametePool( int a, vector<double> l, vector<ChrType> t);         // Makes the chromosomes needed to start the simulation
        void drainChromosomePool();     // Clears the current experiment's chromosomes from the chromosome pool
        void runGenerations(MetaPop* mP);         // Iterate through the generations of the simulation

        // Summary/Statistic Functions
        void makeMarkers(int c, double l);
        void genotypeMarkers(MetaPop* mp);

        // Data Management Functions
        void textMarkers(string dir, string fname);
        const char* makeFileName( string file, string dir );

    private:
        // Experiment parameters
        string experimentName;             // the name of the scenario -- used to make files and directories
        int numSim;                     // the total number of simulations to run

        // Demographic parameters
        int numGens;                    // the number of generation per simulation
        int numDemes;                   // the number of sub-populations
        int numAncestries;              // the number of ancestries contributing to admixture
        vector<int> demeSizes;          // a list of population sizes for each deme in the hybrid zone
        vector<double> ancProportions;  // a list of ancestry proportions for each deme
        vector<int> demeLocations;
        vector<double> migRates;        // a list of migration rates between demes

        // Genomic parameters
        int intOpt;                     // indicates the type of recombination to employ
                                        // 0 for one crossover per chromosome
                                        // 1 for poisson distributed crossovers
                                        // 2 crossover interference
        int numChromosomes;             // the haploid number of chromosomes
        //vector<ChrType> chrTypes;       // the types of the chromosomes - see jgenum.h
        vector<double> chrLengths;      // the lengths of each chromosome - in Morgans
        int numMarkers;
        vector <Marker> markers;        // a list of markers for genotyping

        // Selection parameters
        //bool selectionOn;


        // Data and Summary Parameters
        vector<int> summaries;           // the list of generation at which to generate summaries

};

//Experiment constructor
Experiment::Experiment():

{
    //experimentInput();
}

//GET EXPERIMENT PARAMETERS FROM THE USER
void Experiment::experimentInput()
{
    //Ask the user if the data will be entered manually or from a file
    cout << "Do you want to enter the experiment's parameters manually <M>, or from a file <F>?: ";

    enum inputOpt = { M, F};
    inputOpt manFile << cin;

    switch(inputOpt)
    {
        case(M):
        {

        }
    }



}

// RUN THE EXPERIMENT
void Experiment::runExperiment()
{
    // Make a directory to store output files
    const char* dirName = experimentName.c_str();
    mkdir(dirName, 0775);

    // set the migration rates and interferenceOpt so that the rest of the program behaves according to the experiment
    migrationRates = migRates;
    interferenceOpt = intOpt;

    //Begin a loop that starts with the first simulation
    for( sim = 0 ; sim < numSim ; sim++ )
    {
        //cout << "Running Simulation " << sim << endl;
        // Build the initial chromosomes of the population -- also makes the markers
        makeChromosomes();

        // Set up the hybrid zone
        MetaPop* hybZone = new MetaPop(numDemes, demeSizes, ancProportions, demeLocations);

        // iterate through the generations
        runGenerations(hybZone);

        // after all the generations are done, remove all the unused junctions, and print the rest to a file
        cleanJPool();
        cout << generation  << "\t" << junctionPool.size() << endl;
        //textJPool();

        delete hybZone;

        // when the simulation is finished, clear all the junctions
        markers.clear();

        //hybZone.~MetaPop();
        drainJPool();
        cout << generation  << "\t" << junctionPool.size() << endl;
        drainChromosomePool();
    }

}

// JUNCTION POOL FUNCTIONS
// REMOVE EXTINCT JUNCTIONS FROM THE JUNCTION POOL
void Experiment::cleanJPool()
{
    vector<Junction*>::iterator iterJ;
    for (iterJ = junctionPool.begin(); iterJ < junctionPool.end() ; iterJ++ )
    {
        if ( (*(*iterJ)).getNumOccur() <= 0 )   // maybe need to change this to a while loop instead of if statement
        {
            delete *iterJ;
            junctionPool.erase(iterJ);
        }
    }
}

// REMOVE ALL JUNCTIONS FROM MEMORY
void Experiment::drainJPool()
{
    vector<Junction*>::iterator iterJ;
    for (iterJ = junctionPool.begin(); iterJ < junctionPool.end() ; iterJ++ )
    {
        delete (*iterJ);
    }
    junctionPool.clear();
}

// CLEAR CHROMOSOMES FROM THE CHROMOSOME POOL
void Experiment::drainChromosomePool()
{
    vector<Chromosome*>::iterator iterC;
    for ( iterC = chromosomePool.begin() ; iterC < chromosomePool.end() ; iterC++ )
    {
        delete (*iterC);
    }
    chromosomePool.clear();
}

// PRINT ALL THE JUNCTIONS TO A TEXT FILE
void Experiment::textJPool()
{
    errorCheck << "Junction Pool " << sim << " " << generation << endl;
    vector<Junction*>::iterator iterJ;
    for (iterJ = junctionPool.begin(); iterJ < junctionPool.end() ; iterJ++ )
    {
        errorCheck << (*iterJ) << "\t" ;
        (*(*iterJ)).textJunction();
    }
}

// GENOTYPE THE MARKERS USED FOR SUMMARIES OF THE EXPERIMENT
void Experiment::genotypeMarkers(MetaPop* mp)
{
    vector<Marker>::iterator iterM;
    vector<Deme*> mpDemes = (*mp).getDemes();
    vector<Deme*>::iterator iterD;

    for ( iterM = markers.begin() ; iterM < markers.end() ; iterM++ )
    {
        vector<double> aFreqs;
        for ( iterD = mpDemes.begin() ; iterD < mpDemes.end() ; iterD++ )
        {
            double aF = (*iterM).g_Deme( *iterD );
            aFreqs.push_back(aF);
        }

        (*iterM).setDemeFreqs(aFreqs);
        (*iterM).g_Meta(mp);
        (*iterM).gVar_Meta();
    }

}

// TEXT THE MARKERS USED IN THE EXPERIMENT TO A TEXT FILE
void Experiment::textMarkers(string dir, string fname)
{
    vector<Marker>::iterator iter;
    for ( iter = markers.begin() ; iter < markers.end() ; iter++ )
    {
        stringstream mLab;
        mLab << (*iter).getChromosome() << "_" << (*iter).getPosition() << ".";
        string mLabel = mLab.str();
        string file = dir + mLabel + fname;
        const char* clines_File = file.c_str();
        os_clines.open( clines_File, ios_base::app );
        (*iter).textMarker();
        os_clines.close();
    }
}

// MAKE A FILENAME FOR OUTPUT
const char* Experiment::makeFileName( string f, string dir = "" )
{
    stringstream mLab;
    mLab << "_" << generation << "."; // sim <<
    string mLabel = mLab.str();
    string file;
    if ( dir == "" )
    {
        file = mLabel + f;
    }
    else
    {
        file = dir + '/' + mLabel + f;
    }

    const char* fileCC = file.c_str();
    return fileCC;
}

// MAKE THE CHROMOSOMES NEEDED TO START THE POPULATION
vector<Gamete*> Experiment::makeGametePool(int ancestries, vector<double> lengths, vector<ChrType> types)
{
    vector<ChrType>::iterator iterT;
    vector<double>::iterator iterL;

    int needHetGamete = 0;
    ChrType hetGamType;
    double hetGamLength;

    for( iterT = types.begin() , iterL = lengths.begin() ; iterT < types.end(), iterL < lengths.end() ; iterT++, iterL++ )
    {
        if ( *iterT == X )
        {
            needHetGamete = 1;
            hetGamType = Y;
            hetGamLength = *iterL;
        }
        else if ( *iterT == Z )
        {
            needHetGamete = 2;
            hetGamType = W;
            hetGamLength = *iterL;
        }
    }

    if ( needHetGamete == 1 )
    {
        XYgametePool(ancestries, lengths, types);
    }
    else if ( needHetGamete == 2 )
    {
        ZWgametePool(ancestries, lengths, types);
    }
    else
    {
        hermGametePool(ancestries, lengths, types);
    }

    for ( int i = 0 ; i < ancestries ; i++ )
    {
        int n = 0;
        vector<Chromosome*> gameteGenome;
        Chromosome* hetGamChr;

        if( needHetGamete )
        {
            Junction* cen = new Junction(i, 0, j);
            junctionPool.push_back(cen);
            Junction* tel = new Junction(i, chrLengths[i], j);
            junctionPool.push_back(tel);

            CNode *cent1 = new CNode(cen);
            (*cent1).setProxP(0);
            CNode *telo1 = new CNode(tel);
            (*telo1).setDistP(0);

            (*cent1).setDistP(telo1);
            (*telo1).setProxP(cent1);

            hetGamChr = new SexChromosome(cent1, telo1, n, hetGamLength, hetGamType);
        }


        for(iterT = types.begin(), iterL = lengths.begin() ; iterT < types.end() , iterL < lengths.end() ; iterT++, iterL++ )
        {
            Junction* cen = new Junction(i, 0, j);
            junctionPool.push_back(cen);
            Junction* tel = new Junction(i, chrLengths[i], j);
            junctionPool.push_back(tel);

            CNode *cent1 = new CNode(cen);
            (*cent1).setProxP(0);
            CNode *telo1 = new CNode(tel);
            (*telo1).setDistP(0);

            (*cent1).setDistP(telo1);
            (*telo1).setProxP(cent1);

            Chromosome* newChr = new SexChromosome(cent1, telo1, n, *iterL, *iterT);


            Chromosome* chrAncOne = new Chromosome(cent1, telo1, lengths[i]);

            // Destruct the extra copies of CNodes in memory
            (*cent1).~CNode();
            (*telo1).~CNode();

            chromosomePool.push_back(chrAncOne);
            n++;
        }

        makeMarkers(i, chrLengths[i]);
    }
}

// MAKE MARKERS FOR A CHROMOSOME
void Experiment::makeMarkers(int c, double l)
{
    double inc = l / numMarkers;

    if ( l == 0 )
    {
        Marker m(c,l);
        markers.push_back(m);
    }
    else
    {
        for ( double i = 0 ; i <= l ; i+= inc)
        {
            Marker m(c, i);
            markers.push_back(m);
        }
    }
}

// ITERATE THROUGH THE GENERATIONS
void Experiment::runGenerations(MetaPop* mP)
{
    // make sure the vector of summary generations are in order, and set up an iterator for them
    sort( summaries.begin(), summaries.end() );
    vector<int>::iterator iterS = summaries.begin();;

    for ( generation = 0 ; generation < numGens ; generation++ )
    {
        //cout << "New Generation: " << generation << endl;

        (*mP).newGeneration();

        //cout << generation << " " << (*iterS) << endl;

        if ( generation == (*iterS) )
        {
            // a text output file for monitoring clines
            stringstream gen;
            gen << generation;
            string genNum = gen.str();
            string dir = experimentName + "/";
            string clines_FileName = genNum + ".clines.txt";


            //cout << generation  << " " << g << endl;

            genotypeMarkers(mP);
            //textMarkers(dir, clines_FileName);

            const char* junFileName = makeFileName("jpw", dir);
            (*mP).junByWindow(0.001, junFileName);

            cout << (*iterS) << endl;
            iterS++;

        }

        // clean extinct junctions from memory, print the generaton number, and size of the junction pool so the user can see what's going on
        // cout << generation  << "\t" << junctionPool.size() << "\t";
        cleanJPool();
        //cout << junctionPool.size() << endl;
    }
}
#endif // EXPERIMENT_H
