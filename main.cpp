// Junction Generator -- main.cpp

#include <iostream>
#include <vector>
#include <ctime>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "junction.h"
#include "cnode.h"
#include "chromosome.h"
#include "autosome.h"
#include "sex_chromosome.h"
#include "cytoplasm.h"
#include "individual.h"
#include "gene.h"
#include "jungen_acc.cpp"
#include "landscape.h"
#include "interaction_graphs.h"
#include "deme.h"
#include "metaPop.h"
//#include "marker.h"
//#include "experiment.h"

using namespace std;

vector<Junction*> junctionPool;
int generation;
ofstream errorCheck;
gsl_rng * r;
int interferenceOpt = 1;
int numDemes = 1;
int Junction::numLocations = numDemes;
Landscape* rFitLand;


CNode* junkNode;

int main()
{
    // GSL RNG setup
    r = gsl_rng_alloc(gsl_rng_taus);
    int seed = time(0);
    gsl_rng_set(r, seed);

    // a text output file for checking errors
    string error_FileName = "HybOne/errorCheck.txt";
    const char* error_File = error_FileName.c_str();
    errorCheck.open(error_File);

    generation = 0;

    vector<double> cLens;
    double chrLen[] = {1, 0.5, 0};
    cLens.assign(chrLen, chrLen+3);
    //currentExp.setChrLens(cLens);

    vector<ChrType> cTypes;
    ChrType chrTypes[] = {X, A, M};
    cTypes.assign(chrTypes, chrTypes+3);
    //currentExp.setChrLens(cLens);

    int numAnc = 3;
    vector<Gamete*> gamPool = XYgametePool(numAnc, cLens, cTypes);

    //Gene(int chr, double pos, vector<double> ae);
    vector<double> addEffOne;
    addEffOne.push_back(0);
    addEffOne.push_back(0);
    Gene* geneOne = new Gene(0, 0.1, addEffOne);

    vector<double> addEffTwo;
    addEffTwo.push_back(-0.01);
    addEffTwo.push_back(-0.01);
    Gene* geneTwo = new Gene(0, 0.5, addEffTwo);


    vector<Gene*> genes;
    genes.push_back(geneOne);
    genes.push_back(geneTwo);

    //IntGraph( Gene* locA, Gene* locB, int ancA, int ancB, double s, int d );
    vector<IntGraph*> intGraphs;
    IntGraph* interactionOne = new IntGraph(geneOne, geneTwo, 0, 1, -0.5, 3);
    intGraphs.push_back(interactionOne);

    //Landscape( double a, int na, vector<Gene*> l, vector<IntGraph*> iG );
    rFitLand = new Landscape(1.0, 2, genes, intGraphs);

    //Deme(int s, vector<double> a, int l, vector<Gamete*> gP)
    vector<double> ancProps;
    double ancs[] = {0.6, 0.3, 0.1};
    ancProps.assign(ancs, ancs+numAnc);

    Deme* newDeme = new Deme(100, ancProps, 0, gamPool);

    CNode* dump = CNode::s_newAddy;
    int d = 0;
    while( dump != 0 )
    {
        cout << d << " " << dump << endl;
        dump = (*dump).getDistP();
        d++;
    }

    cout << endl << "done" << endl;

    return 0;
}
