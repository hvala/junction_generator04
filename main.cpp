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

    double selCoeff = 0.5;
    //Phenotype normal;
    //Phenotype hybrid;
    //hybrid.rFit -= selCoeff;

    //phenPool.push_back(normal);
    //phenPool.push_back(hybrid);

    generation = 0;


    // Autosome One
    Junction* junOne = new Junction(0, 0, 0);
    CNode* nodeOne = new CNode(junOne);

    Junction* junTwo = new Junction(0, 1, 0);
    CNode* nodeTwo = new CNode(junTwo);

    junctionPool.push_back(junOne);
    junctionPool.push_back(junTwo);

    (*nodeOne).setDistP(nodeTwo);
    (*nodeTwo).setProxP(nodeOne);

    Chromosome* chrOne = new Autosome(nodeOne, nodeTwo, 0, 1.0, A );

    // Autosome Two
    Junction* junThr = new Junction(0, 0, 1);
    CNode* nodeThr = new CNode(junThr);

    Junction* junFour = new Junction(0, 1, 1);
    CNode* nodeFour = new CNode(junFour);

    (*nodeThr).setDistP(nodeFour);
    (*nodeFour).setProxP(nodeThr);

    junctionPool.push_back(junThr);
    junctionPool.push_back(junFour);

    Chromosome* chrTwo = new Autosome(nodeThr, nodeFour, 0, 1.0, A );

    // Sex Chromosome X and Y One
    Junction* sJunOne = new Junction(1, 0, 0);
    CNode* xNodeOne = new CNode(sJunOne);
    CNode* yNodeOne = new CNode(sJunOne);
    junctionPool.push_back(sJunOne);

    Junction* sJunTwo = new Junction(1, 1, 0);
    CNode* xNodeTwo = new CNode(sJunTwo);
    CNode* yNodeTwo = new CNode(sJunTwo);
    junctionPool.push_back(sJunTwo);

    (*xNodeOne).setDistP(xNodeTwo);
    (*xNodeTwo).setProxP(xNodeOne);

    (*yNodeOne).setDistP(yNodeTwo);
    (*yNodeTwo).setProxP(yNodeOne);

    Chromosome* X_One = new SexChromosome(xNodeOne, xNodeTwo, 1, 1.0, X, 0.001);
    Chromosome* Y_One = new SexChromosome(yNodeOne, yNodeTwo, 1, 1.0, Y, 0.001);

    // Sex Chromosome X and Y Two
    Junction* sJunThr = new Junction(1, 0, 1);
    CNode* xNodeThr = new CNode(sJunThr);
    CNode* yNodeThr = new CNode(sJunThr);
    junctionPool.push_back(sJunThr);

    Junction* sJunFour = new Junction(1, 1, 1);
    CNode* xNodeFour = new CNode(sJunFour);
    CNode* yNodeFour = new CNode(sJunFour);

    junctionPool.push_back(sJunFour);

    (*xNodeThr).setDistP(xNodeFour);
    (*xNodeFour).setProxP(xNodeThr);

    (*yNodeThr).setDistP(yNodeFour);
    (*yNodeFour).setProxP(yNodeThr);

    Chromosome* X_Two = new SexChromosome(xNodeThr, xNodeFour, 1, 1.0, X, 0.001);
    Chromosome* Y_Two = new SexChromosome(yNodeThr, yNodeFour, 1, 1.0, Y, 0.001);

    // set up cytoplasm anc 0
    Junction* cJunOne = new Junction(2, 0, 0);
    junctionPool.push_back(cJunOne);
    CNode* cNodeOne = new CNode(cJunOne);

    Junction* cJunTwo = new Junction(2, 0, 0);
    junctionPool.push_back(cJunTwo);
    CNode* cNodeTwo = new CNode(cJunTwo);

    (*cNodeOne).setDistP(cNodeTwo);
    (*cNodeTwo).setProxP(cNodeOne);

    Chromosome* M_One = new Cytoplasm(cNodeOne, cNodeTwo, 2, 0.0, M, 0);

    // set up cytoplasm anc 1
    Junction* cJunThr = new Junction(2, 0, 1);
    junctionPool.push_back(cJunThr);
    CNode* cNodeThr = new CNode(cJunThr);

    Junction* cJunFour = new Junction(2, 0, 1);
    junctionPool.push_back(cJunFour);
    CNode* cNodeFour = new CNode(cJunFour);

    (*cNodeThr).setDistP(cNodeFour);
    (*cNodeFour).setProxP(cNodeThr);

    Chromosome* M_Two = new Cytoplasm(cNodeThr, cNodeFour, 2, 0.0, M, 0);


    vector<Junction*>::iterator iterJ;
/*    for ( iterJ = junctionPool.begin() ; iterJ < junctionPool.end(); iterJ++ )
    {
        (**iterJ).displayJunction();
    }*/

    //(*chrOne).displayChromosome();
    //(*chrTwo).displayChromosome();

    vector<Chromosome*> f1;
    vector<Chromosome*> f2;

    f1.push_back(chrOne); f1.push_back(X_One); f1.push_back(M_One);
    f2.push_back(chrTwo); f2.push_back(X_Two); f2.push_back(M_Two);

    vector<Chromosome*> m1;
    vector<Chromosome*> m2;

    m1.push_back(chrOne); m1.push_back(Y_One); m1.push_back(M_One);
    m2.push_back(chrTwo); m2.push_back(Y_Two); m2.push_back(M_Two);


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

    (*rFitLand).printLandscape();

    Individual* femOne = new Individual(f1, f2);
    cout << "femOne" << endl;
    (*femOne).displayChromosomes( );
    cout << (*femOne).getRFitness() << endl;

    Individual* femTwo = new Individual(f2, f1);
    cout << "femTwo" << endl;
    (*femTwo).displayChromosomes( );
    cout << (*femTwo).getRFitness() << endl;

    Individual* maleOne = new Individual(f1, m2);
    cout << "maleOne" << endl;
    (*maleOne).displayChromosomes( );
    cout << (*maleOne).getRFitness() << endl;

    Individual* maleTwo = new Individual(m1, f2);
    cout << "maleTwo" << endl;
    (*maleTwo).displayChromosomes(  );
    cout << (*maleTwo).getRFitness() << endl;

    vector<Chromosome*> gamete1 = (*femOne).makeGamete();


    vector<Chromosome*>::iterator iterG;
    cout << "G1" << endl;
    for( iterG = gamete1.begin() ; iterG < gamete1.end() ; iterG++ )
    {
        (*(*iterG)).displayChromosome();
    }

    vector<Chromosome*> gamete2 = (*femTwo).makeGamete();
    cout << "G2" << endl;
    for( iterG = gamete2.begin() ; iterG < gamete2.end() ; iterG++ )
    {
        (*(*iterG)).displayChromosome();
    }

    vector<Chromosome*> gamete3 = (*maleOne).makeGamete();
    cout << "G3" << endl;
    for( iterG = gamete3.begin() ; iterG < gamete3.end() ; iterG++ )
    {
        (*(*iterG)).displayChromosome();
    }

    vector<Chromosome*> gamete4 = (*maleTwo).makeGamete();
    cout << "G4" << endl;
    for( iterG = gamete4.begin() ; iterG < gamete4.end() ; iterG++ )
    {
        (*(*iterG)).displayChromosome();
    }

    cout << "IndTen:" << endl;

    Individual* indTen = new Individual(gamete1, gamete3);

    (*indTen).displayChromosomes();
    cout << (*indTen).getRFitness() << endl;

    cout << endl;

    delete indTen;

    CNode* dump = CNode::s_newAddy;
    int d = 0;
    while( dump != 0 )
    {
        cout << d << " " << dump << endl;
        dump = (*dump).getDistP();
        d++;
    }


    cout << endl << "done" << endl;

    int integer = 50;
    vector<int> key = change_base(50, 2, 8);
    string sKey = vec_ints_to_string(key);
    int decKey = convert_to_decimal(key, 2);
    cout << integer << "\t" << sKey << "\t" << decKey << endl;

    return 0;
}
