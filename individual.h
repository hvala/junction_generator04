
#ifndef INDIVIUDUAL003_H
#define INDIVIUDUAL003_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "chromosome.h"
#include "autosome.h"
#include "sex_chromosome.h"
#include "cytoplasm.h"
#include "junction.h"
#include "gene.h"
#include "landscape.h"

using namespace std;

extern gsl_rng * r;
extern int generation;
extern int interferenceOpt;
extern vector<Junction*> junctionPool;
extern int numChromosomes;
extern Landscape* rFitLand;


class Individual
{
    public:
        // constuctor prototype
        Individual(Gamete* matGamete, Gamete* patGamete, int l, bool m );
        friend class Gene;

        // destructor prototype
        ~Individual();

        // getters
        vector<Chromosome*> getGenome() { return genome; }
        int getSex() { return sex; }
        int getLocation() { return location; }
        double getRFitness() { return phenotypes[0]; }
        double getDFitness() { return phenotypes[1]; }
        double getEFitness() { return phenotypes[2]; }
        double getPhenotype(int i) { return phenotypes[i]; }
        int getLifespan() { return lifespan; }

        // setters
        //void inheritOne(Chromosome *c) { chrOne = c; }
        //void inheritTwo(Chromosome *c) { chrTwo = c; }
        void setLocation(int l) { location = l; }

        // Biological/Genetic Functions/variables
        void determineSex();
        void setRFitness( double r ) { phenotypes[0] = r; }
        void setDFitness( double d ) { phenotypes[1] = d; }
        void setEFitness( double e ) { phenotypes[2] = e; }
        void setPhenotype( int i, double p ) { phenotypes[i] = p; }
        Chromosome* gameteChromosome(Chromosome* c, Chromosome* d);
        Chromosome* gameteHetChromosome(Chromosome* c, Chromosome* d );
        Gamete* makeGamete();
        vector<double> oneCO(double m);
        vector<double> coNoInt(int xo, double l);
        vector<double> coGamInt(double l);
        vector<double> parCrossovers(vector<double> co, double p);
        Chromosome * recombination(Chromosome* c, Chromosome* d, vector<double> coList, ChrType t);
        bool migrated;
        double phenotype(vector<Chromosome*> g, Landscape* l);

        // Summary Statistics Calculators
        int calc_numJunctOne(int x);
        double calc_ATL(int x);
        double calc_HybridIndex();
        double calc_FractionHeterogenic();
        double calc_FractionHomogenic(int a);

        // Data Output/Display
        void displayChromosomes();
        void textChromosomes(int x);

        // Static variables
        static int* gp_type;
        static Gene* genes;

    private:
        vector<Chromosome*> genome;
        vector<double> phenotypes;
        int numChr;
        int sex;
        int lifespan;
        int location;

};

// INDIVIDUAL CONSTRUCTOR
Individual::Individual(Gamete* matGamete, Gamete* patGamete, int l = 0, bool m = false ):
    lifespan(generation), location(l), migrated(m)
{
    vector<Chromosome*> gam1 = (*matGamete).getHapGenome();
    vector<Chromosome*> gam2 = (*patGamete).getHapGenome();
    vector<Chromosome*>::iterator iterG1;
    vector<Chromosome*>::iterator iterG2;

    sex = -1;

    for( iterG1 = gam1.begin(), iterG2 = gam2.begin() ; iterG1 < gam1.end() - 1 ; iterG1++ , iterG2++ )
    {
        ChrType chromoType1 = (**iterG1).getType();

        ChrType chromoType2 = (**iterG2).getType();

        if ( chromoType1 == Y || chromoType2 == Y )
        {
            sex = 1;
        }
        else if ( chromoType1 == W || chromoType2 == W )
        {
            sex = 0;
        }
        else if ( chromoType1 == X && chromoType2 == X )
        {
            sex = 0;
        }
        else if ( chromoType1 == Z && chromoType2 == Z )
        {
            sex = 1;
        }

        genome.push_back(*iterG1);
        genome.push_back(*iterG2);


    }

    ChrType cytoplasmType = (**iterG1).getType();

    switch(cytoplasmType)
    {
        case M:
        {
            genome.push_back(*iterG1);
            break;
        }
        case C:
        {
            genome.push_back(*iterG1);
            break;
        }
        case CP:
        {
            genome.push_back(*iterG2);
            break;
        }

    }

    numChr = genome.size();

    // Determine the phenotypes of the individual
    // for each gene, get the genotype and alter the phenotype accordingly
    // calculate the reproductive fitness
    double rFitness = phenotype(genome, rFitLand);
    phenotypes.push_back(rFitness);

    // calculate the developmental fitness
    // double dFitness = phenotype(genome, dFitLand);
    // phenotypes.push_back(dFitness);

    // calculate the environmental fitness
    // double eFitness = phenotype(genome, eFitLand);
    // phenotypes.push_back(eFitness);



}

// INDIVIDUAL DESTRUCTOR
Individual::~Individual()
{
    vector<Chromosome*>::iterator iterC;
    for ( iterC = genome.begin() ; iterC < genome.end() ; iterC++ )
    {
        delete *iterC;
    }
}

// PRINT THE CHROMOSOMES TO THE TERMINAL
void Individual::displayChromosomes()
{
    int x = genome.size();
    for ( int i = 0 ; i < x ; i++ )
    {
        (*genome[i]).displayChromosome();
    }
}

// PRINT THE CHROMOSOMES TO A TEXT FILE
void Individual::textChromosomes(int x = 2)
{
    for ( int i = 0 ; i < x ; i++ )
    {
        (*genome[i]).textChromosome();
    }
}

// MAKE A GAMETE FROM AN INDIVIDUAL
Gamete* Individual::makeGamete()
{
    vector<Chromosome*> gamete;
    Gamete* newGamete;

    Chromosome* g_Chr;

    ChrType c_type1;
    ChrType c_type2;

    for ( int i = 0 ; i < numChr ; i++ )
    {

        if ( i % 2 != 0 )
        {
            i++;
        }

        if ( i == numChr - 1 )
        {
            gamete.push_back( genome[i] );
        }
        else
        {
            c_type1 = (*genome[i]).getType();
            c_type2 = (*genome[i+1]).getType();

            if( ( c_type1 == A && c_type2 == A ) || ( c_type1 == X && c_type2 == X ) || ( c_type1 == Z && c_type2 == Z ) )
            {
                g_Chr = gameteChromosome(genome[i], genome[i+1] );
                gamete.push_back(g_Chr);
            }
            else if ( ( c_type1 == X && c_type2 == Y ) || ( c_type1 == Y && c_type2 == X ) || ( c_type1 == W && c_type2 == Z ) || ( c_type1 == Z && c_type2 == W )   )
            {
                g_Chr = gameteHetChromosome( genome[i], genome[i+1] );
                gamete.push_back(g_Chr);
            }
            else
            {
                cout << "Mis-Ordered genomic vector. Must be A1,A1,A2,A2,...,An,An,S,S,C" << endl;
                exit(1);
            }
        }
    }

    newGamete = new Gamete(gamete);
    return newGamete;
}

// MAKE A CHROMOSOME FOR THE GAMETE
Chromosome* Individual::gameteChromosome(Chromosome* c, Chromosome* d)
{
    double mu = (*c).getLength();               // the genetic length of the chromosome provides the mean of a Poisson distribution of crossover number
    double strandOpt = gsl_rng_uniform(r);      // strandOpt chooses which of the parental gametes will be chosen, either to be passed on, or as the starting strand for recombination
    Chromosome* gam_Chr;                        // a vector of pointers to the chromosomes that will compose the gamete
    vector<double> crossovers;                  // this vector stores a list of crossovers to be used by the recombination function if needed

    // make the gamete according to the type of recombination being used
    switch(interferenceOpt)
    {
        case 0:
        {
            int x = gsl_rng_uniform(r) + 0.5;
            if(x == 0)
            {
                if ( strandOpt < 0.5 )
                {
                    gam_Chr = (*c).duplicateChr();
                }
                else
                {
                    gam_Chr = (*d).duplicateChr();
                }
            }
            else
            {
                double xo = gsl_rng_uniform(r) * mu;
                crossovers.push_back(xo);
                if ( strandOpt < 0.5 )
                {
                    gam_Chr = recombination(c, d, crossovers, (*d).getType() );
                }
                else
                {
                    gam_Chr = recombination(d, c, crossovers, (*c).getType() );
                }
            }
            break;
        }
        case 1:
        {
            int x = gsl_ran_poisson(r, mu);
            if(x == 0)
            {
                if ( strandOpt < 0.5 )
                {
                    gam_Chr = (*c).duplicateChr();
                }
                else
                {
                    gam_Chr = (*d).duplicateChr();
                }
            }
            else
            {
                crossovers = coNoInt(x, mu);
                if ( strandOpt < 0.5 )
                {
                    gam_Chr = recombination(c, d, crossovers, (*d).getType() );
                }
                else
                {
                    gam_Chr = recombination(d, c, crossovers, (*c).getType() );
                }
            }
            break;
        }
        case 2:
        {
            int x = gsl_rng_uniform(r);
            if (x < 0.5)
            {
                if ( strandOpt < 0.5 )
                {
                    gam_Chr = (*c).duplicateChr();
                }
                else
                {
                    gam_Chr = (*d).duplicateChr();
                }
            }
            else
            {
               crossovers = coGamInt(mu);
                if ( strandOpt < 0.5 )
                {
                    gam_Chr = recombination(c, d, crossovers, (*d).getType() );
                }
                else
                {
                    gam_Chr = recombination(d, c, crossovers, (*c).getType() );
                }
            }
            break;
        }
        default:
        {
            cout << "Recombination interference option not specified correctly (0, 1, or 2)!" << endl;
            exit(1);
            break;
        }
    }

    return gam_Chr;
}


// GET A SEX CHROMOSOME WHEN THE SEX CHROMOSOMES ARE HETEROGAMETIC
Chromosome* Individual::gameteHetChromosome(Chromosome* c, Chromosome* d)
{

    double mu;
    double parB;

    if ( (*c).getType() == Y || (*c).getType() == W )
    {
        mu = (*c).getLength();        // the genetic length of the chromosome provides the mean of a Poisson distribution of crossover number
        parB = (*c).getParB();
    }
    else
    {
        mu = (*d).getLength();
        parB = (*d).getParB();
    }

    double strandOpt = gsl_rng_uniform(r);      // strandOpt chooses which of the parental gametes will be chosen, either to be passed on, or as the starting strand for recombination
    Chromosome* gam_Chr;                        // a vector of pointers to the chromosomes that will compose the gamete
    vector<double> crossovers;                  // this vector stores a list of crossovers to be used by the recombination function if needed

    // allow only one crossover on the par
    int x = gsl_rng_uniform(r) + 0.5;

    if(x == 0)
    {
        if ( strandOpt < 0.5 )
        {
            gam_Chr = (*c).duplicateChr();
        }
        else
        {
            gam_Chr = (*d).duplicateChr();
        }
    }
    else
    {

        double XO = gsl_rng_uniform(r) * mu * parB;
        crossovers.push_back(XO);

        if ( strandOpt < 0.5 )
        {
            gam_Chr = recombination(c, d, crossovers, (*d).getType() );
        }
        else
        {
            gam_Chr = recombination(d, c, crossovers, (*c).getType() );
        }
    }


    return gam_Chr;
}

// RECOMBINE THE TWO CHROMOSOMES OF AN INDIVIDUAL
Chromosome* Individual::recombination(Chromosome* c, Chromosome* d, vector<double> coList, ChrType t)
{
    CNode* curA = 0;     // curA is a pointer to the following CNode on strand that is currently being copied to the recombinant
    CNode* curB = 0;     // curB is the leading CNode
    CNode* oppA = 0;     // the follower on the opposite strand
    CNode* oppB = 0;     // the leader on the opposite strand

    curB = (*c).getCentromere();
    oppB = (*d).getCentromere();



    CNode* recA = new ( (*curB).newAddyAssign() ) CNode( (*curB).getJunction(), location );        // recA is the following CNode pointer on the recombinant, and starts by pointing to the centromere on the current strand

    CNode* recf;

    curB = (*curB).getDistP();                              // advance curB to the next CNode

    Chromosome* recombinant = 0;
    switch(t)
    {
        case A:
        {
            recombinant = new Autosome(recA, 0, (*c).getNumber(), (*c).getLength(), t );      // initialize a recombinant chromosome, and give it the centromere, recA, and leave the telomere as NULL
            break;
        }
        case W:
        {
            recombinant = new SexChromosome(recA, 0, (*c).getNumber(), (*c).getLength(), t, (*c).getParB() );      // initialize a recombinant chromosome, and give it the centromere, recA, and leave the telomere as NULL
            break;
        }
        case X:
        {
            recombinant = new SexChromosome(recA, 0, (*c).getNumber(), (*c).getLength(), t, (*c).getParB() );      // initialize a recombinant chromosome, and give it the centromere, recA, and leave the telomere as NULL
            break;
        }
        case Y:
        {
            recombinant = new SexChromosome(recA, 0, (*c).getNumber(), (*c).getLength(), t, (*c).getParB() );      // initialize a recombinant chromosome, and give it the centromere, recA, and leave the telomere as NULL
            break;
        }
        case Z:
        {
            recombinant = new SexChromosome(recA, 0, (*c).getNumber(), (*c).getLength(), t, (*c).getParB() );      // initialize a recombinant chromosome, and give it the centromere, recA, and leave the telomere as NULL
            break;
        }
        default:
        {
            recombinant = new Autosome(recA, 0, (*c).getNumber(), (*c).getLength(), t);      // initialize a recombinant chromosome, and give it the centromere, recA, and leave the telomere as NULL
            break;
        }
    }

    vector<double>::iterator iterCO;

    for( iterCO = coList.begin() ; iterCO < coList.end() ; iterCO++ )   // for each crossover in the list
    {
        cout << *iterCO << endl;

        while( ( (*curB).getJPosition() <= *iterCO && (*curB).getDistP() != 0 ) )  // until curB goes past the crossover, or reaches the end of the chromosome
        {

            CNode * recB = new ( (*curB).newAddyAssign() ) CNode( (*curB).getJunction(), location );    // make a new CNode for the recombinant, and point it to the junction of curB's CNode
            (*recA).setDistP(recB);                               // make recA point to recB distally
            (*recB).setProxP(recA);                               // make recB point to recA proximally
            curB = (*curB).getDistP();                            // advance curB
            recA = (*recA).getDistP();                            // advance recA
        }

        curA = (*curB).getProxP();                                // advance curA to the node that preceeds curB

        while( (*oppB).getJPosition() <= *iterCO && (*oppB).getDistP() != 0 )   // move oppB past the crossover but not off the chromosome
        {
            oppB = (*oppB).getDistP();
        }

        oppA = (*oppB).getProxP();              // advance oppB to the node preceeding oppB

        int curAnc = (*curA).getJAncestry();    // get the ancestries of the strands at the crossover's position

        int oppAnc = (*oppA).getJAncestry();



        if ( curAnc != oppAnc )                 // if the ancestries differ make a new junction and add it to the strand
        {
            Junction * junct = new  Junction(0, *iterCO, oppAnc);    // make a new junction at the crossovers position, with the opposite strand's ancestry
            junctionPool.push_back(junct);                          // add then new junction to the junction pool
            CNode *recB = new ( (*curB).newAddyAssign() ) CNode(junct, location);                         // make a new CNode that points to the new junction
            (*recA).setDistP(recB);                                 // link recA a to the new Junction distally
            (*recB).setProxP(recA);                                 // link the new junction to recA proximally
            recB = new ( (*curB).newAddyAssign() ) CNode( (*oppB).getJunction(), location );              // make a new node to the next junction which is now on the opposite strand
            recA = (*recA).getDistP();                              // advance recA
            (*recA).setDistP(recB);                                 // link recA and recB as before
            (*recB).setProxP(recA);
            recf = recB;
        }
        else                                                        // otherwise
        {
            CNode *recB = new ( (*curB).newAddyAssign() ) CNode( (*oppB).getJunction(), location );       // skip making and linking a new junciton and crosslink the strands about the new crossover
            (*recA).setDistP(recB);
            (*recB).setProxP(recA);
            recf = recB;
        }

        swap(curA,oppA);            // swap the pointers to reflect the new strand orientation when past the crossover
        swap(curB,oppB);

    }

    while( (*curB).getDistP() != 0 )                        // after the last crossover, finish the strand up until the telomere
    {
        curB = (*curB).getDistP();
        CNode * recB = new ( (*curB).newAddyAssign() ) CNode( (*curB).getJunction(), location);
        recA = (*recA).getDistP();
        (*recA).setDistP(recB);
        (*recB).setProxP(recA);
        recf = recB;
    }

    (*recombinant).setTelomere(recf);

    (*recombinant).displayChromosome();

    return recombinant;
}

// ONE CROSSOVER -- GET ONE RANDOM CROSSOVER
vector<double> Individual::oneCO(double m)
{
    vector<double> co;

    double x = gsl_rng_uniform(r) * m;
    co.push_back(x);

    return co;
}

// POISSON CROSSOVERS -- MAKE A LIST OF CROSSOVERS WITHOUT INTERFERENCE
vector<double> Individual::coNoInt(int c, double m)
{
    vector<double> co;

    for ( int i = 0 ; i < c ; i++)
    {
        double x = gsl_rng_uniform_pos(r) * m;
        co.push_back(x);
    }
    sort(co.begin(), co.end());
    return co;
}

// GAMMA CROSSOVERS -- MAKE A LIST OF CROSSOVERS WITH INTERFERENCE
vector<double> Individual::coGamInt(double m)
{
    vector<double> crossoverList;
    double crossover = m - ( gsl_rng_uniform_pos(r) * m );  // the first position
    crossoverList.push_back(crossover);
    double previousCrossover = crossover;
    double dirChoice = 1 - ( gsl_rng_uniform_pos(r) );      // choose the direction with the greatest oppurtunity for crossing over

    double v = 8;               // The interference parameter, 8 is the mouse genomewide average as determined in Broman et al. 2002.
    double w = 1 / ( 2 * v );    // the rate parameter for the gamma model from above

    do                                                       // put new crossovers that fall on the chromosome onto the list of crossovers
    {
        double gammaRand = gsl_ran_gamma(r, v, w);              // add crossover chosen from gamma distribution  Function: double gsl_ran_gamma (const gsl_rng * r, double a, double b)
        if (dirChoice < 0.5)
        {
            crossover =  gammaRand + previousCrossover;
        }
        else
        {
            crossover =  previousCrossover - gammaRand;
        }

        if ( crossover > 0 && crossover < m )                // if the new crossover is on the chromosome
        {
            crossoverList.push_back(crossover);              // add it to the list
        }
        previousCrossover = crossover;                       // set the new crossover as the previous
    }  while ( crossover > 0 && crossover < m );               // if the new one is off the chromosome, quit adding crossovers

    sort(crossoverList.begin(), crossoverList.end());                // sort the array of positions

    return crossoverList;
}

// REDUCE CROSSOVER POSITIONS TO FIT INTO PSEUDOAUTOSOMAL REGION
vector<double> Individual::parCrossovers(vector<double> co, double p)
{
    vector<double>::iterator iterC;
    for ( iterC = co.begin() ; iterC < co.end() ; iterC++ )
    {
        *iterC = *iterC * p;
    }
    return co;
}

// SUMMARY STATISTICS CALCULATORS

// CALCULATE THE NUMBER OF JUNCTIONS ON ONE CHROMOSOME INDICATED BY X -- 1 = chrOne, 2 = chrTwo, 3 = both chrs
int Individual::calc_numJunctOne(int x)
{
    int numJunct = (*genome[x]).calc_NumJunctions();
    return numJunct;
}

// CALCULATE THE AVERAGE TRACT LENGTH FOR THE GENOME
double Individual::calc_ATL(int x)
{
    double ATL = (*genome[x]).calc_ATL();
    return ATL;
}

// CALCULATE THE HYBRID INDEX OF THE GENOME
double Individual::calc_HybridIndex()
{
    double hybInd = -1;

    return hybInd;
}

// GET THE GENOTYPE KEY FOR A GENE IN AN INDIVIDUAL
double Individual::phenotype(vector<Chromosome*> g, Landscape* l)
{
    vector<int> genotype;
    vector<Gene*> genes = (*l).getLoci();
    vector<Gene*>::iterator iterG;

    for( iterG = genes.begin() ; iterG < genes.end() ; iterG++ )
    {
        int chr =  (**iterG).getChr() * 2;
        int pos =  (**iterG).getPos();

        genotype.push_back( (*g[ chr     ] ).positionAnc(pos) );
        genotype.push_back( (*g[ chr + 1 ] ).positionAnc(pos) );
    }

    string strKey = vec_ints_to_string(genotype);
    double phenotype = (*l).findPhenotype(strKey);

    return phenotype;
}

/*
// CALCULATE THE FRACTION OF THE GENOME THAT IS HETEROGENIC
double Individual::calc_FractionHeterogenic()
{

    double numHetMarkers = 0;

    for ( int i = 0 ; i < nMarkers ; i++ )     // scan along the markers and count the number that are hets
    {

        if ( markers[i].getAncestry() != (*iter2).getAncestry() )
        {
            numHetMarkers++;
        }
        iter2++;
    }
    double numMarkers = (*chrOne).getMarkers().size();
    double fracHeterogenic = numHetMarkers / numMarkers;    // calculate the fraction that are hets
    return fracHeterogenic;
}

// CALCULATE THE FRACTION OF THE GENOME THAT IS HETEROGENIC
double Individual::calc_FractionHomogenic(int a)
{
    vector<Marker>::iterator iter1 = (*chrOne).getMarkers().begin();
    vector<Marker>::iterator iter2 = (*chrTwo).getMarkers().begin();

    double numHetMarkers = 0;

    for ( iter1 ; iter1 < (*chrOne).getMarkers().end() ; iter1++ )
    {
        if ( (*iter1).getAncestry() == a && (*iter2).getAncestry() == a )
        {
            numHetMarkers++;
        }
        iter2++;
    }
    double numMarkers = (*chrOne).getMarkers().size();
    double fracHomogenic = numHetMarkers / numMarkers;
    return fracHomogenic;
}
*/
#endif // INDIVIDUAL003_H
