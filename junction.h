
#ifndef JUNCTION003_H
#define JUNCTION003_H

#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;


extern int generation;
extern ofstream errorCheck;
extern int numDemes;


// COMMENTS ON THE JUNCTION CLASS
//
// Background notes:
// The junction is the most fundamental and important datatype in this program.
// A junction is a switchpoint between two different ancestries that occurs along the length
// of a chromosome. They are the result of recombinatin between chromosomal tracts of different
// ancestry, and are inherited as point mutations. As such, they provide a direct measure of the
// amount of hybridization between two diverged lineages of organisms.
//
// Implementation notes:
// The junction is defined by the following simple, numeric parameters:
// 1) Position in the genome
//     a) The chromsome it lies on and b) its genetic (or physical) position on that chromosome
// 2) Ancestry - the ancestry that follows the junction in the direction distal to the centromere
//
// The above three parameters are required to instantiate a junction object.
//
// In addition, we keep track of other parameters
// generation -- the generation the junction was born in
// numOccur -- the number of times the junction occurs in the (meta-)population
//
// Being a fundamental object, it is simple, and few functions beyond getters and setters are needed

// JUNCTION CLASS DEFINITION
class Junction
{
    public:
        // junction constructor prototype
        Junction(int c, double p, int a);

        static int numLocations;

        // setters
        //void setChromosome(int c) { chromosome = c; }
        //void setPosition(double p) { position = p; }
        //void setAncestry(int a) { ancestry = a; }
        //void setGen(int g) { j_gen = g; }

        // getters
        int getChromosome() { return chromosome; }
        double getPosition() { return position; }
        int getAncestry() { return ancestry; }
        int getGen() { return j_gen; }
        int getNumOccur() { return numOccur; }
        int* getLocationCounts() { return locationCounts; }

        // functions
        void incNumOccur() { numOccur++; }
        void decNumOccur() { numOccur--; }
        void incLocCount(int i) { locationCounts[i]++; }
        void decLocCount(int i) { locationCounts[i]--; }
        void displayJunction();
        void textJunction();
        void deallocCN(int i);

        //static vector<Junction*> junctionPool;

    private:
        int chromosome;         // the chromosome on which the junction resides
        double position;        // the genetic position of the junction
        int ancestry;           // the ancestry of the tract following the junction
        int j_gen;              // the generation in which the junction was born
        int numOccur;           // the number of times the junction is used in a simualtion
        int* locationCounts;    // an array containing the number of times the junction occurs in each deme

};



// JUNCTION CONSTRUCTOR
Junction::Junction(int c, double p, int a):
    chromosome(c),
    position(p),
    ancestry(a),
    j_gen(generation)
{
    numOccur = 0;
    locationCounts = new int[numLocations];
    for ( int j = 0 ; j < numLocations ; j++ )
    {
        locationCounts[j] = 0;
    }
}

// DEALLOCATE A JUNCTION FROM A CNODE -- CALLED IN CNODE::DEALLOCATE()
void Junction::deallocCN(int i)
{
    decLocCount(i);
    decNumOccur();
}

// DISPLAY A JUNCTION IN THE TERMINAL
void Junction::displayJunction()
{
    cout << "(" << chromosome << ", " << position << ", " << ancestry << ", " << numOccur << ")" << endl;
}

// DISPLAY A JUNCTION IN THE TERMINAL
void Junction::textJunction()
{
    errorCheck << "\t" << chromosome << "\t" << position << "\t" << ancestry << "\t" << j_gen << "\t" << numOccur << endl; ;
}
#endif // JUNCTION003_H
