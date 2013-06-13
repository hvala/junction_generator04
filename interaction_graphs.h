#ifndef INTERACTION_GRAPHS_H
#define INTERACTION_GRAPHS_H

#include <map>
#include <algorithm>

class IntGraph
{
    public:
        // Interaction Graph Constructor prototype
        IntGraph( Gene* locA, Gene* locB, int ancA, int ancB, double s, int d );

        // Setters

        // Getters
        double getMaxSel() { return maxSel; }
        Gene* getLocusA() { return loci[0]; }
        Gene* getLocusB() { return loci[1]; }
        int getAncA() { return ancestryA; }
        int getAncB() { return ancestryB; }
        virtual double calc_Selection(vector<int> key, vector<Gene*> loci);

    protected:
        double maxSel;
        vector<Gene*> loci;
        vector<int> ancestries;
        int direction;

        map<string, double> graph;

};

// Iteraction Graph Constructor
IntGraph::IntGraph( Gene* locA, Gene* locB, int ancA, int ancB, double s, int d ):
    maxSel(s), direction(d)
{
    loci.push_back(locA);
    loci.push_back(locB);
    ancestries.push_back(ancA);
    ancestries.push_back(ancB);
}

// Class for a dominant DMI interaction graph
class DomDMI : public IntGraph
{
    public:
        // DomDMI Graph constructor prototype
        DomDMI( Gene* locA, Gene* locB, int ancA, int ancB, double s );

    private:
};

// Dominant Interaction graph constructor
DomDMI::DomDMI( Gene* locA, Gene* locB, int ancA, int ancB, double s, int d ):
    IntGraph( Gene* locA, Gene* locB, int ancA, int ancB, double s, int d )
{

}

// calculate selection for key and genes given
double IntGraph::calc_selection()
{

}

#endif // INTERACTION_GRAPHS_H
