#ifndef INTERACTION_GRAPHS_H
#define INTERACTION_GRAPHS_H

#include <map>
#include <algorithm>
#include <string>
#include <sstream>

#include "jungen_acc.cpp"

class IntGraph
{
    public:
        // Interaction Graph Constructor prototype
        IntGraph( Gene* locA, Gene* locB, int ancA, int ancB, double s, int d );

        // Setters

        // Getters
        double getMaxSel() { return maxSel; }
        Gene* getLocusA() { return locusA; }
        Gene* getLocusB() { return locusB; }
        int getAncA() { return ancA; }
        int getAncB() { return ancB; }
        double calc_Selection(int a1, int a2, int b1, int b2);


    protected:
        double maxSel;
        Gene* locusA;
        Gene* locusB;
        int ancA;
        int ancB;
        int edgeFX;
        map<string, int> graph;
};

// Iteraction Graph Constructor
IntGraph::IntGraph( Gene* locA, Gene* locB, int a, int b, double s, int d ):
    locusA(locA), locusB(locB), ancA(a), ancB(b), maxSel(s), edgeFX(d)
{

}


// calculate selection for key and genes given
double IntGraph::calc_Selection(int a1, int a2, int b1, int b2)
{
    int edges = 0;
    double value = 0;

    if( a1 == ancA && b1 == ancB )  { edges++; }
    if( a1 == ancA && b2 == ancB )  { edges++; }
    if( a2 == ancA && b1 == ancB )  { edges++; }
    if( a2 == ancA && b2 == ancB )  { edges++; }

    //cout << a1 << "\t" << a2 << "\t" << b1 << "\t" << b2 << "\t\t" << edges << endl;

    switch(edgeFX)
    {
        case(0):
            value = recEdges(edges, maxSel);
            break;

        case(1):
            value = domEdges(edges, maxSel);
            break;

        case(2):
            if ( edges > 1 && b1 == b2 )
            {
                value = maxSel;
            }
            break;

        case(3):
            value = addEdges(edges, maxSel);
            break;

        default:
            value = 0;
            cout << "Edge Effect not properly specified!" << endl;
            break;
    }

    return value;
}




#endif // INTERACTION_GRAPHS_H
