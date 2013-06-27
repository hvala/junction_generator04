#ifndef GAMETE_H
#define GAMETE_H

class Gamete
{
    public:
        //Gamete constructor prototype
        Gamete(vector<Chromosome*> h);
        ~Gamete();

        //Getters
        Chromosome* getChromosome(int i) { return hapGenome[i]; }
        vector<Chromosome*> getHapGenome() { return hapGenome; }
        //Chromosome* getCytoplasm() { return cyto; }

        //Setters


        // Data display
        void displayGamete();

    protected:
        vector<Chromosome*> hapGenome;


};

// GAMETE CONSTRUCTOR
Gamete::Gamete(vector<Chromosome*> h):
    hapGenome(h)
{

}

// GAMETE DESTRUCTOR
Gamete::~Gamete()
{

}

// DISPLAY A GAMETE IN THE TERMINAL
void Gamete::displayGamete()
{
    vector<Chromosome*>::iterator iterG;
    cout << "Gamete: " << this << endl;

    for( iterG = hapGenome.begin() ; iterG < hapGenome.end() ; iterG++ )
    {
        //cout << (*( (**iterG).getCentromere() )).getJAncestry() << endl;
        (**iterG).displayChromosome();
    }

}

#endif // GAMETE_H
