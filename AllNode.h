#ifndef ALLNODE_H
#define ALLNODE_H

class AllNode
{
    public:
        AllNode(int s, int* i);
        ~AllNode();

        //


    private:
        AllNode* nextLevel;
        int numStates;
        int* sIndex;

}

// ALLELE NODE CONSTRUCTOR
AllNode::AllNode(int s, int* i = 0):
    numStates(s), sIndex(i)
{
    nextLevel = new AllNode[numStates];
}

#endif //ALLNODE_H
