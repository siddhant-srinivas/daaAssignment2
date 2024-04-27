#include "headers/fileutils.h"
#include "headers/timeutils.h"
#include "headers/rna.h"
using namespace std;

void predictSecondaryStructure(string rna_sequence,string actualDotBracket){ /// function that predicts the secondary structure of the rna molecule
    int rnaLen = rna_sequence.length();  /// calculating the length of the rna sequence
    if(rnaLen <= 5){  ///if the length of the rna sequence is less than or equal to 5 then there are no matching pairs in that sequence
        cout << "There are 0 matching pairs. Please enter another sequence." << endl;
        return;
    }
    int brackCount = 0;    ///index to count the total number of '(' brackets in the dot bracket sequence
    for(char c : actualDotBracket){
        if(c == '(') brackCount++;   ///increasing the count of brackCount for every open bracket found in the dot bracket sequence
    }
    vector<vector<int>> numPairs, splitIndices;
    vector<pair<int,int>> pairs;
    
    string dotBracket = calculateOptimalSequence(rna_sequence,numPairs,splitIndices,pairs);  ///calling the function that calculates the optimal sequence for the secondary structure of the rna

    printOutputs(rna_sequence, rnaLen, actualDotBracket, numPairs[rnaLen-6][rnaLen-1], brackCount, dotBracket, pairs,false);  ///printing the outputs of the program

    writeToFile(pairs,rna_sequence);
}

int main(){
    string rna_sequence;
    string actualDotBracket;

    readFromFile(&rna_sequence,&actualDotBracket);

    computeExecutionTime(rna_sequence,actualDotBracket,&predictSecondaryStructure); ///calling the function that computes the execution time that is taken to predict the secondary structure of the rna

    
    return 0;
}