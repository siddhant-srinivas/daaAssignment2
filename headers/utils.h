#include <iostream>
#include <string>
#include <vector>
#include <queue>
using namespace std;

#define line "\n------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"

void printMatrix(vector<vector<int> > dp){   ///prints the 2d table
    int n = dp.size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            cout << dp[i][j] << " ";
        }
        cout << endl;
    }
}


bool pairMatch(char base1, char base2){
    return (base1 == 'A' && base2 == 'U') ||
           (base1 == 'U' && base2 == 'A') ||
           (base1 == 'C' && base2 == 'G') ||
           (base1 == 'G' && base2 == 'C');
}

void printOutputs(string rna_sequence,int rnaLen,string actualDotBracket,int numOptimalPairs,int brackCount,string dotBracket, vector<pair<int,int>> pairs,bool printPairs){   ///function used to print the final outputs of the program
    cout << line << endl;
    cout << "The given rna sequence is [Mus musculus SSU rRNA]: \n" << rna_sequence << endl;
    cout << line << endl;
    cout << "The length of the nucleotides are: "<< rnaLen << endl;
    cout << line << endl;
    cout << "The number of pairs in the actual secondary structure is: " << brackCount << endl;
    cout << "The number of optimal pairs obtained for the above RNA sequence is: " << numOptimalPairs << endl;
    cout << line << endl;
    cout << "The actual dot bracket structure is: \n" << actualDotBracket << endl << endl;
    cout << "The predicted dot bracket structure is: \n" << dotBracket << endl ;
    cout << line << endl;
    if(printPairs){
            cout << "Using 1-based indexing, the list of indices at which bases are paired is: \n" << endl;
    for(const pair<int,int> pairInst : pairs){
        cout << "[" << pairInst.first + 1<< "]" << "," << "[" << pairInst.second + 1 << "]" << endl;
    }
    }

}