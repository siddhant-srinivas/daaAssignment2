#include <iostream>
#include "utils.h"

vector<vector<int>> findMaxPairs(const string& rna_sequence){   ///function that returns the maximum number of base pairs in the rna sequence
    int n = rna_sequence.length();   ///calculating the length of the rna sequence
    int edge = n - 5;  ///sequence length of the base pairs must be atleast 5 hence edge initilaized to n-5
    vector<vector<int>> dp(n, vector<int>(n, 0)); ///dp table to calculate optimum sequence for rna secondary structure prediction

    for(int i = 0; i < edge; i++){  ///filling in dp table that calculates the maximum number of pairs that can be formed by taking or not taking the current pair
        int row = i;
        int col = n - 1;
        while(row < edge && col >= 5){
            int notTaken = dp[row][col-1];
            int taken = -1;
            for(int t = edge - 1 - row; t < col - 4; t++){
                if(pairMatch(rna_sequence[t], rna_sequence[col]) == false){
                    continue;
                }
                int leftSplit = (t > 0) ? dp[row][t-1] : 0;
                int rightSplit = (t < edge - 1) ? dp[edge-t-2][col-1] : 0;
                taken = max(taken, 1 + leftSplit + rightSplit);
            }
            dp[row][col] = max(taken, notTaken);
            row++;
            col--;

        }
    }
    return dp;
}

vector<vector<int>> findSplits(const string& rna_sequence, vector<vector<int> > & numPairs){ ///function that finds the optimal index at which a sequence is split into subproblems
    int n = rna_sequence.length();   ///calculating the length of the rna sequence
    int edge = n - 5;  ///sequence length of the base pairs must be atleast 5 hence edge initilaized to n-5
    vector<vector<int> > dp(n, vector<int>(n, -1));  ///dp table to calculate optimum sequence for rna secondary structure prediction
    int splitIndex = 0;  ///variable to store where a split occurs in the rna molecule initialized to 0
    for(int i = 0; i < edge; i++){  ///iterates over the elements of the dp table
        int row = i;
        int col = n - 1;
        while(row < edge && col >= 5){
            int notTaken = numPairs[row][col-1];
            int taken = -1;
            for(int t = edge - 1 - row; t < col - 4; t++){
                if(pairMatch(rna_sequence[t], rna_sequence[col]) == false){
                    continue;
                }
                int leftSplit = (t > 0) ? numPairs[row][t-1] : 0;
                int rightSplit = (t < edge - 1) ? numPairs[edge-t-2][col-1] : 0;
                if(1 + leftSplit + rightSplit > taken){
                    taken = 1 + leftSplit + rightSplit;
                    splitIndex = t;
                }
            }
            dp[row][col] = (taken > notTaken) ? splitIndex : dp[row][col-1];  ///filling in the dp table
            row++;
            col--;

        }
    }
    return dp;
}



vector<pair<int, int>> matchingPairs (const string& rna_sequence, vector<vector<int> > &splitIndices){   ///function that returns the index of the pairs in the sequence
    int n = rna_sequence.length();  ///calculating the length of the rna sequence
    int edge = n - 5;   /// sequence length of the base pairs must be atleast 5 hence edge initilaized to n-5
    queue<pair<int,int> > q;
    vector<pair<int,int> > matchingIndices;
    int leftInd = edge - 1;  ///setting leftIndex to n-6th index
    int rightInd = n - 1;   ///setting leftIndex to n-1th index
    q.push(std::make_pair(leftInd, rightInd));  ///setting leftIndex and rightIndex as the first pair
    while(!q.empty()){  ///does bfs on the queue that finds the split for each half until there are no splits left
        leftInd = q.front().first;  
        rightInd = q.front().second;
        q.pop();
        if(leftInd < 0 || leftInd > n - 1 || rightInd < 0 || rightInd > n - 1) continue;
        int split = splitIndices[leftInd][rightInd];
        if(split != -1){
            matchingIndices.push_back(std::make_pair(split, rightInd));
            q.push(std::make_pair(leftInd, split - 1));
            q.push(std::make_pair(edge - 2 - split, rightInd - 1));
        }
    }
    return matchingIndices;   ///return index pairs of the sequence
}



string calculateOptimalSequence(string rna_sequence, vector<vector<int>> &numPairs, vector<vector<int>> &splitIndices, vector<pair<int,int>> &pairs){
    numPairs = findMaxPairs(rna_sequence);  ///find the maximum number of base pairs in the rna sequence
    splitIndices = findSplits(rna_sequence, numPairs);  ///finds the optimal index at which a sequence is split into subproblems
    pairs = matchingPairs(rna_sequence, splitIndices);  ///returns the index of the pairs in the sequence
    
    int rnaLen = rna_sequence.length();  ///length of rna sequence
    string dotBracket;
    for(int i = 0; i < rnaLen; i++){  ///initializes dotBracket with dots to represent unpaired bases
        dotBracket.push_back('.');
    }
    for(const pair<int,int> pairInst : pairs){  ///this loop iterates through each pair in the base pair list and accordingly replaces the dot with '(' or ')'
        int leftBrack = pairInst.first;
        int rightBrack = pairInst.second;
        dotBracket[leftBrack] = '(';
        dotBracket[rightBrack] = ')';
    }
    return dotBracket;  ///retruning the final dot bracket sequence
}