#include<iostream>
#include<vector>
#include<string>
#include<queue>
#include<utility>
using namespace std;

bool pairMatch(char base1, char base2){
    return (base1 == 'A' && base2 == 'U') ||
           (base1 == 'U' && base2 == 'A') ||
           (base1 == 'C' && base2 == 'G') ||
           (base1 == 'G' && base2 == 'C');
}

vector<vector<int> > findMaxPairs(const string& rna_sequence) {
    int n = rna_sequence.length();
    int edge = n - 5;
    vector<vector<int> > dp(n, vector<int>(n, 0));

    for(int i = 0; i < edge; i++){
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

vector<vector<int> > findSplits(const string& rna_sequence, vector<vector<int> > & numPairs) {
    int n = rna_sequence.length();
    int edge = n - 5;
    vector<vector<int> > dp(n, vector<int>(n, -1));
    int splitIndex = 0;
    for(int i = 0; i < edge; i++){
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
            dp[row][col] = (taken > notTaken) ? splitIndex : dp[row][col-1];
            row++;
            col--;

        }
    }
    return dp;
}

void printMatrix(vector<vector<int> > dp){
    int n = dp.size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            cout << dp[i][j] << " ";
        }
        cout << endl;
    }
}

vector<pair<int, int> > matchingPairs (const string& rna_sequence, vector<vector<int> > &splitIndices){
    int n = rna_sequence.length();
    int edge = n - 5;
    queue<pair<int,int> > q;
    vector<pair<int,int> > matchingIndices;
    int leftInd = edge - 1;
    int rightInd = n - 1;
    q.push(std::make_pair(leftInd, rightInd));
    while(!q.empty()){
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
    return matchingIndices;
}

int main(){
    string rna_sequence = "AAGACUAUACUUUCAGGGAUCAUUUCUAUAGUUUGUUACUAGAGAAGUUUCUCUGAACGUGUAUAGCACUGAAAACCACAAAGAAGAGGUGCAGCAUUAUCUCCUAAGUGUAAAGCCGGCUCUUGAUGUUGCUUUGCUGCAACUGCCAUUUGCCAUUGAUGAUCGUUCUUUUCUUCCUUUGGGAGACUGGGAAGGAAAGGAUGCAAUCUGAGUGG";
    string actualDotBracket = "(((((((((((..(((((.....))))).)))))).))))).......(((.((......)))))..........((((..............((((((((((((..(((.......((((..(((((((((......))))).))))..)))).........)))((..(((((((...))))))).))..)))))))))))).......))))";
    cout << "The given rna sequence is: " << rna_sequence << endl;
    //string rna_sequence = "ACCGGUAGU";
    int rnaLen = rna_sequence.length();
    if(rnaLen <= 5){
        cout << "There are 0 matching pairs. Please enter another sequence." << endl;
        return 1;
    }
    int brackCount = 0;
    for(char c : actualDotBracket){
        if(c == '(') brackCount++;
    }
    cout << "The number of nucleotides in the sequence is: " << rnaLen << endl;
    cout << "The number of pairs in the actual secondary structure is: " << brackCount << endl;
    cout << "The actual dot bracket structure is: " << actualDotBracket << endl;
    vector<vector<int> > numPairs = findMaxPairs(rna_sequence);
    //printMatrix(numPairs);
    vector<vector<int> > splitIndices = findSplits(rna_sequence, numPairs);
    //printMatrix(splitIndices);
    vector<pair<int,int> > pairs = matchingPairs(rna_sequence, splitIndices);
    cout << "The number of optimal pairs obtained for the above RNA sequence is: " << numPairs[rnaLen-6][rnaLen-1] << endl;
    string dotBracket;
    for(int i = 0; i < rnaLen; i++){
        dotBracket.push_back('.');
    }
    for(const pair<int,int> pairInst : pairs){
        int leftBrack = pairInst.first;
        int rightBrack = pairInst.second;
        dotBracket[leftBrack] = '(';
        dotBracket[rightBrack] = ')';
    }
    cout << "The predicted dot bracket structure is: " << dotBracket << endl;
    cout << "Using 1-based indexing, the list of indices at which bases are paired is: " << endl;
    for(const pair<int,int> pairInst : pairs){
        cout << "[" << pairInst.first + 1<< "]" << "," << "[" << pairInst.second + 1 << "]" << endl;
    }
    return 0;
}