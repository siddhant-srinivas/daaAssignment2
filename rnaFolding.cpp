#include<iostream>
#include<vector>
#include<string>
#include<queue>
#include<utility>
#include<chrono>
using namespace std::chrono;
using namespace std;

#define line "\n------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"

bool pairMatch(char base1, char base2){
    return (base1 == 'A' && base2 == 'U') ||
           (base1 == 'U' && base2 == 'A') ||
           (base1 == 'C' && base2 == 'G') ||
           (base1 == 'G' && base2 == 'C');
}

vector<vector<int> > findMaxPairs(const string& rna_sequence){   ///function that returns the maximum number of base pairs in the rna sequence
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

vector<vector<int> > findSplits(const string& rna_sequence, vector<vector<int> > & numPairs){ ///function that finds the optimal index at which a sequence is split into subproblems
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

void printMatrix(vector<vector<int> > dp){   ///prints the 2d table
    int n = dp.size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            cout << dp[i][j] << " ";
        }
        cout << endl;
    }
}

vector<pair<int, int> > matchingPairs (const string& rna_sequence, vector<vector<int> > &splitIndices){   ///function that returns the index of the pairs in the sequence
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

void printOutputs(string rna_sequence,int rnaLen,string actualDotBracket,int numOptimalPairs,int brackCount,string dotBracket, vector<pair<int,int>> pairs){   ///function used to print the final outputs of the program
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
    // cout << "Using 1-based indexing, the list of indices at which bases are paired is: \n" << endl;
    // for(const pair<int,int> pairInst : pairs){
    //     cout << "[" << pairInst.first + 1<< "]" << "," << "[" << pairInst.second + 1 << "]" << endl;
    // }
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

    printOutputs(rna_sequence, rnaLen, actualDotBracket, numPairs[rnaLen-6][rnaLen-1], brackCount, dotBracket, pairs);  ///printing the outputs of the program
}

void computeExecutionTime(string rna_sequence,string actualDotBracket){ ///function that computes the execution time that is taken to predict the secondary structure of the rna

    auto start = high_resolution_clock::now();  ///starting the clock

    predictSecondaryStructure(rna_sequence,actualDotBracket);  ///calling the function that predicts the secondary structure of the rna molecule

    auto stop = high_resolution_clock::now();  ///stopping the clock

    auto duration = duration_cast<milliseconds>(stop - start);  ///calculating the execution time of the function that predicts the secondary structure of the rna molecule

    cout << "\nThe time taken for the program to execute is: " << duration.count() << "ms" << endl;
}
int main(){
    string rna_sequence = "UACCUGGUUGAUCCUGCCAGUAGCAUAUGCUUGUCUCAAAGAUUAAGCCAUGCAUGUCUAAGUACGCACGGCCGGUACAGUGAAACUGCGAAUGGCUCAUUAAAUCAGUUAUGGUUCCUUUGGUCGCUCGCUCCUCUCCUACUUGGAUAACUGUGGUAAUUCUAGAGCUAAUACAUGCCGACGGGCGCUGACCCCCCUUCCCGGGGGGGGAUGCGUGCAUUUAUCAGAUCAAAACCAACCCGGUGAGCUCCCUCCCGGCUCCGGCCGGGGGUCGGGCGCCGGCGGCUUUGGUGACUCUAGAUAACCUCGGGCCGAUCGCACGCCCCCCGUGGCGGCGACGACCCAUUCGAACGUCUGCCCUAUCAACUUUCGAUGGUAGUCGCCGUGCCUACCAUGGUGACCACGGGUGACGGGGAAUCAGGGUUCGAUUCCGGAGAGGGAGCCUGAGAAACGGCUACCACAUCCAAGGAAGGCAGCAGGCGCGCAAAUUACCCACUCCCGACCCGGGGAGGUAGUGACGAAAAAUAACAAUACAGGACUCUUUCGAGGCCCUGUAAUUGGAAUGAGUCCACUUUAAAUCCUUUAACGAGGAUCCAUUGGAGGGCAAGUCUGGUGCCAGCAGCCGCGGUAAUUCCAGCUCCAAUAGCGUAUAUUAAAGUUGCUGCAGUUAAAAAGCUCGUAGUUGGAUCUUGGGAGCGGGCGGGCGGUCCGCCGCGAGGCGAGUCACCGCCCGUCCCCGCCCCUUGCCUCUCGGCGCCCCCUCGAUGCUCUUAGCUGAGUGUCCCGCGGGGCCCGAAGCGUUUACUUUGAAAAAAUUAGAGUGUUCAAAGCAGGCCCGAGCCGCCUGGAUACCGCAGCUAGGAAUAAUGGAAUAGGACCGCGGUUCUAUUUUGUUGGUUUUCGGAACUGAGGCCAUGAUUAAGAGGGACGGCCGGGGGCAUUCGUAUUGCGCCGCUAGAGGUGAAAUUCUUGGACCGGCGCAAGACGGACCAGAGCGAAAGCAUUUGCCAAGAAUGUUUUCAUUAAUCAAGAACGAAAGUCGGAGGUUCGAAGACGAUCAGAUACCGUCGUAGUUCCGACCAUAAACGAUGCCGACUGGCGAUGCGGCGGCGUUAUUCCCAUGACCCGCCGGGCAGCUUCCGGGAAACCAAAGUCUUUGGGUUCCGGGGGGAGUAUGGUUGCAAAGCUGAAACUUAAAGGAAUUGACGGAAGGGCACCACCAGGAGUGGAGCCUGCGGCUUAAUUUGACUCAACACGGGAAACCUCACCCGGCCCGGACACGGACAGGAUUGACAGAUUGAUAGCUCUUUCUCGAUUCCGUGGGUGGUGGUGCAUGGCCGUUCUUAGUUGGUGGAGCGAUUUGUCUGGUUAAUUCCGAUAACGAACGAGACUCUGGCAUGCUAACUAGUUACGCGACCCCCGAGCGGUCGGCGUCCCCCAACUUCUUAGAGGGACAAGUGGCGUUCAGCCACCCGAGAUUGAGCAAUAACAGGUCUGUGAUGCCCUUAGAUGUCCGGGGCUGCACGCGCGCUACACUGACUGGCUCAGCGUGUGCCUACCCUACGCCGGCAGGCGCGGGUAACCCGUUGAACCCCAUUCGUGAUGGGGAUCGGGGAUUGCAAUUAUUCCCCAUGAACGAGGAAUUCCCAGUAAGUGCGGGUCAUAAGCUUGCGUUGAUUAAGUCCCUGCCCUUUGUACACACCGCCCGUCGCUACUACCGAUUGGAUGGUUUAGUGAGGCCCUCGGAUCGGCCCCGCCGGGGUCGGCCCACGGCCCUGGCGGAGCGCUGAGAAGACGGUCGAACUUGACUAUCUAGAGGAAGUAAAAGUCGUAACAAGGUUUCCGUAGGUGAACCUGCGGAAGGAUCAUUA";
    string actualDotBracket = "...((((.........))))(((.((((((...(.((..((.....(((.(((..((....((....((..........))...)).))......(((.......(((.(.(..((....(((....................((.....((.(((.....))).))......)).........((((...((((((......))))))...))))((..(((((..............(((.((.(((((((((((((....))))))))).))))...)))))....(......)..)))))......)).(.((((.((((......)))))))))...)))...))).).))).(((....(((....(((((((.........))))))))))......)))...(((.((((....))))....)))))).((.(((..........))).)).(.((....)).)...)))))).).....(.(((....(((.....)))..)))).)...)).).....(((.(((((.(((....))).).))))..)))......((..(...........)..)).........(((((........((((.....((....)).......)))).)))))..)))))).))).........(.(.......(.((...(((.((.....((.(.((((((((((((((((....)))).)..))))))))))).).))....((.(..(((.(((((.(.(((((((......)))))))..).))))))))).))..((((((.(.......((..((.......))((((.......))))...))......).)))..))).......((.(..(.((.((...............)).)).)..)))....))...(((((((..(...((((..(((.((((((((...((........))......)))))))).))).......((....))...))))..)..))).)))..))))...)).)....((((((....(...((((.........))))...).))))))..........(((.((.(((..(.((((((.(((((....))))))))))).)..)))...((....))...)).....))).).).(((......";
    
    
    computeExecutionTime(rna_sequence,actualDotBracket); ///calling the function that computes the execution time that is taken to predict the secondary structure of the rna

    
    return 0;
}