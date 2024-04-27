#include <iostream>
#include <fstream> 
#include<vector>
#include<string>
#include <vector>
using namespace std;

void readFromFile(string* rna_sequence,string* actualDotBracket){
    ifstream ip("files/input.txt");

    getline(ip,*rna_sequence);
    getline(ip,*actualDotBracket);

    ip.close();

}

void writeToFile(vector<pair<int,int>> pairs,string rna_sequence){

    FILE* op = fopen("files/output.txt","w");
    if(op==NULL){
        return;
    }
    for(char c: rna_sequence){
         fprintf(op,"%c",c);
    }
   fprintf(op,"\n");
    for(const pair<int,int> pairInst : pairs){
        fprintf(op,"%d %d\n",pairInst.first,pairInst.second);
    }
    fclose(op);
    return;

}