#include <iostream>
#include<chrono>
using namespace std::chrono;
using namespace std;

void computeExecutionTime(string rna_sequence,string actualDotBracket,void (*func)(string,string)){ ///function that computes the execution time that is taken to predict the secondary structure of the rna

    auto start = high_resolution_clock::now();  ///starting the clock

    func(rna_sequence,actualDotBracket);  ///calling the function that predicts the secondary structure of the rna molecule

    auto stop = high_resolution_clock::now();  ///stopping the clock

    auto duration = duration_cast<milliseconds>(stop - start);  ///calculating the execution time of the function that predicts the secondary structure of the rna molecule

    cout << "\nThe time taken for the program to execute is: " << duration.count() << "ms" << endl;
}