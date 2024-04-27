## Predicting Optimal 2D Structure of RNA Molecule

Before running the code make sure Python 3 is installed on the system, then run:
```
pip install matplotlib mpl_interactions
```

To run the code on windows:

```
./run.cmd
```

To run the code on MacOS or Linux:

```
g++ main.cpp -o rna && ./rna && python ./python/plotcurve.py
```

To modify the input that the C++ file takes in, go to files/input.txt and modify its contents with the RNA Sequence on the First Line and the Actual Dot Bracket Notation of the RNA Molecule on the Second Line. An Example is shown:
```
GUCUACGGCCAUACCACCCUGAACGCGCCCGAUCUCGUCUGAUCUCGGAAGCUAAGCAGGGUCGGGCCUGGUUAGUACUUGGAUGGGAGACCGCCUGGGAAUACCGGGUGCUGUAGGCUUU
(((((((((....((((((((.....((((((............))))..))....)))))).)).(((((......((.((.(((....))))).)).....))))).)))))))))...
```

To deactivate the live interactive visualisation and just save the visualisation to the file, go to python/plotcurve.py and comment out the line:
```python
plt.show()
```

Upon running the script, the output of the program will be printed onto the console. The time taken by the program to run will also be printed on the console and the visualisation of the structure will be saved as `structure.png` inside the files folder.
