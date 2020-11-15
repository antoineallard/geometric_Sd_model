## modelSD

#### Minimal working example

_Each command should be run from the root directory of the project.

Run the python 3 script to generate two sample files for the hidden variables.
```
python3 examples/generate_hidden_variables.py
```

Compile the C++ command line program.
```
g++ -O3 -std=c++11 examples/generate_edgelist_from_modelSD.cpp -o generate_edgelist_from_modelSD
```

Generate an edgelist according to the S^1 model.
```
./generate_edgelist_from_modelSD -n -t -b 2.5 graph01_expo_S1_hidvar.dat  
```

Generate an edgelist according to the S^2 model.
```
./generate_edgelist_from_modelSD -n -d 2 -t -b 2.5 graph01_expo_S2_hidvar.dat  
```
