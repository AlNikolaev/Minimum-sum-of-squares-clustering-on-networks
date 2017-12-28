# Minimum-sum-of-squares-clustering-on-networks

main.cpp is the source code of the algorithm which was proposed in Nikolaev, A., Mladenović, N., & Todosijević, R. (2017). J-means and I-means for minimum sum-of-squares clustering on networks. Optimization Letters, 11(2), 359-376.

# Execution

After creation of an executor run it with the following five arguments:
1. The full path to the input file. The information about the format of input file can be found here: http://people.brunel.ac.uk/~mastjjb/jeb/orlib/pmedinfo.html   
2. The time limit (in seconds).
3. Strategy for generating of an initial solution. The value of argument 3 should be 's', 'v', 'e', 'f', 'r'.
4. Strategy for generating new prototypes during working of shaking procedure. The value of argument 4 should be 'v', 'e', 'f', 'r'. 
5. Strategy for removing k-means degeneracy. The value of argument 5 should be 'm', 'v', 'e', 'f', 'r'. 

**According to the computational results the best 3-5 arguments are: v f m**

**Working example for Windows**:
MinSumOfSquares.exe path\pmed1.txt 5 v f m

# Citation

If you use this implementation in your research please cite the following paper:

Nikolaev, A., Mladenović, N. & Todosijević, R. Optim Lett (2017) 11: 359. https://doi.org/10.1007/s11590-015-0974-4
