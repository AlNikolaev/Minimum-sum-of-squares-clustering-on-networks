# Minimum-sum-of-squares-clustering-on-networks

main.cpp is the source code of the algorithm which was proposed in Nikolaev, A., Mladenović, N., & Todosijević, R. (2017). J-means and I-means for minimum sum-of-squares clustering on networks. Optimization Letters, 11(2), 359-376.

# Execution

After creation of an executor run it with the following five arguments:
1. The full path to the input file. The information about the format of input file can be found here: http://people.brunel.ac.uk/~mastjjb/jeb/orlib/pmedinfo.html   
2. The number of clusters
3. Strategy for finding initial solution ().
4. ...
5. ...

**According to the paper the best arguments are: v f m**

**Working example for Windows**:
MinSumOfSquares.exe path\pmed1.txt 5 v f m

# Citation

If you use this implementation in your research please cite the following paper:

Nikolaev, A., Mladenović, N. & Todosijević, R. Optim Lett (2017) 11: 359. https://doi.org/10.1007/s11590-015-0974-4
