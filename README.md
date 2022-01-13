# CV Quantum Complex Networks v1

This code is used to generate covariance matrices Gaussian cluster states associated with random graphs of CZ gates. In this particular case, the random graphs are generated with python-iGraph according to the Barabasi-Albert (BA) or Watts-Strogatz (WS) model. The number of nodes in the cluster state and the number of different random realisation can be controled as parameters when the code is called (see BA-Graph-State.sh and WS-Graph-State.sh for examples of scripts that pass the parameters).

For each cluster state, we subtract 0, 1, 2, ..., 10 photons in one particular mode and produce the photon-subtracted state's photon-number correlations <n_i n_j> - <n_i><n_j> between all pairs of modes. This results in a weighted correlation network. This calculation is done by summing over all perfect matchings in a highly optimised way. Because all photons are subtracted in the same mode, many perfect matchings give the same contribution. We go over all different contributions by hand. 

Similarly to the correlation matrix, we also evaluate the covariance matrix of quadrature operators after each photon subtraction.

## Output
For each realisation of a random cluster state, the adjacency matrix of the associated graph is stored as a .txt file. For each realisation of a random cluster state, and each number of subtracted photons, the quadrature covariance matrix and photon-number correlation matrix are stored as .csv files. The code automatically creates directories and a readme file for each run

## Python
This code was written to run with python3.

## Associated research
The code lies at the basis of the numerical results in the following publication: arXiv:2012.15608

## A note on the V2 file
This V2 of the code uses the analytical result that local photon subtraction can only change correlations up to a certain distance from the node in which photons were subtracted. For WS networks with a mean connectivity that is orders of magnitude smaller than the total number of nodes in the network, the number of correlations that are affected by the photon subtraction is much smaller than the total number of nodes. In this version of the code, we first calculate all Gaussian correlations and subsequently only update the correlations that are affected by the photon subtraction.

This V2 version of the code calculates correlations in the initial Gaussian state through an analytical formula that was not used in the previous version of the code.

Both improvements considerable accelerate the calculations, allowing us to study much larger networks. V1 of the code does not make any assumption and calculates Gaussian correlations and correlations in the photon subtracted state in the same way by using pair partitions. In this sense, V1 allows us to test theoretical the results used in numerically.
