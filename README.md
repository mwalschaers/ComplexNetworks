# CV Quantum Complex Networks

This code is used to generate covariance matrices Gaussian cluster states associated with random graphs of CZ gates. In this particular case, the random graphs are generated with python-iGraph according to the Barabasi-Albert (BA) or Watts-Strogatz (WS) model. The number of nodes in the cluster state and the number of different random realisation can be controled as parameters when the code is called (see BA-Graph-State.sh and WS-Graph-State.sh for examples of scripts that pass the parameters).

For each cluster state, we subtract 0, 1, 2, ..., 10 photons in one particular mode and produce the photon-subtracted state's photon-number correlations <n_i n_j> - <n_i><n_j> between all pairs of modes. This results in a weighted correlation network. This calculation is done by summing over all perfect matchings in a highly optimised way. Because all photons are subtracted in the same mode, many perfect matchings give the same contribution. We go over all different contributions by hand. 

Similarly to the correlation matrix, we also evaluate the covariance matrix of quadrature operators after each photon subtraction.

## Output
For each realisation of a random cluster state, the adjacency matrix of the associated graph is stored as a .txt file. For each realisation of a random cluster state, and each number of subtracted photons, the quadrature covariance matrix and photon-number correlation matrix are stored as .csv files. The code automatically creates directories and a readme file for each run

## Python
This code was written to run with python3.

## Associated research
The code lies at the basis of the numerical results in the following publication: arXiv:2012.15608

