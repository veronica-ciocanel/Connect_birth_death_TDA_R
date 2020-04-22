# Connect_birth_death_TDA_R
Sample R code for connecting and visualizing paths of birth-death pairs through time in persistence diagrams generated from time-series point cloud data.

This code accompanies a preprint on "Topological data analysis approaches to uncovering the timing of ring structure onset in filamentous networks" by Ciocanel, Juenemann, Dawes, and McKinley.

We provide a sample dataset for the x, y, and z locations of actin monomer beads for each of 200 time frames in a MEDYAN simulation with standard parameters (see reference [1]). This dataset corresponds to extracting roughly 10% of the beads along each simulated actin filament to generate the point cloud at each time.

The main code connect_birth_death.R loads the sample dataset, generates the birth-death pairs in the persistence diagram for the point cloud at each time frame, and connects these pairs through time. It also extracts the most significant path corresponding to a 1- or 2-dimensional hole (as specified by the user with parameter maxBetti) in the simulation data. The code provides visualizations of all paths connected through time in persistence diagram space, of the most significant path connected through time in persistence diagram space, and of the persistence of the most significant path as a function of time. The supporting functions are included in the main code. 

The persistence diagrams are generated using the ripsDiag function in the TDA package in R, which calculates the Vietoris-Rips filtration built on top of a point cloud (see reference [2]). In particular, we use this function with the GUDHI C++ library for computing persistence diagrams (see reference [3]).

Sample Python code for the same method of connecting paths of birth-death pairs (for 1 or 2-dimensional holes) in the persistence diagrams generated from time-series data (represented as point clouds naturally or extracted from filamentous networks) is also provided in a Jupyter Notebook. Here, the persistence diagrams are generated using the Ripser.py library in Python (see reference [4]), which uses the Ripser C++ library (see reference [5]).

[1] Popov K, Komianos J, Papoian GA (2016). MEDYAN: Mechanochemical Simulations of
Contraction and Polarity Alignment in Actomyosin Networks. PLoS Comput Biol 12(4): e1004877. doi:10.1371/journal.pcbi.1004877
[2] Fasy BT, Kim J, Lecci F, Maria C (2014a). Introduction to the R package TDA. arXiv preprint arXiv:14111830
[3] Maria C, Boissonnat JD, Glisse M, Yvinec M (2014). The GUDHI library: Simplicial  complexes and persistent homology. In: International Congress on Mathematical Software, Springer, pp 167–174
[4] Tralie C, Saul N, Bar-On R. Ripser.py (2018): A lean persistent homology library for python. Journal of Open Source Software. 2018 Sep 13;3(29):925. Software available at https://github.com/scikit-tda/Ripser.py
[5] Bauer U. Ripser (2014). A lean C++ code for the computation of Vietoris–Rips persistence barcodes. Software available at https://github.com/Ripser/ripser 
