# dApdA
Data and codes for NAR paper on dApdA

codes: contains the codes used to anlyze the data in the manuscript and SI, except to generate the g(r) and solvent orientation around the phosphorous atom in the dApdA backbone, which were generated using the GROMACS analysis functions gmx rdf and gmx sorient, respectively.
- codes/CDcalc_extended_dipole.py: Python code to calculate the CD for each structure in the simulation trajectory using the extended dipole method for the coupling between transition dipole moments on the two adenine bases in dApdA.
- codes/msm.py: Python code to perform the MSM analysis. Requires the PyEMMA library (http://emma-project.org/latest/) be installed.

ndx: Contains a zipped folder with the .ndx files assigning each frame in the dApdA trajectory at 0.1 M to a macrostate in the PCCA+ decomposition of the dApdA free-energy landscape. These are required to obtain the decomposition of the CD shown in Figure 6 of the manuscript.
