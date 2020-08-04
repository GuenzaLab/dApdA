# dApdA
Data and codes for NAR paper on dApdA

codes: contains the codes used to anlyze the data in the manuscript and SI, except to generate the g(r) and solvent orientation around the phosphorous atom in the dApdA backbone, which were generated using the GROMACS analysis functions gmx rdf and gmx sorient, respectively.
- codes/avg_struct_statistics.sh: Utility shell script to get the distance, twist, roll, and tilt coordinates for the dApdA average structures in each of the five macrostates from the MSM analysis.
- codes/gen_unformatted_traj.py: Python script to convert the atomic coordinates of dApdA used to calculate the CD into a an 'unformatted' trajectory stored as a numpy array. The rationale is two-fold: 1) this method should be largely impervious to any changes in PDB and G96 file formats in future versions of GROMACS and 2) loading the numpy array into Python takes significantly less time than a human-readable PDB or G96 file, especially if the trajectory is large.
- codes/CDcalc_extended_dipole.py: Python code to calculate the CD for each structure in the simulation trajectory using the extended dipole method for the coupling between transition dipole moments on the two adenine bases in dApdA.
- codes/msm.py: Python code to perform the MSM analysis. Requires the PyEMMA library (http://emma-project.org/latest/) be installed.
- codes/param_calc.f08: FORTRAN code to generate the distance, twist, roll, and tilt coordinates of dApdA for the average structures inside each macrostate. 
- codes/run_python.pbs: SLURM submission script (but can also be run using 'bash' or 'sh' in the shell) to run the CD analysis. This script is specifically formatted to run on the Comet supercomputer (https://www.sdsc.edu/support/user_guides/comet.html); if run on a local machine, two changes should be made. First, the 'module load gromacs' line should be commmented out. Second, one must ensure that GROMACS (http://manual.gromacs.org/) is installed. If GROMACS is installed without MPI parallelization, the GROMACS executable ('gmx_mpi' in the script) should be changed to 'gmx'.

codes/notebooks: Jupyter notebooks used in the MSM, CD, and solvent structure and orientation analyses. 
ndx: Contains a zipped folder with the .ndx files assigning each frame in the dApdA trajectory at 0.1 M to a macrostate in the PCCA+ decomposition of the dApdA free-energy landscape. These are required to obtain the decomposition of the CD shown in Figure 6 of the manuscript.
