# coding: utf-8

def gen_protinfo(PROTNAME,G96,TOP):
    import subprocess
    import numpy as np

    protname = str(PROTNAME)
    N = int(subprocess.check_output('grep -c "DA" '+str(TOP),shell=True))
    NFRS = int(subprocess.check_output('grep -c "TIMESTEP" '+str(G96),shell=True))
    NATOMS = int(subprocess.check_output('grep -c "ATOM" '+str(TOP),shell=True))

    array = np.array([protname,N,NFRS,NATOMS],dtype=str)
    np.savetxt("protname.txt",array.T,fmt="%s")

def make_unformatted_traj(G96):
	import numpy as np
	import sys

	txt = np.genfromtxt('protname.txt',dtype=str)
	protname = txt[0]
	N = int(txt[1])
	nfrs = int(txt[2])
	natoms = int(txt[3])

	print(G96)
	print(protname,N,nfrs,natoms)

	dummy = True                    
	linr = []
	with open(G96) as f:
	    for line in f:
	        value = line.split()[0].isalpha()
	        if (dummy == False) and (value == False):
	            linr.append(line_prev.split())
	            linr.append(line.split())
	            for i in range(N-2):
	                linr.append(f.readline().split())
	        dummy = value
	        line_prev = line

	np.save('unformatted_traj.npy',np.array(linr,dtype=float))


gen_protinfo('Atoms-anly','Atoms-anly.g96','top.pdb') 
make_unformatted_traj('Atoms-anly.g96')