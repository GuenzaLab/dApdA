#!/usr/bin/env python

import pyemma
import numpy as np
import pyemma.coordinates as coor
import os
import matplotlib.pyplot as plt
import pyemma.msm as msm
import pyemma.plots as mplt
import subprocess
import scipy.interpolate

plt.rcParams.update({'figure.autolayout': True})


print('PYEMMA Version: ',pyemma.__version__)

kBT=0.596   # kcal/mol
DTfile = "DisTwist.dat"

dist, phi = np.loadtxt(DTfile,usecols=(0,1), unpack=True)
data=[np.column_stack([dist,phi])] # Note that the xtc files are saved every 0.2 ps.


nbins=50
hist,x,y = np.histogram2d(data[0][:,0],data[0][:,1], bins=nbins)
sump=0.0
Prob=np.zeros((hist.shape[0],hist.shape[1]))

for i in range(hist.shape[0]):
    for j in range(hist.shape[1]):
        if hist[i,j] == 0:
           Prob[i,j] = 1.0/hist.sum()
        else:
           Prob[i,j] = hist[i,j]/hist.sum()
        sump = sump + Prob[i,j]
print('(Normalization Check) \sigma_i^N P_i :',sump)
FES=-kBT*np.log(Prob.T)
FES=FES-FES.min()

# 1- Clustering
#cl = coor.cluster_uniform_time(data=data, k=100, stride=10)
#cl = coor.cluster_kmeans(data=data, k=250, stride=10)
# for later use we save the discrete trajectories and cluster center coordinates:
#dtrajs = cl.dtrajs
#cc_x = cl.clustercenters[:,0]
#cc_y = cl.clustercenters[:,1]

if os.path.isfile('clusterenters_kmeans.npy'):
    dtrajs=np.load('dtrajs_kmeans.npy')
    dtrajs=np.ravel(dtrajs)
    dummy=np.load('clustercenters_kmeans.npy')
    cc_x=dummy[:,0]
    cc_y=dummy[:,1]
else:
    cl = coor.cluster_kmeans(data=data, k=100, stride=25)
    # for later use we save the discrete trajectories and cluster center pyemma.coordinatesdinates:
    dtrajs = cl.dtrajs
    cc_x = cl.clustercenters[:,0]
    cc_y = cl.clustercenters[:,1]
    np.save('dtrajs_kmeans.npy',dtrajs)
    np.save('clustercenters_kmeans.npy',np.column_stack([cc_x,cc_y]))

# 2- Lag time: Note that the xtc files are saved every 0.2 ps.

# Making the Markov model
M = msm.estimate_markov_model(dtrajs, 2500)
print('fraction of states used = ', M.active_state_fraction)
print('fraction of counts used = ', M.active_count_fraction)
print('transition matrix',M.transition_matrix)  # doctest: +SKIP

# Decomposing the Markov State Models into 5 metastates by PCCA+
nmeta=5
M.pcca(nmeta)

extent = [x[0], x[-1], y[0], y[-1]]
pcca_sets = M.metastable_sets

#Define the nested list holding the macrostate trajectories for each state
macrostate_list=[]
for i in range(nmeta):
        macrostate_list.append([])
        print(macrostate_list)


#Make the macrostate trajectory from the microstate trajectory
metastable_traj=M.metastable_assignments[dtrajs]
for frame,i in enumerate(metastable_traj):
    #add 1 to each frame for now, since Python indexes lists
    #starting at 0, but the index for the first frame in the GROMACS .xtc
    #trajectory is 1.
    macrostate_list[i].append(frame+1)

#Check to make sure I caught all the states
summa=0
for i in range(len(macrostate_list)):
        summa += len(macrostate_list[i])
        print(len(macrostate_list[i]))
print(summa)


#Write frames to file
#Use the subprocess module to remove any existing .ndx files with the same name;
#if it doesn't exist, create one instead.
for i in range(1,nmeta+1):
    if os.path.exists('macrostate_'+str(i)+'.ndx'):
         subprocess.run("rm -rfv macrostate_"+str(i)+".ndx",shell=True)
    else:
         subprocess.run("touch macrostate_"+str(i)+".ndx",shell=True)
    with open("macrostate_"+str(i)+".ndx", 'w+') as f:
         f.write('[ State '+str(i)+' ] \n')
         for k in macrostate_list[i-1]:
             f.write(str(k)+"\n")


#Interpolate on a grid with 100 bins (PyEMMA default)
nbins=100

#Grid the data and interpolate metastable values based on where the trajectory
#is on the surface. NOTE: method='nearest' rounds to the nearest integer, which
#is what we want since we are trying the pass off the memberships as crisp
#assignments to a single metastable state.
xall=dist
yall=phi
x,y=np.meshgrid(np.linspace(xall.min(),xall.max(),nbins),np.linspace(yall.min(),yall.max(),nbins))
z = scipy.interpolate.griddata(np.column_stack([xall,yall]),metastable_traj,(x,y),method='nearest')

#Make the plot: macrostate decomposition with solid borders
ext=[xall.min(),xall.max(),yall.min(),yall.max()]

#Set the min and max values of the contour plot
vmin = np.min(metastable_traj[metastable_traj > -np.inf])
vmax = np.max(metastable_traj[metastable_traj < np.inf])

#fig, ax = plt.subplots(1, 1, figsize=(8, 5))
#ax.contourf(FES,zorder=-3,cmap='Spectral',extent=ext)
#ax.contour(FES,zorder=-2,colors='k',extent=ext)

#Note: feel free to adjust the alpha value. Lower alpha values make the
#overlaid metastable partition more transparent while higher values make it
#more opaque.
ax=plt.gca()
#im1=ax.contourf(x, y, z,levels=None,vmin=vmin,vmax=vmax,alpha=0.7,cmap='Spectral',extent=ext)
#cbar=plt.colorbar(im1,ax=ax)
#cbar.set_label('State')

#Make sure the colorbar is spaced how we want it to be spaced
#cmin, cmax = cbar.get_clim()
#f = (cmax - cmin) / float(nmeta)
#n = np.arange(nmeta)
#cbar.set_ticks((n + 0.5) * f)
#cbar.set_ticklabels(n+1)

#contour plot
fig=plt.contourf(FES, nbins, cmap='gnuplot', extent=ext)
plt.xlabel('R, nm', fontsize=36 )
plt.ylabel('$\phi$, rad', fontsize=36 )
plt.axis([0.3,1.5,-3.14,3.14], fontsize=36 )
plt.xticks(fontsize= 20)
plt.yticks(fontsize= 20)
plt.gcf().subplots_adjust(bottom=0.15)
cbar=plt.colorbar()
cbar.set_label('Free Energy, $k_{\mathrm{B} }T$', fontsize=28)

ax.contour(x, y, z,colors='k',extent=ext)
ax.set_xlim((ext[0],ext[1]))
ax.set_ylim((ext[2],ext[3]))
plt.savefig('MACRO-FES.pdf') #Probably should have a better name
plt.savefig('MACRO-FES.png', format='png',dpi=1000)
plt.savefig('MACRO-FES.eps', format='eps', dpi=1000)

