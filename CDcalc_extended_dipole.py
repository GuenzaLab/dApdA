# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import sympy as sp

print(os.getcwd())

def J2wn(J):
    #convert cm-1 to J(m2*Kg*s-2)
    C0=2.998E8; #speed of light m/s
    h=6.626E-34; #plank constant
    return J/(h*C0*100)

def wn2J(V):
    C0=2.998E8; #speed of light m/s
    h=6.626E-34; #plank constant
    return h*C0*100*V;

def wntofre(x):
    #convert wavenumber(cm-1) to frequency(s-1)
    C0=2.998*1.0E8; #m*s-1 speed of light
    #?:wavelength(nm)
    #?:frequency(s-1)
    #3k:wavenumber(cm-1)
    #?*10^(-9)=C0/?; k=10^7/?
    #?=10^9*C0/?; ?=10^9*C0/?
    return 1.0E2*C0*x

def fretonm(x):
    #convert wavelength(nm) to frequency(s-1)
    C0=2.998*1.0E8; #m*s-1 speed of light
    #?:wavelength(nm)
    #?:frequency(s-1)
    #k:wavenumber(cm-1)
    #?*10^(-9)=C0/?; k=10^7/?
    #?=10^9*C0/?; ?=10^9*C0/?
    return (1.0E9*C0)/x;

def nmtown(x):
    return 1.0E7/x;

def f2debye(f,vo):
    #f is dimensionless osscilation strength
    #vo is transition center in wavenumber
    me=9.10938356E-31; #Kg
    e=-1.602E-19; #C
    c=2.9979E8; #m/s
    D=3.336000000000000E-30;
    eo=8.854E-12; #kg-1m-3s2C2
    NA=6.022E23;
    vo=wntofre(vo);
    pre=4*me*c*eo*np.log(10)/NA/e**2/10;
    #pre=4.39E-9;
    return np.sqrt(1.022E-61/vo*(f/pre))/D;

def debye2SI(D):
    #D dipole in debye
    return D*3.336E-30; #C*m

#Use the strange sign convention

def calcRk_sign(eigvec,R12,Emat,Coupmat,mumat):
    hbar=1.054E-34; #J s hbar
    e0=8.854E-12; #permittivity of a vacuum Kg-1m-3s2C2
    diel=2.0;
    R12=R12/(10**10);  #Ang

    Rhat = R12/np.linalg.norm(R12)

    nbasis = len(Emat)
    
    n1 = nbasis//2
    n2 = nbasis - n1
    Rk = 0
    degen = 0
    nondegen = 0
    crossmat = np.zeros((nbasis,nbasis))
    for i in range(nbasis):
        for j in range(i,nbasis):
            crossmat[i,j] = np.dot(np.cross(mumat[j,:],mumat[i,:]),R12)
            crossmat[j,i] = -crossmat[i,j]
    
    for i in range(1,nbasis+1):
        for j in range(1,nbasis+1):
            # i indexes the A chromophore, j the C chromophore
            if (i <= n1) and (j > n1):
                #Non-degenerate case
                if Emat[i-1] != Emat[j-1]:
                    V = Coupmat[i-1,j-1]
                    pre1ac = -Emat[i-1]*Emat[j-1]/(hbar*(Emat[i-1]**2-Emat[j-1]**2))
                    pre3ac = np.dot(np.cross(mumat[i-1,:],mumat[j-1,:]),R12)
                    
                    if Emat[i-1] > Emat[j-1]:
                        Rk += eigvec[i-1]*eigvec[j-1]*(pre1ac*V*pre3ac)
                        nondegen += eigvec[i-1]*eigvec[j-1]*(pre1ac*V*pre3ac)
                    else:
                        Rk += eigvec[i-1]*eigvec[j-1]*(-pre1ac*V*pre3ac)
                        nondegen += eigvec[i-1]*eigvec[j-1]*(-pre1ac*V*pre3ac)
                #Degenerate case
                else:
                    if (i + n1 == j) or (j - n1 == i):
                        #Factor of 1/2\hbar b/c we are not counting the degenerate terms for
                        #the C-A coupling later.
                        pre1 = Emat[i-1]/(2*hbar)
                        pre2 = np.dot(np.cross(mumat[j-1,:],mumat[i-1,:]),R12)
                        Rk += eigvec[i-1]*eigvec[j-1]*pre1*pre2
                        degen += eigvec[i-1]*eigvec[j-1]*pre1*pre2
                    else:
                        #print("I actually came here.")
                        V = Coupmat[i-1,j-1]
                        pre1ac = -Emat[i-1]*Emat[j-1]/(hbar*(Emat[i-1]**2-Emat[j-1]**2))
                        pre3ac = np.dot(np.cross(mumat[i-1,:],mumat[j-1,:]),R12)

                        if Emat[i-1] > Emat[j-1]:
                            Rk += eigvec[i-1]*eigvec[j-1]*(pre1ac*V*pre3ac)
                            nondegen += eigvec[i-1]*eigvec[j-1]*(pre1ac*V*pre3ac)
                        else:
                            Rk += eigvec[i-1]*eigvec[j-1]*(-pre1ac*V*pre3ac)
                            nondegen += eigvec[i-1]*eigvec[j-1]*(-pre1ac*V*pre3ac)
            # i indexes the C chromophore, j the A chromophore
            if (i > n1) and (j <= n1):
                if Emat[i-1] != Emat[j-1]:
                    V = Coupmat[i-1,j-1]
                    pre1ac = -Emat[i-1]*Emat[j-1]/(hbar*(Emat[i-1]**2-Emat[j-1]**2))
                    pre3ac = np.dot(np.cross(mumat[i-1,:],mumat[j-1,:]),R12)
                    
                    if Emat[i-1] > Emat[j-1]:
                        Rk += eigvec[i-1]*eigvec[j-1]*(-pre1ac*V*pre3ac)
                        nondegen += eigvec[i-1]*eigvec[j-1]*(-pre1ac*V*pre3ac)
                    else:
                        Rk += eigvec[i-1]*eigvec[j-1]*(pre1ac*V*pre3ac)
                        nondegen += eigvec[i-1]*eigvec[j-1]*(pre1ac*V*pre3ac)
                    
                #Degenerate case
                #These terms should NEVER be evaluated. Not sure why they are in the original code.
                else:
                    #These degenerate terms have already been calculated above
                    if (i - n1 == j) or (j + n1 == i):
                        pass
                    else:
                        print(i,j,i+n1,j-n1)
                        #print("I actually came here.")
                        V = Coupmat[i-1,j-1]
                        pre1ac = -Emat[i-1]*Emat[j-1]/(hbar*(Emat[i-1]**2-Emat[j-1]**2))
                        pre3ac = np.dot(np.cross(mumat[i-1,:],mumat[j-1,:]),R12)
                        #print(pre1ac,pre3ac)
                        if Emat[i-1] > Emat[j-1]:
                            Rk += eigvec[i-1]*eigvec[j-1]*(-pre1ac*V*pre3ac)
                            nondegen += eigvec[i-1]*eigvec[j-1]*(-pre1ac*V*pre3ac)
                        else:
                            Rk += eigvec[i-1]*eigvec[j-1]*(pre1ac*V*pre3ac)
                            nondegen += eigvec[i-1]*eigvec[j-1]*(pre1ac*V*pre3ac)
    return Rk,degen,nondegen

def calcRk_degen(eigvec,R12,Emat,Coupmat,mumat):
    hbar=1.054E-34; #J s hbar
    e0=8.854E-12; #permittivity of a vacuum Kg-1m-3s2C2
    diel=2.0;
    R12=R12/(10**10);  #Ang

    Rhat = R12/np.linalg.norm(R12)

    nbasis = len(Emat)
    
    n1 = nbasis//2
    n2 = nbasis - n1
    Rk = 0
    degen = 0
    nondegen = 0
    crossmat = np.zeros((nbasis,nbasis))
    for i in range(nbasis):
        for j in range(i,nbasis):
            crossmat[i,j] = np.dot(np.cross(mumat[j,:],mumat[i,:]),R12)
            crossmat[j,i] = -crossmat[i,j]
    
    for i in range(1,nbasis+1):
        for j in range(1,nbasis+1):
            # i indexes the A chromophore, j the C chromophore
            if (i <= n1) and (j > n1):
                #Non-degenerate case
                if Emat[i-1] != Emat[j-1]:
                    pass
                    #V = Coupmat[i-1,j-1]
                    #pre1ac = -Emat[i-1]*Emat[j-1]/(hbar*(Emat[i-1]**2-Emat[j-1]**2))
                    #pre3ac = np.dot(np.cross(mumat[i-1,:],mumat[j-1,:]),R12)
                    
                    #if Emat[i-1] > Emat[j-1]:
                    #    Rk += eigvec[i-1]*eigvec[j-1]*(pre1ac*V*pre3ac)
                    #    nondegen += eigvec[i-1]*eigvec[j-1]*(pre1ac*V*pre3ac)
                    #else:
                    #    Rk += eigvec[i-1]*eigvec[j-1]*(-pre1ac*V*pre3ac)
                    #    nondegen += eigvec[i-1]*eigvec[j-1]*(-pre1ac*V*pre3ac)
                #Degenerate case
                else:
                    if (i + n1 == j) or (j - n1 == i):
                        #Factor of 1/2\hbar b/c we are not counting the degenerate terms for
                        #the C-A coupling later.
                        pre1 = Emat[i-1]/(2*hbar)
                        pre2 = np.dot(np.cross(mumat[j-1,:],mumat[i-1,:]),R12)
                        Rk += eigvec[i-1]*eigvec[j-1]*pre1*pre2
                        degen += eigvec[i-1]*eigvec[j-1]*pre1*pre2
                    else:
                        pass
                        #print("I actually came here.")
                        #V = Coupmat[i-1,j-1]
                        #pre1ac = -Emat[i-1]*Emat[j-1]/(hbar*(Emat[i-1]**2-Emat[j-1]**2))
                        #pre3ac = np.dot(np.cross(mumat[i-1,:],mumat[j-1,:]),R12)

                        #if Emat[i-1] > Emat[j-1]:
                        #    Rk += eigvec[i-1]*eigvec[j-1]*(pre1ac*V*pre3ac)
                        #    nondegen += eigvec[i-1]*eigvec[j-1]*(pre1ac*V*pre3ac)
                        #else:
                        #    Rk += eigvec[i-1]*eigvec[j-1]*(-pre1ac*V*pre3ac)
                        #    nondegen += eigvec[i-1]*eigvec[j-1]*(-pre1ac*V*pre3ac)
            # i indexes the C chromophore, j the A chromophore
            if (i > n1) and (j <= n1):
                if Emat[i-1] != Emat[j-1]:
                    pass
                    #V = Coupmat[i-1,j-1]
                    #pre1ac = -Emat[i-1]*Emat[j-1]/(hbar*(Emat[i-1]**2-Emat[j-1]**2))
                    #pre3ac = np.dot(np.cross(mumat[i-1,:],mumat[j-1,:]),R12)
                    
                    #if Emat[i-1] > Emat[j-1]:
                    #    Rk += eigvec[i-1]*eigvec[j-1]*(-pre1ac*V*pre3ac)
                    #    nondegen += eigvec[i-1]*eigvec[j-1]*(-pre1ac*V*pre3ac)
                    #else:
                    #    Rk += eigvec[i-1]*eigvec[j-1]*(pre1ac*V*pre3ac)
                    #    nondegen += eigvec[i-1]*eigvec[j-1]*(pre1ac*V*pre3ac)
                    
    return Rk,degen,nondegen

#Parameters from Nordén's JACS paper (Table 4)

angleA=[66,19,-15,-21,-64];  # Degrees
miuA=[f2debye(0.047,36710),f2debye(0.24,38820),f2debye(0.027,43370),f2debye(0.14,46840),f2debye(0.12,48320)];  # Debye
qiuA=1.602*10**(-19)*np.array([0.087,0.280,0.090,0.193,0.124]);  # C
liuA=[3.96,2.70,2.66,2.72,3.87];       # Ang
transcenterA=[36700,38800,43370,46840,48320];  # cm-1
numtranA=5;
#widthA=[1550,1550,700,1050,1000];
widthA=[1550,1550,700,1050,1550];

#Replaces axisdetermine.m
def make_parameter_table(coords,angle,qiu,liu,numtran,file_name):
    import csv
    import sympy as sp

    factorw = np.ones(numtran)
    width = 1550*factorw
    angle = np.deg2rad(angle)

    #coords contains the (x,y,z) coordinates of
    #the three atoms in the plane of the nucleobase
    #required to define the dipole transition moments.
    
    A = coords[0,:]
    B = coords[1,:]
    C = coords[2,:]

    #(x,y,z) of the C4 and C5 atoms for ApA
    origin = 0.5*(A + B)
    vAB = (B-A)/np.linalg.norm(B-A)
    vAC = (C-A)/np.linalg.norm(C-A)
    vBC = (C-B)/np.linalg.norm(C-B)
    normABC = np.cross(vAC,vBC)/np.linalg.norm(np.cross(vAC,vBC))

    XO = origin[0]
    YO = origin[1]
    ZO = origin[2]
    XA = A[0]
    YA = A[1]
    ZA = A[2]
    XC = C[0]
    YC = C[1]
    ZC = C[2]

    AAA = normABC[0]
    BBB = normABC[1]
    CCC = normABC[2]
    OAdis=np.linalg.norm(A-origin)

    xqp = np.zeros(numtran)
    yqp = np.zeros(numtran)
    zqp = np.zeros(numtran)
    xqn = np.zeros(numtran)
    yqn = np.zeros(numtran)
    zqn = np.zeros(numtran)

    for i in range(numtran):
        #print(i)

        XE = sp.Symbol("XE")
        YE = sp.Symbol("YE")
        ZE = sp.Symbol("ZE")
        eqn1 = (XE-XO)**2 + (YE-YO)**2 + (ZE-ZO)**2 - (liu[i]/2.0)**2 
        eqn2 = AAA*(XE-XO) + BBB*(YE-YO) + CCC*(ZE-ZO) 
        eqn3 = OAdis*liu[i]*np.cos(angle[i])/2.0 - (XA-XO)*(XE-XO) - (YA-YO)*(YE-YO) - (ZA-ZO)*(ZE-ZO) 
        sol = sp.solve((sp.Matrix([eqn1,eqn2,eqn3]),sp.Matrix([0,0,0])),[XE,YE,ZE])

        #print(sol)

        qpxx = np.zeros(2)
        qpyy = np.zeros(2)
        qpzz = np.zeros(2)
        dissq = np.zeros(2)

        #Make sure the solutions are actually solutions
        for k in [0,1]:
            XE = sol[k][0]
            YE = sol[k][1]
            ZE = sol[k][2]

            #print((XE-XO)**2 + (YE-YO)**2 + (ZE-ZO)**2 - (liu[i]/2.0)**2)
            #print(AAA*(XE-XO) + BBB*(YE-YO) + CCC*(ZE-ZO))
            #print(OAdis*liu[i]*np.cos(angle[i])/2.0 - (XA-XO)*(XE-XO) - (YA-YO)*(YE-YO) - (ZA-ZO)*(ZE-ZO))
            #print("\n")

            qpxx[k] = float(XE)
            qpyy[k] = float(YE)
            qpzz[k] = float(ZE)
            dissq[k] = (qpxx[k]-XC)**2 + (qpyy[k]-YC)**2 + (qpzz[k]-ZC)**2

        if angle[i] > 0 :
            if dissq[0] > dissq[1]:
                xqp[i] = qpxx[0]
                yqp[i] = qpyy[0]
                zqp[i] = qpzz[0]
            else:
                xqp[i] = qpxx[1]
                yqp[i] = qpyy[1]
                zqp[i] = qpzz[1]

            xqn[i] = 2*XO - xqp[i]
            yqn[i] = 2*YO - yqp[i]
            zqn[i] = 2*ZO - zqp[i]

        else:
            if dissq[1] > dissq[0]:
                xqp[i] = qpxx[0]
                yqp[i] = qpyy[0]
                zqp[i] = qpzz[0]
            else:
                xqp[i] = qpxx[1]
                yqp[i] = qpyy[1]
                zqp[i] = qpzz[1]

            xqn[i] = 2*XO - xqp[i]
            yqn[i] = 2*YO - yqp[i]
            zqn[i] = 2*ZO - zqp[i]

    return xqp,yqp,zqp,xqn,yqn,zqn,origin
            
def get_axes(coords,angleA,muA,transcenterA,numtranA):
    
    factorw = np.ones(numtranA)
    width = 1550*factorw
    angle = np.deg2rad(angleA)
    
    #coords contains the (x,y,z) coordinates of
    #the three atoms in the plane of the nucleobase
    #required to define the dipole transition moments.
    
    A = coords[0,:]
    B = coords[1,:]
    C = coords[2,:]
    
    #(x,y,z) of the C4 and C5 atoms for ApA
    origin = 0.5*(A + B)
    vAB = (B-A)/np.linalg.norm(B-A)
    vAC = (C-A)/np.linalg.norm(C-A)
    vBC = (C-B)/np.linalg.norm(C-B)
    normABC = np.cross(vAC,vBC)/np.linalg.norm(np.cross(vAC,vBC))
    yaxis = vAB/np.linalg.norm(vAB)
    zaxis = -normABC
    xaxis = np.cross(yaxis,zaxis)/np.linalg.norm(np.cross(yaxis,zaxis))
    
    if np.dot(np.cross(yaxis,zaxis),normABC) < 0:
        xaxis = -xaxis
        
    mu = debye2SI(muA)
    mux = np.sin(angle)*mu
    muy = np.cos(angle)*mu
    muxyz = np.zeros((numtranA,3))
    
    for i in range(numtranA):
        muxyz[i,:] = mux[i]*xaxis+muy[i]*yaxis

    return xaxis,yaxis,zaxis,muxyz

def rotate(vec,theta,point):
    import numpy as np
    
    #Ensure unit vector
    vec = vec/np.linalg.norm(vec)
    c = np.cos(np.deg2rad(theta))
    s = np.sin(np.deg2rad(theta))
    x = vec[0]
    y = vec[1]
    z = vec[2]
    
    #Source: https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    rot_mat = np.zeros((3,3))
    rot_mat[0,0] = c + vec[0]**2*(1-c)
    rot_mat[0,1] = x*y*(1-c)-z*s
    rot_mat[0,2] = x*z*(1-c) + y*s
    rot_mat[1,0] = y*x*(1-c) + z*s
    rot_mat[1,1] = c + y**2*(1-c)
    rot_mat[1,2] = y*z*(1-c) - x*s
    rot_mat[2,0] = z*x*(1-c) - y*s
    rot_mat[2,1] = z*y*(1-c) + x*s
    rot_mat[2,2] = c + z**2*(1-c)
    
    #Rotate the point
    
    rotated_point = np.dot(rot_mat,point)
    
    return rotated_point

def Coupling(qiua,qiuc,r13,r14,r23,r24):
    import numpy as np
    hbar=1.054E-34; # J.s 
    e0=8.854E-12; # permittivity of a vacuum Kg-1m-3s2C2
    diel=2.0;

    r13 = r13/(10**10);
    r14 = r14/(10**10);
    r23 = r23/(10**10);
    r24 = r24/(10**10);
    
    V = 0;
    V = V + qiua*qiuc/r13;
    V = V - qiua*qiuc/r14;
    V = V - qiua*qiuc/r23;
    V = V + qiua*qiuc/r24;
    
    V = V / (4.0*np.pi*e0*diel)

    return J2wn(V)

angle=[66,19,-15,-21,-64];  # Degrees
miuA=[f2debye(0.047,36710),f2debye(0.24,38820),f2debye(0.027,43370),f2debye(0.14,46840),f2debye(0.12,48320)];  # Debye
qiu=1.602*10**(-19)*np.array([0.087,0.280,0.090,0.193,0.124]);  # C
liu=[3.96,2.70,2.66,2.72,3.87];       # Ang
transcenterA=[36700,38800,43370,46840,48320];  # cm-1
numtran=5;
#widthA=[1550,1550,700,1050,1000];
widthA=[1550,1550,700,1050,1550];

#Hard-code for now; could also get by, say, reading the number of rows 
#in the .csv file and subtracting 1 for each .csv file (to remove
#the line due to the header). But, not going to do that right now.

nbasis = 10

#Some useful bits of code have been commented out because I don't want the extra overhead
#when performing the CD calculations for tens of thousands of structures (or more). They can
#be uncommented safely, if so desired.

stick_spectra_list = []
stick_center_list = []
basis_array_list = []
coupling_list = []
tau_list = []
assignment_list = []

assignment = []
#Generalize for an arbitrary number of frames
#Convert the trajectory into units of Angstroms
traj = 10.0*np.load('unformatted_traj.npy')
NFRS = traj.shape[0]//6
#NFRS = 100
print(NFRS,"\n")

x = wntofre(np.linspace(nmtown(200),nmtown(320),100))
basis_array = np.zeros((nbasis,len(x)))

CD = np.zeros_like(x)
stick_spectra = np.zeros(nbasis)

summalist = []
for k in range(0,6*NFRS,6):
    #print(k)
    if k % 1000 == 0: print('FRAME ',k)

    file_name = 'A1'
    coords = traj[k:k+3,:]
    xqp1,yqp1,zqp1,xqn1,yqn1,zqn1,origin1 = make_parameter_table(coords,angle,qiu,liu,numtran,file_name)

    file_name = 'A2'
    coords = traj[k+3:k+6,:]
    xqp2,yqp2,zqp2,xqn2,yqn2,zqn2,origin2 = make_parameter_table(coords,angle,qiu,liu,numtran,file_name)

    R12 = origin2 - origin1
    Rhat = R12/np.linalg.norm(R12);

    Emat = np.zeros(nbasis)
    mumat = np.zeros((nbasis,3))
    Hmat = np.zeros((nbasis,nbasis))

    for i in range(nbasis):
        Emat[i] = transcenterA[i%5]

    coordqp = np.zeros((nbasis,3))
    coordqn = np.zeros((nbasis,3))
    qp = np.zeros(nbasis)
    for i in range(nbasis):
        if i < nbasis//2:
            coordqp[i,:] = np.array((xqp1[i],yqp1[i],zqp1[i]))
            coordqn[i,:] = np.array((xqn1[i],yqn1[i],zqn1[i]))
            qp[i] = qiu[i]
        else:
            coordqp[i,:] = np.array((xqp2[i-nbasis//2],yqp2[i-nbasis//2],zqp2[i-nbasis//2]))
            coordqn[i,:] = np.array((xqn2[i-nbasis//2],yqn2[i-nbasis//2],zqn2[i-nbasis//2]))
            qp[i] = qiu[i-nbasis//2]  

    for i in range(nbasis):
        mumat[i,:] = -(qp[i]*coordqp[i,:] - qp[i]*coordqn[i,:])/(10**10)

    for i in range(nbasis):
        for j in range(i,nbasis):
            if i != j:
                if i <= nbasis//2 and j - i >= nbasis//2:
                    r13 = np.linalg.norm(coordqp[j,:] - coordqp[i,:])
                    r14 = np.linalg.norm(coordqn[j,:] - coordqp[i,:])
                    r23 = np.linalg.norm(coordqp[j,:] - coordqn[i,:])
                    r24 = np.linalg.norm(coordqn[j,:] - coordqn[i,:])
                    Hmat[i,j] = Coupling(qp[i],qp[j],r13,r14,r23,r24)
                    Hmat[j,i] = Hmat[i,j]
            else:
                Hmat[i,i] = Emat[i]

    for i in range(nbasis):
        for j in range(i,nbasis):
            if i == j: Hmat[i,j] = Emat[i]
            else:
                if i < nbasis//2:
                    if j >= nbasis//2:
                        r13 = np.linalg.norm(coordqp[j,:] - coordqp[i,:])
                        r14 = np.linalg.norm(coordqn[j,:] - coordqp[i,:])
                        r23 = np.linalg.norm(coordqp[j,:] - coordqn[i,:])
                        r24 = np.linalg.norm(coordqn[j,:] - coordqn[i,:])
                        Hmat[i,j] = Coupling(qp[i],qp[j],r13,r14,r23,r24)
                        Hmat[j,i] = Hmat[i,j]
                else:
                    if j < nbasis//2:
                        r13 = np.linalg.norm(coordqp[j,:] - coordqp[i,:])
                        r14 = np.linalg.norm(coordqn[j,:] - coordqp[i,:])
                        r23 = np.linalg.norm(coordqp[j,:] - coordqn[i,:])
                        r24 = np.linalg.norm(coordqn[j,:] - coordqn[i,:])
                        Hmat[i,j] = Coupling(qp[i],qp[j],r13,r14,r23,r24)
                        Hmat[j,i] = Hmat[i,j]


    eigval,eigvec = np.linalg.eigh(Hmat)
    #np.savetxt('Hmat.dat',Hmat,fmt="%7.6e")
    #np.savetxt('eigenvectors_python.dat',eigvec,fmt="%7.6e")
    #np.savetxt('eigvals_python.dat',eigval)
    #np.save(str(struct)+'/eigval_python.npy',eigval)
    #np.save(str(struct)+'/eigvec_python.npy',eigvec)
    #np.save(str(struct)+'/coupling_python.npy',Hmat)

    #Look at a single frame for now
    #NFRS = 1

    #Define arrays
    CDsig = np.zeros(nbasis)
    widthd = np.zeros(nbasis)
    sumyCD = np.zeros_like(x)

    Coupmat = Hmat - np.diag(np.diagonal(Hmat)) 
    for n in range(1):
        for m in range(n+1,2):
            #Calculate the rotational strength
            Rk = np.zeros(nbasis)
            Ek = np.zeros(nbasis)
            #degenlist = []
            #nondegenlist = []
            for i in range(nbasis):
                eigveckk = eigvec[:,i]
                Rk[i],degen,nondegen = calcRk_sign(eigveckk,R12,wn2J(Emat),wn2J(Coupmat),mumat)
                Ek[i] = eigval[i] 
                #degenlist.append(degen)
                #nondegenlist.append(nondegen)

            yd = np.zeros_like(x)
            for i in range(nbasis):
                widthd[i] = 0.2*8065.5

            for i in range(nbasis):
                CDsig[i] = Rk[i]*wntofre(Ek[i])/(7.659E-54)
                y = CDsig[i]*np.exp(-0.5*((x-wntofre(Ek[i]))/wntofre(widthd[i]))**2)/(np.sqrt(2*np.pi)*wntofre(widthd[i]))
                #stick_spectra[i] += CDsig[i]/(np.sqrt(2*np.pi)*wntofre(widthd[i]))
                #basis_array[i,:] += y
                yd += y

            sumyCD += yd
    summalist.append(sumyCD)

    CD = CD + sumyCD

CD = CD / NFRS

np.savetxt('RESULTS_python.dat',np.column_stack([fretonm(x),CD]))
np.savetxt('RESULTS_python_all.dat',np.column_stack([fretonm(x),CD,np.array(summalist).std(0)]))


#Check

#print(len(summalist))
#chk = np.zeros_like(CD)
#for i in range(len(summalist)):
#    chk += summalist[i]
#chk = chk / len(summalist)

#np.savetxt('RESULTS_summa.dat',np.column_stack([fretonm(x),chk,np.array(summalist).std(0)]))
#np.save('CDlist.npy',np.array(summalist))
'''colors = plt.cm.winter(np.linspace(0,1,len(summalist)))
colors1 = ['k','r','b','c','magenta']
colors2 = plt.cm.Spectral(np.linspace(0,1,basis_array.shape[0]//2))
for num,i in enumerate(summalist):
    plt.plot(fretonm(x),i,lw=3,c=colors[num])
    
plt.plot(fretonm(x),CD,lw=3,c='k')

stdlist = np.array(summalist).std(0)
plt.fill_between(fretonm(x),CD,CD+stdlist,alpha=0.5,c='grey')
plt.fill_between(fretonm(x),CD - stdlist,CD,alpha=0.5,c='grey')
plt.xticks(size=12)
plt.yticks(size=12)
plt.xlabel('wavelength, nm',fontsize=16)
plt.ylabel('CD',fontsize=16)
#plt.legend()
#plt.savefig('spectral_decomp_python_coupling_all.pdf',dpi=300)
#plt.savefig('spectral_decomp_python_coupling_all.png',dpi=300)
plt.show()
plt.close()'''
