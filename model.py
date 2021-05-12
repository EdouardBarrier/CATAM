# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 15:50:49 2021

@author: Edouard Barrier
"""
import numpy as np
from matplotlib import pyplot as plt

pi=np.pi
num=5000
realisation=np.zeros((num,10))
progress=0
reached=False
Rt=100
vmax=1
tried=0
rng=np.random.default_rng() #having it here instead of in the while loop 
#saves a huge amount of time as creating the random generator is very expensive

def f(et):
    a=(1.5/(pi**3))
    b=((3-4*et)*((2*et)**0.5))/(1-2*et)
    c=3*np.arcsinh((2*et / (1-2*et))**0.5)
    return (a*(b-c))

def g(et):
    a=(8*pi**2)
    b=((1-2*et)**0.5)*(3-14*et-8*(et**2))/(12*(et**2))                
    c=((1-6*et+16*(et**2))/((2*et)**2.5))*np.arccos(-((1-2*et)**0.5))
    return (a*(b-pi+c))
                          
while (reached==False):
    taken=True
    tried+=1
    ran1=rng.random()
    ran2=rng.random()
    ran3=rng.random()
    ran4=rng.random()
    ran5=rng.random()
    ran6=rng.random()
    xhi=rng.random()
    R=Rt*rng.random()
    V=vmax*rng.random()
    azi=2*pi*rng.random()
    polar=pi*rng.random()
    vazi=2*pi*rng.random()
    vpolar=pi*rng.random()
    eta= 0.5* (1 -(R/(R+1))**2 -V**2)
    if (eta<0):
        taken=False
        continue
    func=f(eta)
    quan=(R**2)*(V**2)*func
    if (xhi>(quan/0.002884)): #monte carlo done here
        taken=False
    if (taken==True):
        realisation[progress,0]=progress #accept all the data
        realisation[progress,1]=R
        realisation[progress,2]=V
        realisation[progress,3]=azi
        realisation[progress,4]=polar
        realisation[progress,5]=eta
        realisation[progress,6]=func
        realisation[progress,7]=quan
        realisation[progress,8]=vazi
        realisation[progress,9]=vpolar
        progress+=1
    if ((progress<num)==False):
        reached=True
        
"""
nbine=200  #DENSITY OF STATES QUESTION 3
ebins=np.zeros((nbine,3))
width=0.5/nbine #eta range
for i in range((num)):#sort into eta bins
    location= int(realisation[i,5]//width)
    ebins[location,1]+=1
for i in range((nbine)):#get /de term in, plot analytic expression
    ebins[i,1]=ebins[i,1]/(width*num)
    ebins[i,0]=(i+0.5)*width
    et1=ebins[i,0]
    ebins[i,2]=f(et1)*g(et1)
xval=ebins[:,0]
numer=ebins[:,1]
analy=ebins[:,2]
plt.title('Figure 2: Energy Distribution, N=50000')
plt.xlabel('\u03B5')
plt.ylabel('dM(\u03B5)/d\u03B5')
plt.plot(xval,numer,'b',label='Numerical values')
plt.plot(xval,analy,'r',label='Analytic Expression')
plt.legend(loc="upper right")
plt.show()
"""

"""
space=np.zeros((num,3)) #SPATIAL DISTRIBUTION QUESTION 4
for i in range((num)):
    rad=realisation[i,1]
    azimuthal=realisation[i,3] 
    polar=realisation[i,4]
    space[i,0]=rad*np.sin(polar)*np.cos(azimuthal) #x direction
    space[i,1]=rad*np.sin(polar)*np.sin(azimuthal) #y direction
    space[i,2]=rad*np.cos(polar) #z direction
fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
ax.scatter(space[:,0],space[:,1],space[:,2],c='b')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title('Figure 3: 3D realisation of the stars')
plt.show()
"""

"""
vspace=np.zeros((num,3)) #VELOCITY DISTRIBUTION QUESTION 4
for i in range((num)):
    vel=realisation[i,2]
    azimuthal=realisation[i,8] 
    polar=realisation[i,9]
    vspace[i,0]=vel*np.sin(polar)*np.cos(azimuthal) #x direction
    vspace[i,1]=vel*np.sin(polar)*np.sin(azimuthal) #y direction
    vspace[i,2]=vel*np.cos(polar) #z direction
fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
ax.scatter(vspace[:,0],vspace[:,1],vspace[:,2],c='b')
ax.set_xlabel('Vx')
ax.set_ylabel('Vy')
ax.set_zlabel('Vz')
plt.title('Figure 4: Stars in velocity space')
plt.show()
"""

"""
nbinr=500  #RADIUS BINS QUESTION 5
rbins=np.zeros((nbinr,3))
width=Rt/nbinr #r range
for i in range((num)):
    location=int(realisation[i,1]//width)
    rbins[location,1]+=1
for i in range((nbinr)):
    rbins[i,0]=(i+0.5)*width
    r=rbins[i,0]
    vol=(4*pi/3)*((r+0.5*width)**3 - (r-0.5*width)**3)
    rbins[i,1]=rbins[i,1]/(vol*num)
    rbins[i,2]=(3/(4*pi))*((1+r)**(-4))
xval=rbins[:,0]
rval=rbins[:,1]
logxval=np.log10(rbins[:,0])
logrval=np.log10(rbins[:,1])
logeval=np.log10(rbins[:,2])
plt.plot(logxval,logrval,'b',label='Numerical')
plt.plot(logxval,logeval,'r',label='Analytic')
plt.legend(loc="upper right")
plt.title('Figure 6: log-log plot of density')
plt.xlabel('log(r)')
plt.ylabel('log(density)')
plt.show()
"""

"""
track=np.zeros((num,2)) #PLOT TWO VARIABLES AGAINST EACH OTHER
track[:,0]=realisation[:,5]
track[:,1]=realisation[:,7]
track=track[np.argsort(track[:,0])] #sorts it according to that column
plt.plot(track[:,0],track[:,1])
plt.title('Figure 17: \u03B5 vs Probability')
plt.xlabel('\u03B5')
plt.ylabel('Relative Probability')
plt.show()
"""

"""
nbinv=200  #VELOCITY BINS
vbins=np.zeros((nbinv,3)) #1 for x, 2 for numerical, 3 for analytic
width=vmax/nbinv #v range
for i in range((num)):
    location=int(realisation[i,2]//width)
    vbins[location,1]+=1
for i in range((nbinv)):
    a=0.18
    v=(i+0.5)*width
    vbins[i,0]=v
    vbins[i,1]=vbins[i,1]/(num*width)
    vbins[i,2]=np.sqrt(2/pi)*(((v**2)*np.exp(-v**2/(2*a**2)))/(a**3))
plt.plot(vbins[:,0],vbins[:,1],'b',label='Numerical Velocity Distribution')
plt.plot(vbins[:,0],vbins[:,2],'r',label='Maxwell Distribution')
plt.legend(loc="upper right")
plt.title('Figure 14: Velocity distribution')
plt.xlabel('v')
plt.ylabel('Probability distribution')
plt.show()
"""

"""
nbind=200 #DISPERSION QUESTION 6
dbins=np.zeros((nbind,4))
width=Rt/nbind #we need to do dispersion per radius bin
for i in range((num)):
    location=int(realisation[i,1]//width)
    dbins[location,1]=dbins[location,1]+(realisation[i,2])**2
    dbins[location,2]+=1
for i in range((nbind)):
    if dbins[i,2]>0:#avoids divide by 0
        dbins[i,1]=dbins[i,1]/dbins[i,2]
    r=(i+0.5)*width
    dbins[i,0]=r
    dbins[i,3]=(1+6*r)/(10 * (1+r)**2)
plt.plot(dbins[:,0],dbins[:,1],'b',label='Numerical dispersion')
plt.plot(dbins[:,0],dbins[:,3],'r',label='Analytic dispersion')
plt.legend(loc="upper right")
plt.title('Figure 7: Velocity dispersion for different radius bins')
plt.xlabel('r')
plt.ylabel('$\sigma$^2(r)')
plt.show()
"""

"""
nbinl=200 #ANGULAR MOMENTUM QUESTION 6
amval=np.zeros((nbinl,4))#r in 0, Lx/Ly/Lz in 1/2/3
width=Rt/nbinl
for i in range((num)):
    location=int(realisation[i,1]//width)
    x=(realisation[i,1])*np.sin((realisation[i,4]))*np.cos((realisation[i,3])) #x direction
    y=(realisation[i,1])*np.sin((realisation[i,4]))*np.sin((realisation[i,3])) #y direction
    z=(realisation[i,1])*np.cos((realisation[i,4])) #z direction
    vx=(realisation[i,2])*np.sin((realisation[i,9]))*np.cos((realisation[i,8]))#r sintheta cosphi
    vy=(realisation[i,2])*np.sin((realisation[i,9]))*np.sin((realisation[i,8]))
    vz=(realisation[i,2])*np.cos((realisation[i,9]))
    amval[location,1]=amval[location,1]+((y*vz - z*vy)/r)/num #num is to get the m in there, normalising to M=1
    amval[location,2]=amval[location,2]+((z*vx - x*vz)/r)/num
    amval[location,3]=amval[location,3]+((x*vy - y*vx)/r)/num
for i in range((nbinl)):
    amval[i,0]=(i+0.5)*width
plt.plot(amval[:,0],amval[:,1],'b',label='Angular momentum in x direction')
plt.plot(amval[:,0],amval[:,2],'r',label='Angular momentum in y direction')
plt.plot(amval[:,0],amval[:,3],'g',label='Angular momentum in z direction')

plt.legend(loc="upper right")
plt.title('Figure 8: Angular Momentum at different radii')
plt.xlabel('r')
plt.ylabel('L')
plt.show()
"""

"""
nbina=200 #ANISOTROPY QUESTION 7
anival=np.zeros((nbina,5)) #r in 0, ani in 1, v^2 in 2, vr^2 in 3, n in 4
width=Rt/nbina
for i in range((num)): #add velocity components to correct bins as required
    location=int(realisation[i,1]//width)
    r=realisation[i,1]
    v=realisation[i,2]
    xdir=np.sin((realisation[i,4]))*np.cos((realisation[i,3])) #x direction
    ydir=np.sin((realisation[i,4]))*np.sin((realisation[i,3])) #y direction
    zdir=np.cos((realisation[i,4])) #z direction
    vx=v*np.sin((realisation[i,9]))*np.cos((realisation[i,8]))#r sintheta cosphi
    vy=v*np.sin((realisation[i,9]))*np.sin((realisation[i,8]))
    vz=v*np.cos((realisation[i,9]))
    anival[location,2]=anival[location,2] + v**2
    anival[location,3]=anival[location,3] + (xdir*vx + ydir*vy + zdir*vz)**2
    anival[location,4]+=1
for i in range((nbina)):#normalising, finding anisotropy
    anival[i,0]=(i+0.5)*width
    if anival[i,4]>0: #avoids divide by 0
        anival[i,2]=anival[i,2]/anival[i,4] #mean v^2
        anival[i,3]=anival[i,3]/anival[i,4] #mean vr^2
        anival[i,1]= 1 - (anival[i,2]-anival[i,3])/(2*anival[i,3]) #anisotropy
plt.plot(anival[:,0],anival[:,1])
plt.title('Figure 9: Anisotropy at different radii')
plt.xlabel('r')
plt.ylabel('Anisotropy')
plt.ylim([-1,1]) 
plt.show()
"""

""" #POTENTIAL QUESTION 8
#have the m(r)/r bit, now need to add the mass outside
phi=np.zeros((num,3)) #r 0, numerical potential 1, analytic potential 2,
phi[:,0]=np.sort(realisation[:,1]) #array of length num
phioutside=0
for i in range((num)):
    j=num-i-1 #goes from 4999 to 0, work through array backwards
    r=phi[j,0]
    phi[j,1]=-j/(num*r) + phioutside #-GM(r)/R, M(r)=j/num, for the outernmost star there is no mass outside
    phi[j,2]=(-1/2)*(1 - (r/(r+1))**2)
    phioutside=phioutside - (1/num)*(1/r) #as the outer mass component is a sum of m/r
plt.plot(phi[:,0],phi[:,1],'b', label='Numerical potential')
plt.plot(phi[:,0],phi[:,2],'r',label='Analytic potential')
plt.legend(loc="lower right")
plt.title('Figure 10: Numerical vs Analytic Potential')
plt.xlabel('r')
plt.ylabel('$\Phi$(r)')
plt.show()
"""

"""
nbind=200 # ESCAPE VELOCITY QUESTION 8
esbins=np.zeros((nbind,4))
width=Rt/nbind #we need to do dispersion per radius bin
for i in range((num)):
    location=int(realisation[i,1]//width)
    esbins[location,1]=esbins[location,1]+(realisation[i,2])**2
    esbins[location,2]+=1
for i in range((nbind)):
    if esbins[i,2]>0:#avoids divide by 0
        esbins[i,1]=esbins[i,1]/esbins[i,2]
    r=(i+0.5)*width
    esbins[i,0]=r
    esbins[i,3]=1-(r/(r+1))**2 #use analytic potential
plt.plot(esbins[:,0],esbins[:,1],'b',label='Dispersion')
plt.plot(esbins[:,0],esbins[:,3],'r',label='Escape velocity')
plt.legend(loc="upper right")
plt.title('Figure 13: Escape velocity for high radius, larger realisation')
plt.xlabel('r')
plt.ylabel('Velocity')
plt.xlim([50,100])
plt.ylim([-0.01,0.1])
plt.show()
"""