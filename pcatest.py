#!/bin/python

import pylab as p
import numpy as n
from numpy.linalg import svd
from scipy.signal import sawtooth

xs=n.linspace(0,6.*n.pi,200)

def line(xs):
    '''This function just defines a straight line'''
    ys=[]
    for x in xs:
        ys.append(x/8+0.3)
    return ys

# Specify number of 'spectra' to generate
nspec=15

speclist=[]

# This is all plotting, ignore it.
figure1=p.figure(figsize=(20,6))
ax=p.subplot2grid(shape=(3,3),loc=(0,0))
p.plot(line(xs))
ax=p.subplot2grid(shape=(3,3),loc=(1,0))
p.plot(n.sin(xs))
ax=p.subplot2grid(shape=(3,3),loc=(2,0))
p.plot(sawtooth([x*2 for x in xs]))
ax=p.subplot2grid(shape=(3,3),loc=(0,1),rowspan=3)

for i in range(0,nspec):
    # For each spectrum, generate three random numbers. 
    # Distributions are flat, centred on zero, but have differing amplitudes.
    randval1=0.2*(n.random.rand()-0.5)
    randval2=0.6*(n.random.rand()-0.5)
    randval3=0.1*(n.random.rand()-0.5)
    
    # Multiply three different input components (sin wave, sawtooth and stright line) by random values
    sinespec=[randval1 * s for s in n.sin(xs)]
    sinespec2=[randval3 * s for s in sawtooth([x*2 for x in xs])]
    linespec=[randval2 * l for l in line(xs)]

    # Generate some noise
    noisespec=[(n.random.rand()-0.5)*0.1 for x in xs]
    
    # Add all components together, with noise, and append to spectrum list
    speclist.append([a+b+c+d for a,b,c,d in zip(sinespec,linespec,noisespec,sinespec2)])
    p.plot(speclist[-1])

specarray=n.array(speclist)
specarray=n.transpose(specarray)

specfile='specfile.dat'
n.savetxt(specfile,specarray)

# All the analysis is in this line. This decomposes the spectrum array into constituent components.
a,b,c=svd(specarray)

# Matrix a now contains all the output components, vector b their strengths. Ignore c. 
a=n.transpose(a)

# Plot the components
ax=p.subplot2grid(shape=(3,3),loc=(0,2))
p.plot(a[0])
ax=p.subplot2grid(shape=(3,3),loc=(1,2))
p.plot(a[1])
ax=p.subplot2grid(shape=(3,3),loc=(2,2))
p.plot(a[2])
p.show()

#print 'OUTPUT FUNCTIONS'
#for a0,a1,a2 in zip(a[0],a[1],a[2]):
	#print a0,a1,a2


total=0
for val in b:
    total+=val**2

vals=[val**2/total for val in b]
#print '\nEIGENVALUES'
#for val in vals:
	#print val

# p.plot(range(1,len(vals)+1),n.log10(vals),ls='None',marker='x')
# p.show()