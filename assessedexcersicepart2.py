# -*- coding: utf-8 -*-
"""AssessedExcersicePart2.ipynb
Name: Kristinn Petursson
project: Micromagnetism assessed excercise part 2
date: march 20, 2023

Code made by using Excercise 3 as a blueprint

Original file is located at
    https://colab.research.google.com/drive/1eM7ng7uAHfmdaOiXUZpCY6ON_z4gIZBB

# **Part 2**

## Setup

We start by confirming our runtime has access to a GPU (if not, correct this using Runtime -> Change Runtime Type) and installing mumax3.
"""

! echo "This machine runs" $(uname)
! echo
! nvidia-smi

# Download the mumax3 binary
!wget https://mumax.ugent.be/mumax3-binaries/mumax3.10_linux_cuda10.1.tar.gz
!tar -xvf mumax3.10_linux_cuda10.1.tar.gz
!rm mumax3.10_linux_cuda10.1.tar.gz
!rm -rf mumax3.10 && mv mumax3.10_linux_cuda10.1 mumax3.10

#update the PATH environment variable
import os
os.environ['PATH'] += ":/content/mumax3.10"

"""## Simulation

We begin by specifying the shape and materiaol parameters of our system. To simulate an "island shape" system, we will use "periodic boundary conditions" by surrounding our system with copies of itself in a $4\times 4$ grid. 

We then excite our system to find its resonance frequency. Our strategy will be to first allow our film to relax in a static $<11>$ magnetic field $B$ before applying a field pulse in $z$ with a "sinc" profile
\begin{equation}
B_z (t) = A \frac{\sin{\left(2\pi f_c (t-t_0)\right)}}{2\pi f_c (t-t_0)}
\end{equation}
where $A$ is the peak amplitude of the pulse, $f_c$ is our cut-off frequency, $t_0$ is the time at which our pulse reaches its peak magnitude.

In this simulation we will set $f_c$ to 25 GHz, so that all modes up to this frequency are excited by the pulse. We will need to make sure that we are sampling the magnetisation at a rate that allows us to capture processes occurring up to this frequency.
"""

!mumax3 PyFilmPulse.mx3

!mumax3-convert -png -arrows 10 PyFilmPulse.out/*.ovf

"""## Postprocessing

We will now Fourier transform the magnetisation data we saved to our table to see at which frequencies our thin film responded to the pulse.
"""

import numpy as np
from matplotlib import pyplot as plt

# Import data from table.txt
folder = 'PyFilmPulse.out'

data = np.loadtxt(folder+'/table.txt',skiprows=1)
times = np.array([d[0] for d in data])
mx = np.array([d[1] for d in data])
my = np.array([d[2] for d in data])
mz = np.array([d[3] for d in data])
Bz = np.array([d[7] for d in data])

# Fourier transform the z-component of the magnetisation

spectrum = np.fft.fftshift(np.fft.fft(mz))
intensity = np.abs(spectrum)**2
timeStep = times[1]-times[0]
freqs = np.fft.fftshift(np.fft.fftfreq(len(times),timeStep)) # Converts our time steps to frequency bins for the FFT.
resFreq = np.abs(freqs[np.argmax(intensity)]) # Frequency for which the Fourier intensity is at its maximum value.

# Plot the z-component of the magnetic field as a function of time
# Note that our table starts at t=50 ns since we relaxed our system before we started saving data.
plt.plot(times/1e-9 - 10,Bz/1e-3,label='$B_z$')
plt.xlim(0,5)
#plt.xlim(0,times[-1]/1e-9-10) # Uncomment this line to plot every time step that was saved.
plt.xlabel('$t$ (ns)', size=18)
plt.ylabel('$B_z$ (mT)', size=18)
plt.show()

# Plot the z-component of the magnetisation as a function of time
plt.plot(times/1e-9 - 10,mz,label='$m_z$')
plt.xlim(0,10)
#plt.xlim(0,times[-1]/1e-9-10)
plt.xlabel('$t$ (ns)', size=18)
plt.ylabel('$m_z$', size=18)
plt.show()

# Plot the Fourier intensity as a functino of frequency. This is the spectrum for our thin film.
plt.plot(freqs/1e9,intensity/intensity.max())
#plt.xlim(0,10)
plt.ylim(0,1)
plt.xlabel('Frequency (GHz)', size=18)
plt.ylabel('Intensity (a.u.)', size=18)
plt.show()
#creating a file used for plotting a scatter plot once all simulations are done
filename = 'PyFilmPulse.mx3'
with open(filename, 'r') as file:
    lines = file.readlines()
initLine = str(lines[68])
magneticField = initLine[6:-1]
f = open("resonanceFreq.txt","a")
f.write('%s %s\r' % ((magneticField), str(resFreq/1e9)))
print('Resonance frequency = '+str(resFreq/1e9)+' GHz')

"""## Kittel Formula

The condition for ferromagnetic resonance in a ferromagnetic thin film with in-plane magnetisation is given by the Kittel formula which relates the resonance frequency $f$ to the magnetic field $B$ and the magnetisation of the film $M$

\begin{equation}
f(B) = \frac{|\gamma|}{2\pi}\sqrt{B(B+\mu_0 M)}
\end{equation}

Here we change the value if the magnetic field B from -100mT, -50mT, 0mT, 50mT, and 100mT
"""

# This cell edits PyFilmPulse.mx3 to alter the
# value of the static external field Bx.
##############################################
# Set the new value in Tesla here.
newField = 100e-3
##############################################

filename = 'PyFilmPulse.mx3'
with open(filename, 'r') as file:
    lines = file.readlines()

initLine = str(lines[68])

lines[68] = 'Bx := '+str(newField)+'\n'

with open(filename, "w") as file:
    file.writelines(lines)

print('Changed '+initLine+' to '+'Bx := '+str(newField))
#making scatter plot:

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import Rbf

#This program shows how to write data in a text file.

f = open("resonanceFreq.txt","r")

data = np.loadtxt("resonanceFreq.txt")
b = np.array( [d[0] for d in data] )
freq = np.array([d[1] for d in data] )

f = open("resonanceFreq110.txt","r")

data = np.loadtxt("resonanceFreq110.txt")
b2 = np.array( [d[0] for d in data] )
freq2 = np.array([d[1] for d in data] )

b_new = np.linspace(b.min(), b.max(),500)

rbf = Rbf(b2, freq2, function = 'thin_plate', smooth = 0.001)
y2_smooth = rbf(b_new)
plt.plot(b2, freq2, 'ro')
plt.plot(b_new, y2_smooth, 'r', label='m = uniform(1,1,0)')

rbf = Rbf(b, freq, function = 'thin_plate', smooth = 0.001)
y_smooth = rbf(b_new)
plt.plot(b, freq, 'bo')
plt.plot(b_new, y_smooth, 'b', label='m = uniform(-1,-1,0)')
plt.legend()
plt.xlabel('Magnetic field, B (mT)', size=18)
plt.ylabel('Frequency, f (GHz)', size=18)

plt.show()