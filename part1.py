#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 15:32:03 2023

@author: kristinnpetursson
Project: Mini-Project "Micromagnetics"
Module: MATE70002 - Theory and Simulation
Date: 19 March, 2023
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from shapely.geometry import LineString
from math import sqrt

# Function: intersectionXCoordinance - finding the x coordinance by interpolating
# arrayX: an array of values that are plotten on x axis
# arrayY: an array of values that are plotted on y axis
# valueYClosest0: a float, closest value of arrayY to 0
# return: the x-coordinance of the intersection between x- axis (y=0) and hysterisis
def intersectionXCoordinance(arrayX, arrayY , indexYClosest0):
    valueY1 = float(arrayY[indexYClosest0])
    valueX1 = float(arrayX[indexYClosest0])
    indexNext = indexYClosest0 +1 
    indexPrevious = indexYClosest0 - 1
    valueY2 = 0
    valueX2 = 0
    #4 possibilities, ++, --, +-, -+
    if((arrayY[indexNext] > 0) & (valueY1 > 0)):
        #find value of the previous
        #print("both higher than 0")
        valueY2 = float(arrayY[indexNext])
        valueX2 = float(arrayX[indexNext])
        return interpolating(valueX1, valueX2, valueY1, valueY2)
    
    elif((arrayY[indexNext] < 0) & (valueY1 < 0)):
        #find value of the next
        #print("both lower than 0")
        valueY2 = float(arrayY[indexPrevious])
        valueX2 = float(arrayX[indexPrevious])
        return interpolating(valueX1, valueX2, valueY1, valueY2)
    
    elif( (arrayY[indexNext] <  valueY1) or (arrayY[indexNext] >  valueY1) ): 
        #downward trending or upwards
        valueY2 = float(arrayY[indexNext])
        valueX2 = float(arrayX[indexNext])
        return interpolating(valueX2, valueX1, valueY2, valueY1)
        
    else: #upward trending
       return print("problem with data")    

#function: interpolating
# return: the x value of interseption
def interpolating(x1, x2, y1, y2):
    result =  np.around( float(x1) - y1*(x2-x1)/(y2-y1), 4)
    print(result)
    return result

folder = "magnetiteParticle.out"
data = np.loadtxt('/Users/kristinnpetursson/Downloads/Miniproj/Part 1 - full gap angle 45/table.txt',skiprows=1)

mx = np.array( [d[1] for d in data] )
my = np.array([d[2] for d in data] )
mtotal = np.around(( (sqrt(2)/2)*mx + (sqrt(2)/2)*my) , 4)

Bextx = np.array([d[4] for d in data])
Bexty = np.array([d[5] for d in data] )
Btotal = np.around(((sqrt(2)/2)*Bextx + (sqrt(2)/2)*Bexty), 4)

#Delete beginning of graph when system is "ramping" up
index = np.argmax(mtotal > 0.55)
m_newTotal = mtotal[index:-1]
b_newTotal = 100*Btotal[index:-1]

#find intersections
closestToZero1 = np.abs(m_newTotal - 0).argmin()
x1 = intersectionXCoordinance(b_newTotal, m_newTotal, closestToZero1)
#Fnding lowest value second intersection with 0
sizeArray = int(2*len(m_newTotal)/3) #the last third part  of array.
newArrayB = b_newTotal[sizeArray:-1]
newArrayM = m_newTotal[sizeArray:-1]
closestZero2 = np.abs(newArrayM - 0).argmin() #closest value to 0
x2 = intersectionXCoordinance(newArrayB, newArrayM, closestZero2)

#Intersection of m_sat and hysterisis
closestToSat = np.abs(m_newTotal - 0.8).argmin()
x3 = m_newTotal[closestToSat]

#plotting
plt.axvline(x=x2, color='gray', linestyle='--')
plt.axhline(y=0, color='gray', linestyle='--')
plt.axvline(x=x3, color='red', linestyle='--')
plt.axhline(y=0.8, color='red', linestyle='--')
plt.plot(b_newTotal ,m_newTotal)

#display text on plot
stringX1 = '\u03BC' + '\u2080' + '$H_C$ = ' + str(x1) + ' mT'
stringX2 = '\u03BC' + '\u2080' + '$H_C$ = ' + str(x2) + ' mT'
stringX3 = '$m_r / m_{sat}$'
plt.text(-10.5, 0.10, stringX1, color='black')
plt.text(4.75, 0.10, stringX2 , color='black')
plt.text(-3.5, 0.7, stringX3, color='red')
#labels and showing plot
plt.xlabel('$\u03bc_0H$ (mT)', size=18)
plt.ylabel('$m / m_{sat} (mT)$', size=18)
plt.ylim(-1,1)
plt.show()

