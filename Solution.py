# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 21:10:16 2021

@author: phild
"""
import numpy as np
from numpy import array, linspace, vstack, transpose, concatenate, zeros, sign, diff
from numpy.random import uniform

from scipy.optimize import brentq

import matplotlib.pyplot as plt
plt.close('all')


def plotsegs(BP,C,fmt = 'r'):
    plt.axis('equal')
    plt.plot(*C.transpose(),'ko')
    segments = BP
    for i in transpose(segments,(1,0,2)):
        plt.plot(*EvaluateCubicBezier(i),fmt)

eps = 10**-8

B = array([[1,0,0,0],[-3,3,0,0],[3,-6,3,0],[-1,3,-3,1]])
cubic_powers = array(range(0,4))
def EvaluateCubicBezier(P,t = None):
    if t is None:
        t = linspace(0,1,100)
    t = np.atleast_2d(t).transpose()
    return (t**cubic_powers @ B @ array(P)).transpose()

def DeCasteljau(P,weight = 0.5):
    #performs de casteljau's algorithm along the first axis of the input array
    M = [array(P)]
    while len(M[-1]) != 1:
        M.append(M[-1][1:]*weight + M[-1][:-1]*(1-weight))
    P1 = array([i[0] for i in M])
    P2 = array([i[-1] for i in reversed(M)])
    return (P1,P2)

def CountDifferenceSignChanges(P):
    D = sign(diff(distances[:,0]))
    return np.sum(D[:-1] != D[1:],axis = 0)

def DifferentiateBezier(P):
    return len(P)*(P[1:] - P[:-1])

def TransformCubicDistanceSquared(P0,P1,P2,P3,C):
    #For a cubic bezier and any point, find the control points of the sextic 
    #bezier that describes the squared difference between the point and the cubic
    T = array( [C**2 - 2*C*P0 + P0**2,
                C**2 + C*(-P0 - P1) + P0*P1,
                C**2 + C*(-2/5*P0 - 6/5*P1 - 2/5*P2) + (2/5)*P0*P2 + (3/5)*P1**2,
                C**2 + C*(-1/10*P0 - 9/10*P1 - 9/10*P2 - 1/10*P3) + (1/10)*P0*P3 + (9/10)*P1*P2,
                C**2 + C*(-2/5*P1 - 6/5*P2 - 2/5*P3) + (2/5)*P1*P3 + (3/5)*P2**2,
                C**2 + C*(-P2 - P3) + P2*P3,
                C**2 - 2*C*P3 + P3**2] )
    return np.sum(T,axis = 2).transpose()

def ExtractBeziers(P):
    #Get cubic beziers from a uniform cubic B-spline
    P1 = (P[1:] + 2*P[:-1])/3
    P2 = (2*P[1:] + P[:-1])/3
    Ptemp = (P1[1:] + P2[:-1])/2
    P0 = vstack((P[0],Ptemp))
    P3 = vstack((Ptemp,P[-1]))
    return array([P0,P1,P2,P3])

def FilterCandidates(segments,distances,C):
    #filters candidates first by circular clipping, then by whether the difference polygon contains sign changes
    endpts = concatenate((segments[0],np.atleast_2d(segments[3,-1])))
    threshold = min(np.sum((endpts - C)**2,axis = 1))
    
    #index to only control polygons inside the test point radius
    boolean_indexer = np.any(distances <= threshold + eps,axis = 0)
    segments_by_threshold = segments[:,boolean_indexer,:]
    distances_by_threshold = distances[:,boolean_indexer]
    
    #index to only control polygons with one or more sign changes in the differences between sequential points (indicating local optima)
    direction_changes = sign(diff(distances_by_threshold[:-1,:],axis = 0)) != sign(diff(distances_by_threshold[1:,:],axis = 0))
    boolean_indexer = np.sum(direction_changes,axis = 0) >= 1
    segments_by_minimum = segments_by_threshold[:,boolean_indexer,:]
    distances_by_minimum = distances_by_threshold[:,boolean_indexer]

    return (segments_by_minimum, distances_by_minimum)

C = uniform(-100,100,(1,2))
P = uniform(-100,100,(100,2))


BP = ExtractBeziers(P)
segments = ExtractBeziers(P)
distances = TransformCubicDistanceSquared(*BP,C).transpose()

maxiter = 50
for i in range(0,maxiter):
    segments, distances = FilterCandidates(segments, distances, C)
    if segments.shape[1] == 1 and CountDifferenceSignChanges(distances) == 1:
        DiffPts = DifferentiateBezier(distances[:,0])
        def distance_function_derivative(t):
            return DeCasteljau(DiffPts,t)[0][-1]
        param = brentq(distance_function_derivative,0,1)
        closest_pt = EvaluateCubicBezier(segments[:,0],param)
        break
    if segments.shape[1] == 0:
        #first or last point is the optimum
        closest_pt = min([BP[0,0,:2],BP[-1,-1,:2]],key = lambda m:np.linalg.norm(m-C.flatten()))
        break
    segments = concatenate(DeCasteljau(segments),axis = 1)
    distances = concatenate(DeCasteljau(distances),axis = 1)
else:
    raise RuntimeError('max iterations reached')


plt.figure()
plotsegs(BP,C)
plt.plot(*closest_pt,'ko')
plotsegs(segments,C,'k')
