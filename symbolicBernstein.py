# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 15:41:50 2021

@author: Phillip Dix
"""

from sympy import *

def BernsteinBasis(n,v,x):
    return binomial(n,v)*x**v*(1-x)**(n-v)

Px = var('P0x P1x P2x P3x')
Py = var('P0y P1y P2y P3y')

T = var('T0 T1 T2 T3 T4 T5 T6')

var('t Cx Cy')

cubicx = sum([BernsteinBasis(3,i,t)*j for i,j in enumerate(Px)])
cubicy = sum([BernsteinBasis(3,i,t)*j for i,j in enumerate(Py)])

E1 = Eq((cubicx - Cx)**2 + (cubicy - Cy)**2,sum([BernsteinBasis(6,i,t)*j for i,j in enumerate(T)]))

newpoints = solve(E1,T)

pts = [collect(expand(i),t) for i in newpoints.values()]

a,b,c,d,e,f = symbols('a b c d e f', Wild = True)

M = [poly(i,[Cx,Cy]) for i in pts]

M1 = [(i.coeff_monomial(Cx**2),i.coeff_monomial(Cx),
       i.coeff_monomial(Cy**2),i.coeff_monomial(Cy),j.subs(((Cx,0),(Cy,0)))) for i,j in zip(M,pts)]

Base = Matrix([Cx**2,Cx,Cy**2,Cy,1]).transpose()

TestPts1 = Base * Matrix(M1).transpose()

testvals = [simplify(i-j) for i,j in zip(pts,TestPts1)]

assert all([i == 0 for i in testvals])

print(pycode(newpoints))