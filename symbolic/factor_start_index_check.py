#!/usr/bin/env python3
"""Symbolic check of the quartic factor recurrence start index.

This script reconstructs the exact factor sequence D_s from the displayed
factor decomposition and origin formula, without using the order-12
recurrence. It then convolves D_s with the announced denominator.
"""
import sympy as sp

Q,T=sp.symbols('Q T')

# Degree-two input.
iota={1:sp.Integer(0),2:Q*(Q-1),3:Q**2*(Q-1)*(2*Q+1)/2}
for s in range(4,25):
    iota[s]=sp.expand(Q**2*iota[s-1]+Q**2*iota[s-2]-Q**4*iota[s-3])
nu={s:sp.cancel(iota[s]/Q**(s-1)) for s in range(1,25)}
M2=(Q**2-Q)/2
Gamma={0:sp.Integer(0)}
for s in range(1,25):
    Gamma[s]=sp.expand(Q*iota[s]+M2*Q**(2*(s-1)))

# Nonzero repeated-curve state.
Lambda={1:sp.Integer(1)}
for n in range(2,25):
    Lambda[n]=sp.expand(Q**2*(Q-1)*Q**(3*(n-2))/2+Q*Lambda[n-1])

# Repeated-basic square state.
def square_Q2(r):
    if r==0:
        return sp.Integer(1)
    return sp.expand(1+(Q**2-1)*sum(Q**(2*(r-2*j-1)) for j in range(0,(r-1)//2+1))/2)

def weighted_square(r):
    S=square_Q2(r)
    return sp.expand(Q*S if r%2==0 else Q*(S-1)+1)

omega={1:sp.Integer(1)}
for k in range(2,25):
    omega[k]=weighted_square(k-2)

CIIneq=(Q-1)*(Q*(Q-1)-2)/8
CIIeq=(Q-1)/2
CID=Q*(Q-1)/2
CDDneq=(Q-1)/2

D={0:sp.Integer(0),1:sp.Integer(0),2:sp.Integer(0)}

def origin(s):
    if s<=4:
        return sp.Integer(0)
    k=s-4
    return sp.expand(
        CIIneq*Q**(3*k)
        + CIIeq*Q**(k+1)*omega[k]
        + CID*Q**(k+2)*iota[k]
        + CDDneq*Q**(k+2)*nu[k]**2
        + Q**3*D[k]
    )

def regular(s):
    return Q*(Q-1)*(Q-2)*Q**(3*s-8)/2

# Reconstruct from the exact factor decomposition, not from the recurrence.
for s in range(3,25):
    D[s]=sp.factor(
        regular(s)
        +(Q-1)*Q*Lambda[s-2]
        +(Q-1)*Q**(s-1)*Gamma[s-3]
        +origin(s)
    )

den=sp.expand(
    (1-Q**3*T**4)*(1-Q*T)*(1-Q**3*T)
    *(1-Q**2*T**2)*(1-Q**4*T**2)*(1-Q**6*T**2)
)
c=[sp.expand(den.coeff(T,j)) for j in range(13)]

def convolution(s):
    return sp.factor(sum(c[j]*D.get(s-j,0) for j in range(13) if s-j>=0))

expected_T12=Q**16*(Q-1)*(2*Q**2-Q+1)/8
assert sp.factor(convolution(12)-expected_T12)==0
for s in range(13,25):
    assert convolution(s)==0

print('factor start-index symbolic check: PASS')
print('[T^12](D_D(T) G_D(T)) =', sp.factor(convolution(12)))
print('convolution vanishes for s=13,...,24 reconstructed from exact decomposition')
