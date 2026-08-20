#!/usr/bin/env python3
import sympy as sp
Q,T=sp.symbols('Q T')
D=sp.expand((1-Q**3*T**4)*(1-Q*T)*(1-Q**3*T)*(1-Q**2*T**2)*(1-Q**4*T**2)*(1-Q**6*T**2))
coeff=[sp.expand(D.coeff(T,j)) for j in range(13)]
expected=[
1,
-(Q**3+Q),
-(Q**6+Q**2),
Q**9+2*Q**7+2*Q**5+Q**3,
-Q**3,
-Q**13-2*Q**11-2*Q**9-Q**7+Q**6+Q**4,
Q**14+Q**10+Q**9+Q**5,
Q**15+Q**13-Q**12-2*Q**10-2*Q**8-Q**6,
-Q**16,
Q**16+2*Q**14+2*Q**12+Q**10,
-Q**17-Q**13,
-Q**18-Q**16,
Q**19]
assert all(sp.expand(a-b)==0 for a,b in zip(coeff,expected))
print('factor denominator expansion: PASS')
for j,c in enumerate(coeff): print(f'c_{j} = {c}')
