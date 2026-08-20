#!/usr/bin/env python3
# Independent finite-level and boundary checks used by the manuscript supplement.

# 0) Valuation min-gap sublemma (unifies checks 1 and 2 below):
#    for c in {0,1}, mu,nu <= lam:  min(c+lam,2nu) - min(2mu,c+2nu) >= -lam.
for c in (0,1):
    for lam in range(0,40):
        for mu in range(0,lam+1):
            for nu in range(0,lam+1):
                assert min(c+lam,2*nu)-min(2*mu,c+2*nu) >= -lam, (c,lam,mu,nu)
# reduction steps used to derive checks 1 and 2 from the sublemma:
for lam in range(0,40):
    for x in range(0,lam+1):
        assert min(1+lam,2*x) >= min(lam,2*x)          # used for (8.19)
        assert min(lam,2*x)   >= min(1+lam,2*x)-1       # used for (8.5)
print('valuation min-gap sublemma (+reductions): PASS')

# 1) Ramified congruence minimum inequality.
for va in range(9):
    for delta in range(9):
        lam=va+delta
        for beta in range(lam+1):
            ve=lam-beta
            lhs=min(lam,2*beta)-min(2*delta,1+2*beta)
            assert lhs >= -1-lam, (va,delta,beta,ve,lam,lhs)
print('ramified minimum inequality: PASS')

# 2) Repeated-basic normalization minimum inequality.
for lam in range(9):
    for alpha in range(lam+1):
        for beta in range(lam+1):
            lhs=min(1+lam,2*alpha)-min(2*beta,2*alpha)
            assert lhs >= -lam, (lam,alpha,beta,lhs)
print('repeated-basic minimum inequality: PASS')

# 3) Exact Q^2 overlap fiber at the boundary s=3 over Z/5^3.
p=5; s=3; mod=p**s
# Sources: a,lambda mod p^(s-1)=25; c mod p^(s-2)=5; lambda unit.
# Targets: P mod25, L mod5, R mod5. R is always 0 at s=3.
from collections import Counter
fib=Counter()
for a in range(p**(s-1)):
    for la in range(p**(s-1)):
        if la%p==0: continue
        for c in range(p**(s-2)):
            P=(la+p*(c-a*a))%(p**(s-1))
            L=(a*(la-p*c))%(p**(s-2))
            R=(p*c*la)%(p**(s-2))
            fib[(P,L,R)] += 1
vals=set(fib.values())
assert vals=={p**2}, vals
print(f'overlap s=3 exact oriented fiber: PASS ({len(fib)} targets, fiber={p**2})')
