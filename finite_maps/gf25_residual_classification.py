#!/usr/bin/env python3
from collections import defaultdict

# GF(25)=F_5[t]/(t^2-2), since 2 is a nonsquare modulo 5.
# (Same nonsquare choice as ground_truth/local_direct_enumeration.cpp;
#  the residual-type multiplicities are independent of this choice.)
P=5

def dec(x): return (x%P, x//P)
def enc(a,b): return (a%P)+P*(b%P)
def add(x,y):
    a,b=dec(x); c,d=dec(y); return enc(a+c,b+d)
def neg(x):
    a,b=dec(x); return enc(-a,-b)
def sub(x,y): return add(x,neg(y))
def mul(x,y):
    a,b=dec(x); c,d=dec(y)
    return enc(a*c+2*b*d,a*d+b*c)
def sq(x): return mul(x,x)

Q=25
els=list(range(Q))
zero=0

def qtype(a,b):
    roots=[x for x in els if add(add(sq(x),mul(a,x)),b)==0]
    if len(roots)==0: return ('I',None)
    if len(roots)==1: return ('D',roots[0])
    return ('S',None)

qt={(a,b):qtype(a,b) for a in els for b in els}
sets=[set() for _ in range(5)]
ambiguous={}
for a in els:
    na=neg(a)
    a2=sq(a)
    for b in els:
        t1,r1=qt[(a,b)]
        for c in els:
            t2,r2=qt[(na,c)]
            cat=None
            if t1=='I' and t2=='I':
                cat=1 if (a==0 and b==c) else 0
            elif {t1,t2}=={'I','D'}:
                cat=2
            elif t1=='D' and t2=='D':
                cat=4 if r1==r2 else 3
            if cat is None: continue
            target=(sub(add(b,c),a2),mul(a,sub(c,b)),mul(b,c))
            for j,s in enumerate(sets):
                if target in s and j!=cat:
                    ambiguous[target]=(j,cat)
            sets[cat].add(target)

h=(Q-1)//2
expected=[(Q-1)*(Q*(Q-1)-2)//8,h,Q*(Q-1)//2,h,1]
actual=[len(s) for s in sets]
print('GF(25) residual product counts:', actual)
print('GF(25) expected counts:', expected)
assert actual==expected and not ambiguous
print('GF(25) FIVE-TYPE CLASSIFICATION PASSED')
