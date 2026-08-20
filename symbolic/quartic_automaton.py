#!/usr/bin/env python3
"""
End-to-end driver for the tame quartic local count J_{Q,s}(4), p>=5.
Implements Section 10 of the manuscript in two algebraically equivalent implementations:
  (1) the scalar finite-state recurrences (Section 10.2 pseudocode);
  (2) the fixed 9x9 root matrix T_Q (eq. TQ) and 12x12 factor companion
      matrix C_Q (eq. companion) via binary powering.
All arithmetic is performed exactly with fractions.Fraction.
No manuscript proof depends on this file.
"""
from fractions import Fraction as Fr

def sequences(Qint, SMAX):
    Q=Fr(Qint)
    io={1:Fr(0),2:Q*(Q-1),3:Q**2*(Q-1)*(2*Q+1)/2}
    for s in range(4,SMAX+2): io[s]=Q**2*io[s-1]+Q**2*io[s-2]-Q**4*io[s-3]
    nu={s:io[s]/Q**(s-1) for s in io}
    MQ2=(Q*Q-Q)/2
    Gam={0:Fr(0)}
    for s in range(1,SMAX+2): Gam[s]=Q*io[s]+MQ2*Q**(2*(s-1))
    alpha={s:io[s]-Q**2*io[s-1] for s in range(2,SMAX+2)}
    Lam={1:Fr(1)}
    for n in range(2,SMAX+2): Lam[n]=Q**2*(Q-1)/2*Q**(3*(n-2))+Q*Lam[n-1]
    def SqQ(r):
        if r==0: return Fr(1)
        return 1+(Q-1)/2*sum(Q**(r-2*j-1) for j in range(0,(r-1)//2+1))
    def SqQ2(r):
        if r==0: return Fr(1)
        return 1+(Q*Q-1)/2*sum(Q**(2*(r-2*j-1)) for j in range(0,(r-1)//2+1))
    def tildeSqQ2(r):
        b=SqQ2(r); return Q*b if r%2==0 else Q*(b-1)+1
    om={1:Fr(1)}
    for k in range(2,SMAX+2): om[k]=tildeSqQ2(k-2)
    return dict(Q=Q,io=io,nu=nu,MQ2=MQ2,Gam=Gam,alpha=alpha,Lam=Lam,SqQ=SqQ,om=om)

def U_scalar(Qint,SMAX,seq):
    Q=seq['Q']; io,alpha,MQ2,SqQ=seq['io'],seq['alpha'],seq['MQ2'],seq['SqQ']
    def A(s): return (Q*Q-Q+1)*(Q**4-1)*Q**(4*s-10)-(Q-1)*Q**(3*s-5)*SqQ(s-3)
    def E(s): return Q*MQ2*Q**(2*(s-1))*alpha[s]+(Q*(Q-1)/2)*(io[s]**2-Q**4*io[s-1]**2)
    MQ4=(Q**4-Q**2)/4
    H1=MQ4+((MQ2+1)*MQ2/2)
    U={1:Fr(0),2:(Q-1)*Q**3,3:(Q-1)*Q**5*(Q*Q+1),4:A(4),5:A(5)+Q**6*H1}
    F={}
    for s in range(6,SMAX+2):
        F[s]=A(s)-Q**4*A(s-1)+Q**6*E(s-4)
        U[s]=Q**4*U[s-1]+Q**7*U[s-4]-Q**11*U[s-5]+F[s]
    return U,F,A,E

def D_scalar(Qint,SMAX,seq):
    Q=seq['Q']; io,nu,Gam,Lam,om=seq['io'],seq['nu'],seq['Gam'],seq['Lam'],seq['om']
    CneqII=(Q-1)*(Q*(Q-1)-2)/8; CeqII=(Q-1)/2
    CID=Q*(Q-1)/2; CneqDD=(Q-1)/2
    def reg(s):
        return (Q*(Q-1)*(Q-2)/2*Q**(3*s-8)+(Q-1)*Q*Lam[s-2]+(Q-1)*Q**(s-1)*Gam[s-3])
    def Ffac(s):
        return (reg(s)+CneqII*Q**(3*(s-4))+CeqII*Q**(s-3)*om[s-4]
                +CID*Q**(s-2)*io[s-4]+CneqDD*Q**(s-2)*(nu[s-4]**2))
    D={0:Fr(0),1:Fr(0),2:Fr(0)}
    D[3]=reg(3); D[4]=reg(4)
    for s in range(5,SMAX+2):
        D[s]=Q**3*D[s-4]+Ffac(s)
    return D,Ffac

# ---- matrix path (integer matrices, exact-Fraction initial vectors) ----
def matmul(A,B):
    n,m,p=len(A),len(B),len(B[0])
    return [[sum(A[i][k]*B[k][j] for k in range(m)) for j in range(p)] for i in range(n)]
def matvec(A,v): return [sum(A[i][k]*v[k] for k in range(len(v))) for i in range(len(A))]
def matpow(A,e):
    n=len(A); R=[[Fr(1) if i==j else Fr(0) for j in range(n)] for i in range(n)]
    while e:
        if e&1: R=matmul(R,A)
        A=matmul(A,A); e>>=1
    return R
def cvec(Q):
    return [Fr(1),-(Q**3+Q),-(Q**6+Q**2),Q**9+2*Q**7+2*Q**5+Q**3,-Q**3,
            -(Q**13+2*Q**11+2*Q**9+Q**7)+Q**6+Q**4,Q**14+Q**10+Q**9+Q**5,
            Q**15+Q**13-Q**12-2*Q**10-2*Q**8-Q**6,-Q**16,
            Q**16+2*Q**14+2*Q**12+Q**10,-(Q**17+Q**13),-(Q**18+Q**16),Q**19]
def CQ_matrix(Q):
    c=cvec(Q); top=[-c[j] for j in range(1,13)]
    return [top]+[[Fr(1) if j==i else Fr(0) for j in range(12)] for i in range(11)]
def TQ_matrix(Q):
    cQ=Q**6+Q**4; dQ=Q**10
    Z=Fr(0); O=Fr(1)
    return [[Q**4,Z,Z,Q**7,-Q**11,Z,cQ,Z,-dQ],[O,Z,Z,Z,Z,Z,Z,Z,Z],
            [Z,O,Z,Z,Z,Z,Z,Z,Z],[Z,Z,O,Z,Z,Z,Z,Z,Z],[Z,Z,Z,O,Z,Z,Z,Z,Z],
            [Z,Z,Z,Z,Z,Z,cQ,Z,-dQ],[Z,Z,Z,Z,Z,O,Z,Z,Z],[Z,Z,Z,Z,Z,Z,O,Z,Z],
            [Z,Z,Z,Z,Z,Z,Z,O,Z]]

def _asint(x):
    x=Fr(x); assert x.denominator==1, f"non-integer {x}"; return x.numerator

def J_matrix_path(Qint,s):
    SMAX=max(s,14); Q=Fr(Qint); seq=sequences(Qint,SMAX)
    U,F,_,_=U_scalar(Qint,SMAX,seq); D,_=D_scalar(Qint,SMAX,seq)
    if s>=9:
        s9=[U[9],U[8],U[7],U[6],U[5],F[9],F[8],F[7],F[6]]
        Us=matvec(matpow(TQ_matrix(Q),s-9),s9)[0]
    else: Us=U[s]
    if s>=12:
        d12=[D[12-i] for i in range(12)]
        Ds=matvec(matpow(CQ_matrix(Q),s-12),d12)[0]
    else: Ds=D[s]
    return Us-Q**(s-1)*Ds, Us, Ds

def J_scalar_path(Qint,s):
    SMAX=max(s,6); Q=Fr(Qint); seq=sequences(Qint,SMAX)
    U,_,_,_=U_scalar(Qint,SMAX,seq); D,_=D_scalar(Qint,SMAX,seq)
    return U[s]-Q**(s-1)*D[s], U[s], D[s]

def Jp3_quartic(Qint):
    Q=Fr(Qint); r=2
    theta=Fr(Q+2,1)/2*Q**(r-1)-1-(Q-1)/2*sum(Q**(r-2*j-1) for j in range(1,(r-1)//2+1))
    return _asint((Q-1)*Q**(4*r-1)+(Q-1)*Q**(2*r-1)*theta)

if __name__=="__main__":
    for Q in (5,7,25):
        for s in range(2,41):
            js=J_scalar_path(Q,s); jm=J_matrix_path(Q,s)
            assert js==jm,(Q,s,js,jm)
    print("scalar-path == matrix-path for Q in {5,7,25}, s in [2,40]: PASS")
    GT={(5,2):(500,0,500),(5,3):(325000,170,320750),(5,4):(203812500,24850,200706250),
        (7,2):(2058,0,2058),(7,3):(5042100,777,5004027),(25,2):(375000,0,375000)}
    for (Q,s),(Uv,Dv,Jv) in GT.items():
        js,Us,Ds=J_scalar_path(Q,s)
        assert (_asint(Us),_asint(Ds),_asint(js))==(Uv,Dv,Jv),((Q,s),(Us,Ds,js),(Uv,Dv,Jv))
    print("driver (U_s,D_s,J) == brute-force ground truth (incl. new s=4): PASS")
    for Q in (5,7,25):
        assert _asint(J_scalar_path(Q,3)[0])==Jp3_quartic(Q)
    print("driver J_{Q,3}(4) == characteristic-p^3 closed form: PASS")
    print("\nQ  s   J_{Q,s}(4)")
    for s in (2,3,4,5,8,12,13,20):
        print(f"5  {s:<3} {_asint(J_matrix_path(5,s)[0])}")
