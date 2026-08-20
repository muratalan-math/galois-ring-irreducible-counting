#include <bits/stdc++.h>
using namespace std;

struct E { int a,b; };
struct Ring {
    int p,s,m,mod,rel0,rel1,base;
    Ring(int P,int S,int M):p(P),s(S),m(M){
        mod=1; for(int i=0;i<s;i++) mod*=p;
        base=mod/p;
        if(m==1){rel0=0; rel1=0;}
        else if(p==2){ rel0=-1; rel1=-1; } // t^2+t+1=0
        else if(p==5){ rel0=2; rel1=0; }   // t^2-2 irreducible mod 5
        else { cerr<<"unsupported quadratic residue field\n"; exit(2); }
    }
    int norm(long long x) const { x%=mod; if(x<0)x+=mod; return (int)x; }
    E add(E x,E y)const{return{norm(x.a+y.a),norm(x.b+y.b)};}
    E neg(E x)const{return{norm(-x.a),norm(-x.b)};}
    E mul(E x,E y)const{
        if(m==1) return {norm(1LL*x.a*y.a),0};
        long long a=1LL*x.a*y.a + 1LL*x.b*y.b*rel0;
        long long b=1LL*x.a*y.b + 1LL*x.b*y.a + 1LL*x.b*y.b*rel1;
        return {norm(a),norm(b)};
    }
    E powe(E x,int n)const{E r{1,0}; while(n){if(n&1)r=mul(r,x); x=mul(x,x); n>>=1;} return r;}
    vector<E> pIdeal() const {
        vector<E> v;
        if(m==1){ for(int a0=0;a0<base;a0++) v.push_back({p*a0,0}); }
        else { for(int b0=0;b0<base;b0++) for(int a0=0;a0<base;a0++) v.push_back({p*a0,p*b0}); }
        return v;
    }
    int idxP(E x) const {
        if(x.a%p || x.b%p){ cerr<<"not in pR\n"; exit(3); }
        int aa=(x.a/p)%base;
        if(m==1) return aa;
        int bb=(x.b/p)%base;
        return aa + base*bb;
    }
    long long Q()const{ long long q=p; for(int i=1;i<m;i++)q*=p; return q; }
};

static inline uint64_t enc4(int a3,int a2,int a1,int a0,int N){
    return (((uint64_t)a3*N+a2)*N+a1)*N+a0;
}

uint64_t count_irred_root_only(const Ring& R,int deg){
    auto P=R.pIdeal(); int N=P.size();
    uint64_t total=1; for(int i=0;i<deg;i++) total*=N;
    vector<uint8_t> red(total,0);
    if(deg==2){
        for(auto x:P){ E x2=R.mul(x,x); for(int i1=0;i1<N;i1++){
            E c0=R.neg(R.add(x2,R.mul(P[i1],x)));
            uint64_t id=(uint64_t)i1*N+R.idxP(c0); red[id]=1;
        }}
    } else if(deg==3){
        for(auto x:P){ E x2=R.mul(x,x),x3=R.mul(x2,x); for(int i2=0;i2<N;i2++) for(int i1=0;i1<N;i1++){
            E val=R.add(x3,R.add(R.mul(P[i2],x2),R.mul(P[i1],x)));
            E c0=R.neg(val);
            uint64_t id=((uint64_t)i2*N+i1)*N+R.idxP(c0); red[id]=1;
        }}
    }
    uint64_t rcount=0; for(auto z:red) rcount+=z;
    return total-rcount;
}

uint64_t count_irred_quartic(const Ring& R, uint64_t *root_count=nullptr, uint64_t *qfac_count=nullptr){
    auto P=R.pIdeal(); int N=P.size();
    uint64_t total=(uint64_t)N*N*N*N;
    vector<uint8_t> red(total,0), roots(total,0), qfac(total,0);
    // Root-bearing quartics: choose x,c3,c2,c1, solve c0.
    for(auto x:P){ E x2=R.mul(x,x),x3=R.mul(x2,x),x4=R.mul(x2,x2);
        for(int i3=0;i3<N;i3++) for(int i2=0;i2<N;i2++) for(int i1=0;i1<N;i1++){
            E val=x4;
            val=R.add(val,R.mul(P[i3],x3));
            val=R.add(val,R.mul(P[i2],x2));
            val=R.add(val,R.mul(P[i1],x));
            int i0=R.idxP(R.neg(val));
            roots[enc4(i3,i2,i1,i0,N)]=1;
        }
    }
    // Products of two monic local quadratics.
    for(int iu=0;iu<N;iu++) for(int iv=0;iv<N;iv++)
      for(int ju=0;ju<N;ju++) for(int jv=0;jv<N;jv++){
        E u=P[iu],v=P[iv], U=P[ju],V=P[jv];
        E c3=R.add(u,U);
        E c2=R.add(R.add(v,V),R.mul(u,U));
        E c1=R.add(R.mul(u,V),R.mul(U,v));
        E c0=R.mul(v,V);
        qfac[enc4(R.idxP(c3),R.idxP(c2),R.idxP(c1),R.idxP(c0),N)]=1;
      }
    uint64_t nr=0,nq=0,nred=0;
    for(uint64_t i=0;i<total;i++){nr+=roots[i]; nq+=qfac[i]; red[i]=roots[i]|qfac[i]; nred+=red[i];}
    if(root_count)*root_count=nr; if(qfac_count)*qfac_count=nq;
    return total-nred;
}

int main(int argc,char**argv){
    if(argc!=5){cerr<<"usage: direct_low_degree p s m degree\n"; return 2;}
    int p=stoi(argv[1]),s=stoi(argv[2]),m=stoi(argv[3]),deg=stoi(argv[4]);
    Ring R(p,s,m); long long Q=R.Q();
    uint64_t ans;
    if(deg<=3) ans=count_irred_root_only(R,deg);
    else if(deg==4){uint64_t nr,nq;ans=count_irred_quartic(R,&nr,&nq);
        auto P=R.pIdeal(); uint64_t N=P.size(); uint64_t total=N*N*N*N;
        uint64_t rootfree=total-nr;            // U_s
        uint64_t union_rc=total-ans;           // |roots ∪ qfac|
        uint64_t rootfree_qfac=union_rc-nr;    // B_s = Q^{s-1} D_s = |qfac \ roots|
        cout<<"root_bearing="<<nr<<" quadratic_factor="<<nq
            <<" rootfree_U="<<rootfree<<" rootfree_qfac_B="<<rootfree_qfac<<"\n";}
    else {cerr<<"degree <=4 only\n"; return 2;}
    cout<<"p="<<p<<" s="<<s<<" m="<<m<<" Q="<<Q<<" degree="<<deg<<" J_direct="<<ans<<"\n";
}
