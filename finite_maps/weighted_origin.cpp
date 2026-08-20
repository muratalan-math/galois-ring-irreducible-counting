#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <vector>
using i64=long long;
static int ipow(int p,int e){int r=1;while(e-->0)r*=p;return r;}
static int md(i64 x,int m){x%=m;if(x<0)x+=m;return (int)x;}
static int qtype(int a,int b,int p,int&root){int ct=0;root=-1;for(int x=0;x<p;++x)if(md((i64)x*x+(i64)a*x+b,p)==0){++ct;root=x;}if(ct==0)return 0;if(ct==1)return 1;return 2;}
struct Key{int P,L,R;bool operator==(Key const&o)const{return P==o.P&&L==o.L&&R==o.R;}};
struct Hash{size_t operator()(Key const&k)const noexcept{return ((size_t)k.P*1000003u+k.L)*1000033u+k.R;}};
static std::array<long long,5> run(int p,int k){
 int mP=ipow(p,k+2),mL=ipow(p,k+1),mR=ipow(p,k),mPoly=ipow(p,k);
 // E modulo p^(k+1) is enough; a,M modulo p^(k+2).
 std::unordered_map<Key,unsigned char,Hash> masks;
 masks.reserve(k==1?5000:400000);
 for(int a=0;a<mP;++a){int ar=a%p;
  for(int M=0;M<mP;++M){
   int P=md(2LL*M-1LL*a*a,mP);
   for(int E=0;E<mL;++E){
    int L=md(2LL*a*E,mL); int R=md(1LL*M*M-1LL*E*E,mR);
    int b=md(M-E,p), c=md(M+E,p), r1=-1,r2=-1;
    int t1=qtype(ar,b,p,r1),t2=qtype(md(-ar,p),c,p,r2),cat=7;
    if(t1==0&&t2==0)cat=(ar==0&&b==c)?1:0;
    else if((t1==0&&t2==1)||(t1==1&&t2==0))cat=2;
    else if(t1==1&&t2==1)cat=(r1==r2)?4:3;
    masks[{P,L,R}]|=(unsigned char)(1u<<cat);
   }
  }
 }
 std::array<long long,5> out{};
 for(auto const&e:masks){auto kx=e.first;bool root=false;int P=kx.P%mPoly,L=kx.L%mPoly,R=kx.R%mPoly;
   for(int z=0;z<mPoly;++z){i64 z2=1LL*z*z%mPoly; if((z2*z2+1LL*P*z2+1LL*L*z+R)%mPoly==0){root=true;break;}}
   if(root)continue; unsigned m=e.second; if(__builtin_popcount((unsigned)m)!=1){std::cerr<<"ambiguity\n";std::exit(2);} int cat=__builtin_ctz((unsigned)m); if(cat<5)++out[cat];
 }
 return out;
}
int main(){const char*nm[]={"II-ne","II-eq","ID","DD-ne","DD-same"};
 std::array<long long,5> e1={1125,50,0,0,0},e2={140625,1250,125000,20000,0};
 for(int k: {1,2}){auto a=run(5,k);std::cout<<"weighted origin Q=5 k="<<k<<"\n";for(int i=0;i<5;++i)std::cout<<"  "<<nm[i]<<": "<<a[i]<<"\n";if(a!=(k==1?e1:e2))return 3;}
 std::cout<<"ALL WEIGHTED-ORIGIN ENUMERATIONS PASSED\n";
}
