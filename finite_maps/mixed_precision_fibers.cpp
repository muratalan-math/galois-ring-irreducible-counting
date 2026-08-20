#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using i64 = long long;
static int ipow(int p,int e){int r=1; while(e-->0) r*=p; return r;}
static int modn(i64 x,int m){x%=m; if(x<0)x+=m; return (int)x;}
static int invmod(int a,int m){int t=0,nt=1,r=m,nr=a;while(nr){int q=r/nr;int tt=t-q*nt;t=nt;nt=tt;int rr=r-q*nr;r=nr;nr=rr;}if(r!=1)std::abort();return modn(t,m);}
static i64 key3(int x,int mx,int y,int my,int z){return (i64(x)*my+y)*1'000'000LL+z;}
static i64 key2(int x,int y,int my){return i64(x)*my+y;}

// Fixed coprime residual pair over F_5:
// Z^2+2 (irreducible) and Z^2 (double), represented by a=0,b=2,c=0.
static void coprime_projection_check(int v){
    const int p=5,k=1,N=k+2+v;
    const int mod=ipow(p,N);
    const int modP=ipow(p,k+2), modL=ipow(p,k+1+v), modR=ipow(p,k);
    std::unordered_map<i64,int> hist;
    hist.reserve(ipow(p,3+v)*2);
    const int lifts=ipow(p,N-1);
    for(int A=0;A<lifts;++A){int a=p*A;
      for(int B=0;B<lifts;++B){int b=2+p*B;
        for(int C=0;C<lifts;++C){int c=p*C;
          int P=modn((i64)b+c-(i64)a*a,modP);
          int L=modn((i64)a*(c-b),modL);
          int R=modn((i64)b*c,modR);
          i64 key=(i64(P)*modL+L)*modR+R;
          ++hist[key];
        }
      }
    }
    const int expectedTargets=ipow(p,3+v);
    const int expectedFiber=ipow(p,3+2*v);
    std::map<int,int> fh;
    for(auto &e:hist) ++fh[e.second];
    std::cout<<"coprime v="<<v<<": targets="<<hist.size()<<", fiber histogram";
    for(auto &e:fh) std::cout<<" "<<e.first<<":"<<e.second;
    std::cout<<"\n";
    if((int)hist.size()!=expectedTargets || fh.size()!=1 || fh.begin()->first!=expectedFiber) std::exit(2);
}

struct PairHash { size_t operator()(const std::uint64_t &x) const noexcept {return std::hash<std::uint64_t>{}(x);} };

static std::unordered_set<std::uint64_t,PairHash> generalized_image(int p,int T,int e,int r,int v){
    int mx=ipow(p,r+1+v), my=ipow(p,r);
    int inv4 = my==1?0:invmod(4,my);
    int pe2 = my==1?0:ipow(p,2*e)%my;
    std::unordered_set<std::uint64_t,PairHash> out;
    out.reserve((size_t)mx*8);
    for(int A=0;A<mx;++A){
      i64 A2=(i64)A*A;
      i64 A4=A2*A2;
      for(int B=0;B<mx;++B){
        int x=modn((i64)A*B,mx);
        int y=0;
        if(my!=1) y=modn((i64)B*B-(i64)T*A2-(i64)pe2*inv4%my*(A4%my),my);
        out.insert((std::uint64_t)x*(std::uint64_t)my+(unsigned)y);
      }
    }
    return out;
}

static std::unordered_set<std::uint64_t,PairHash> square_image(int p,int T,int r){
    int mx=ipow(p,r+1), my=ipow(p,r);
    std::unordered_set<std::uint64_t,PairHash> out;
    out.reserve((size_t)mx*8);
    // theta^2=-T; trace-zero coefficient is 2AB, scalar is B^2-T A^2.
    for(int A=0;A<mx;++A){
      for(int B=0;B<mx;++B){
        int x=modn(2LL*A*B,mx);
        int y=my==1?0:modn((i64)B*B-(i64)T*A*A,my);
        out.insert((std::uint64_t)x*(std::uint64_t)my+(unsigned)y);
      }
    }
    return out;
}

static void repeated_check(int p,int T,int e,int r,int v){
    auto base=generalized_image(p,T,e,r,0);
    auto ref=generalized_image(p,T,e,r,v);
    int my=ipow(p,r), baseMx=ipow(p,r+1), refMx=ipow(p,r+1+v);
    std::unordered_map<std::uint64_t,int,PairHash> fibers;
    fibers.reserve(base.size()*2);
    for(auto key:ref){
      int y=(int)(key%my); int x=(int)(key/my);
      int xb=x%baseMx;
      std::uint64_t bk=(std::uint64_t)xb*my+y;
      ++fibers[bk];
    }
    std::map<int,int> fh; for(auto&e:fibers)++fh[e.second];
    int expected=ipow(p,v);
    std::cout<<"N e="<<e<<", r="<<r<<", v="<<v<<": base="<<base.size()<<", refined="<<ref.size()<<", fibers";
    for(auto&e:fh)std::cout<<" "<<e.first<<":"<<e.second; std::cout<<"\n";
    if(fibers.size()!=base.size() || fh.size()!=1 || fh.begin()->first!=expected) std::exit(3);
    if(v==0){
      auto sq=square_image(p,T,r);
      std::cout<<"  square normal-form image="<<sq.size()<<"\n";
      if(sq.size()!=base.size()) std::exit(4);
    }
}

int main(){
    coprime_projection_check(0);
    coprime_projection_check(1);
    // p=5, T=2 has -T=3 nonsquare.
    for(int r=0;r<=3;++r) repeated_check(5,2,1,r,0);
    for(int e: {1,2}) for(int r=0;r<=2;++r) for(int v: {1,2}) repeated_check(5,2,e,r,v);
    std::cout<<"ALL MIXED-PRECISION LOCAL MAP CHECKS PASSED\n";
}
