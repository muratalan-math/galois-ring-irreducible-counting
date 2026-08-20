#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <unordered_set>
using i64=long long;
static int pw(int p,int e){int r=1;while(e-->0)r*=p;return r;}
static int md(i64 x,int m){x%=m;if(x<0)x+=m;return (int)x;}
struct K{int t,x,y;bool operator==(K const&o)const{return t==o.t&&x==o.x&&y==o.y;}};
struct H{size_t operator()(K const&k)const noexcept{return ((size_t)k.t*1000003u+k.x)*1000033u+k.y;}};
long long countL(int p,int n,int v,int t0=1){int mt=pw(p,n),mx=pw(p,n+v);std::unordered_set<K,H>s;s.reserve(100000);
 for(int t=t0;t<mt;t+=p)for(int A=0;A<mx;A+=p)for(int B=0;B<mx;++B){int x=md(2LL*A*B,mx);int y=md(1LL*t*A*A-1LL*p*B*B,mt);s.insert({t,x,y});}
 return s.size();}
long long Lformula(int p,int n){long long L=1;for(int j=2;j<=n;++j){long long q3=1;for(int z=0;z<3*(j-2);++z)q3*=p;L=1LL*p*p*(p-1)/2*q3+p*L;}return L;}
int main(){for(auto cfg: {std::pair<int,int>{5,3},std::pair<int,int>{7,2}}){int p=cfg.first,N=cfg.second;for(int n=1;n<=N;++n)for(int v=0;v<=1;++v){long long got=countL(p,n,v),exp=Lformula(p,n);for(int i=0;i<v;++i)exp*=p;std::cout<<"p="<<p<<" n="<<n<<" v="<<v<<" got="<<got<<" expected="<<exp<<"\n";if(got!=exp)return 2;}}
 std::cout<<"ALL REFINED NONZERO-CURVE CHECKS PASSED\n";}
