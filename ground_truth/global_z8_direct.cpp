#include <bits/stdc++.h>
using namespace std;
int md(long long x){x%=8;if(x<0)x+=8;return (int)x;}
uint64_t enc4(int c3,int c2,int c1,int c0){return (((uint64_t)c3*8+c2)*8+c1)*8+c0;}
uint64_t enc3(int c2,int c1,int c0){return ((uint64_t)c2*8+c1)*8+c0;}
uint64_t enc2(int c1,int c0){return (uint64_t)c1*8+c0;}
int main(){
 cout<<"I_2,3(1)=8\n";
 {vector<uint8_t> r(64);for(int x=0;x<8;x++)for(int c1=0;c1<8;c1++){int c0=md(-1LL*x*x-c1*x);r[enc2(c1,c0)]=1;}uint64_t z=0;for(auto q:r)z+=q;cout<<"I_2,3(2)="<<64-z<<"\n";}
 {vector<uint8_t> r(512);for(int x=0;x<8;x++)for(int c2=0;c2<8;c2++)for(int c1=0;c1<8;c1++){int c0=md(-1LL*x*x*x-1LL*c2*x*x-c1*x);r[enc3(c2,c1,c0)]=1;}uint64_t z=0;for(auto q:r)z+=q;cout<<"I_2,3(3)="<<512-z<<"\n";}
 {vector<uint8_t> r(4096);for(int x=0;x<8;x++)for(int c3=0;c3<8;c3++)for(int c2=0;c2<8;c2++)for(int c1=0;c1<8;c1++){int c0=md(-1LL*x*x*x*x-1LL*c3*x*x*x-1LL*c2*x*x-c1*x);r[enc4(c3,c2,c1,c0)]=1;}
  for(int a=0;a<8;a++)for(int b=0;b<8;b++)for(int c=0;c<8;c++)for(int d=0;d<8;d++) r[enc4(md(a+c),md(b+d+a*c),md(a*d+b*c),md(b*d))]=1;
  uint64_t z=0;for(auto q:r)z+=q;cout<<"I_2,3(4)="<<4096-z<<"\n";}
}
