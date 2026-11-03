#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
using V=std::array<double,3>;
V add(V a,V b,double c=1){for(int k=0;k<3;k++)a[k]+=c*b[k];return a;}
V cross(V a,V b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
double dot(V a,V b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
double norm(V a){return std::sqrt(dot(a,a));}
struct Bond{int j;double a,c;};
int main(int argc,char**argv){try{
 if(argc<7)throw std::runtime_error("model mode dt_fs duration_ps sample_ps output [initial_binary] [perturb] [field_y_T sot_T on_ps off_ps scalar_sample_ps]");
 std::ifstream in(argv[1]);int N,W,L;double J,D,K,alpha,gamma;
 in>>N>>W>>L>>J>>D>>K>>alpha>>gamma; if(!in)throw std::runtime_error("model header");
 std::vector<V> e(N),m(N);for(int i=0;i<N;i++)for(auto*p:{&e[i],&m[i]})for(auto&v:*p)in>>v;
 std::string mode=argv[2];std::vector<std::vector<Bond>> bonds(N);int nb;in>>nb;
 for(int b=0;b<nb;b++){int i,j;double a,c;in>>i>>j>>a>>c;
  if(mode=="reciprocal" || mode=="triangle"){a=a==0?0:.5;c=c==0?0:std::copysign(.5,c);}
  else if(mode!="current")throw std::runtime_error("mode");
  bonds[i].push_back({j,a*J,c*D});
 }
 if(!in)throw std::runtime_error("model read");
 // K in the input is the field coefficient 2*Ksim/mu. Count occupied bonds.
 std::vector<double> anisotropy(N,K);
 if(mode=="triangle")for(int i=0;i<N;i++){
  int z=0;for(auto b:bonds[i])if(b.a!=0)z++;
  anisotropy[i]*=z/4.0;
 }
 if(argc>7 && std::string(argv[7])!="-"){std::ifstream init(argv[7],std::ios::binary);init.read((char*)m.data(),N*sizeof(V));if(!init)throw std::runtime_error("initial read");}
 double perturb=argc>8?std::stod(argv[8]):0;
 for(int i=0;i<N;i++){if(perturb){V p={std::sin(1.7*(i+1)),std::cos(2.3*(i+1)),std::sin(.9*(i+1))};m[i]=add(m[i],p,perturb);double n=norm(m[i]);for(auto&v:m[i])v/=n;}}
 const double dt=std::stod(argv[3])*1e-15;long long steps=llround(std::stod(argv[4])*1e-12/dt),sample=std::max(1LL,llround(std::stod(argv[5])*1e-12/dt));
 const bool driven=argc>9;if(driven && argc!=14)throw std::runtime_error("drive arguments");
 const double by=driven?std::stod(argv[9]):0,bt=driven?std::stod(argv[10]):0;
 const long long on=driven?llround(std::stod(argv[11])*1e-12/dt):0,off=driven?llround(std::stod(argv[12])*1e-12/dt):0;
 const long long scalar_sample=driven?std::max(1LL,llround(std::stod(argv[13])*1e-12/dt)):sample;
 bool active=driven && on<=0 && off>0;long long active_updates=0;
 std::vector<int> cls(N),degree(N);std::array<int,3> counts{};int edge_count=0,interior_count=0;
 for(int i=0;i<N;i++){
  cls[i]=e[i][1]<-.5?1:(e[i][2]<0?2:0);counts[cls[i]]++;
  for(auto b:bonds[i])if(b.a!=0)degree[i]++;
  if(degree[i]==2)edge_count++;else interior_count++;
 }
 const double scale=gamma/(1+alpha*alpha);double pre=0,post=0;
 auto field=[&](const std::vector<V>&s,std::vector<V>&h){
  for(int i=0;i<N;i++){h[i]={0,0,0};for(auto b:bonds[i]){for(int k=0;k<3;k++)h[i][k]-=b.a*s[b.j][k];h[i][1]+=b.c*s[b.j][2];h[i][2]-=b.c*s[b.j][1];}h[i]=add(h[i],e[i],anisotropy[i]*dot(s[i],e[i]));if(active)h[i][1]+=by;}
 };
 std::vector<V> h(N),k1(N),k2(N),k3(N),k4(N),tmp(N),old=m;
 auto rhs=[&](const std::vector<V>&s,std::vector<V>&k){field(s,h);for(int i=0;i<N;i++){V c=cross(s[i],h[i]),cc=cross(s[i],c);for(int j=0;j<3;j++)k[i][j]=-scale*(c[j]+alpha*cc[j]);
  if(active){V p=cross(s[i],{1,0,0}),pp=cross(s[i],p);for(int j=0;j<3;j++)k[i][j]+=scale*bt*(-pp[j]+alpha*p[j]);}
 }};
 std::string out=argv[6];std::ofstream csv(out+".csv"),bin(out+".bin",std::ios::binary);
 std::ofstream obs,times;if(driven){
  obs.open(out+"-observables.csv");times.open(out+"-snapshot-times.csv");obs<<std::setprecision(17);times<<std::setprecision(17)<<"time_ps\n";
  obs<<"time_ps,site_x,site_y,site_z,balanced_x,balanced_y,balanced_z";
  for(auto name:{"a","b","c","edge","interior"})for(auto axis:{"x","y","z"})obs<<','<<name<<'_'<<axis;
  obs<<",max_sublattice_spread,max_torque_T\n";
 }
 csv<<std::setprecision(17)<<"time_ps,max_torque_T,rms_torque_T,energy_over_mu_per_spin,drift_since_sample,max_abs_x,mean_x,mean_y,mean_z\n";
 double maxt=0,rms=0,energy=0,drift=0,maxx=0;V mean{};
 auto record=[&](long long step){field(m,h);maxt=rms=energy=drift=maxx=0;mean={0,0,0};
  for(int i=0;i<N;i++){double t=norm(cross(m[i],h[i]));maxt=std::max(maxt,t);rms+=t*t;drift=std::max(drift,norm(add(m[i],old[i],-1)));maxx=std::max(maxx,std::abs(m[i][0]));mean=add(mean,m[i],1.0/N);
   energy-=.5*anisotropy[i]*std::pow(dot(m[i],e[i]),2);
   if(active)energy-=by*m[i][1];
   for(auto b:bonds[i])energy+=.5*(b.a*dot(m[i],m[b.j])-b.c*cross(m[i],m[b.j])[0]);
  }
  rms=std::sqrt(rms/N);energy/=N;
  csv<<step*dt/1e-12<<','<<maxt<<','<<rms<<','<<energy<<','<<drift<<','<<maxx<<','<<mean[0]<<','<<mean[1]<<','<<mean[2]<<'\n';
  bin.write((char*)m.data(),N*sizeof(V));old=m;
  if(driven)times<<step*dt/1e-12<<'\n';
 };
 auto observe=[&](long long step){
  std::array<V,3> sub{};V site{},balanced{},edge{},interior{};
  for(int i=0;i<N;i++){sub[cls[i]]=add(sub[cls[i]],m[i],1.0/counts[cls[i]]);site=add(site,m[i],1.0/N);
   if(degree[i]==2)edge=add(edge,m[i],1.0/edge_count);else interior=add(interior,m[i],1.0/interior_count);}
  for(auto v:sub)balanced=add(balanced,v,1.0/3);
  if(interior_count==0)interior.fill(std::numeric_limits<double>::quiet_NaN());
  double spread=0,torque=0;field(m,h);
  for(int i=0;i<N;i++){spread=std::max(spread,norm(add(m[i],sub[cls[i]],-1)));torque=std::max(torque,norm(cross(m[i],h[i])));}
  obs<<step*dt/1e-12;for(auto v:{site,balanced,sub[0],sub[1],sub[2],edge,interior})for(auto x:v)obs<<','<<x;
  obs<<','<<spread<<','<<torque<<'\n';
 };
 record(0);
 if(driven)observe(0);
 for(long long s=1;s<=steps;s++){
  active=driven && s-1>=on && s-1<off;active_updates+=active;
  rhs(m,k1);for(int i=0;i<N;i++)tmp[i]=add(m[i],k1[i],dt/2);
  rhs(tmp,k2);for(int i=0;i<N;i++)tmp[i]=add(m[i],k2[i],dt/2);
  rhs(tmp,k3);for(int i=0;i<N;i++)tmp[i]=add(m[i],k3[i],dt);
  rhs(tmp,k4);
  for(int i=0;i<N;i++){for(int k=0;k<3;k++)m[i][k]+=dt*(k1[i][k]+2*k2[i][k]+2*k3[i][k]+k4[i][k])/6;double n=norm(m[i]);if(!std::isfinite(n))throw std::runtime_error("nonfinite");pre=std::max(pre,std::abs(n-1));for(auto&v:m[i])v/=n;post=std::max(post,std::abs(norm(m[i])-1));}
  active=driven && s>=on && s<off;
  // Resolve pulse edges and the first 10 ps after them without dense full-spin storage.
  bool event=driven && ((s>=on && s<=on+llround(10e-12/dt)) || (s>=off && s<=off+llround(10e-12/dt))) && s%scalar_sample==0;
  if(s%sample==0 || s==steps || event)record(s);
  if(driven && (s%scalar_sample==0 || s==steps))observe(s);
 }
 std::ofstream final(out+"-final.bin",std::ios::binary);final.write((char*)m.data(),N*sizeof(V));
 std::ofstream meta(out+".json");meta<<std::setprecision(17)<<"{\"N\":"<<N<<",\"dt_fs\":"<<dt/1e-15<<",\"steps\":"<<steps<<",\"time_ps\":"<<steps*dt/1e-12<<",\"max_torque_T\":"<<maxt<<",\"rms_torque_T\":"<<rms<<",\"energy_over_mu_per_spin\":"<<energy<<",\"preprojection_error\":"<<pre<<",\"norm_error\":"<<post<<",\"drift_since_sample\":"<<drift;
 if(driven)meta<<",\"field_y_T\":"<<by<<",\"sot_T\":"<<bt<<",\"active_updates\":"<<active_updates;
 meta<<"}\n";
 std::cout<<out<<" N="<<N<<" torque="<<maxt<<" norm="<<post<<"\n";
 }catch(std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
