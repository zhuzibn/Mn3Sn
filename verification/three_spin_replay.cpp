// Diagnostic only: fixed deterministic triangle, independently checked in MATLAB.
// args: output params mode dt_fs seed_deg field_Oe current duration_ns on_ns width_ns sample_ps [seed.csv]
#include <array>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
using Vec=std::array<double,3>;
using State=std::array<Vec,3>;
struct Parameters { double j,d,k,gamma,alpha,b; };
static Vec cross(const Vec&a,const Vec&b) {
    return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
}
static double dot(const Vec&a,const Vec&b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
static State axes;
static State fields(const State&m,const Parameters&p,double B,bool reference) {
    State h{};
    for(int i=0;i<3;i++) {
        int j=(i+1)%3,k=(i+2)%3;
        double projection=dot(m[i],axes[i]);
        if(reference) {
            // Ascending neighbor order and separate multiplications in field_eta.m.
            int a=i==0?1:0, b=i==2?1:2;
            for(int q=0;q<3;q++) h[i][q]=p.k*projection*axes[i][q]+(q==1?B:0)
                +(-p.j*m[a][q]-p.j*m[b][q]);
        } else {
            for(int q=0;q<3;q++) h[i][q]=-p.j*(m[j][q]+m[k][q])
                +p.k*projection*axes[i][q]+(q==1?B:0);
        }
        h[i][1]+=p.d*(-m[j][2]+m[k][2]);
        h[i][2]+=p.d*(m[j][1]-m[k][1]);
    }
    return h;
}
static State rhs(const State&m,const State&h,const Parameters&p,double b) {
    State r{}; const Vec polarization{1,0,0};
    double scale=p.gamma/(1+p.alpha*p.alpha);
    for(int i=0;i<3;i++) {
        Vec mh=cross(m[i],h[i]),mmh=cross(m[i],mh),mp=cross(m[i],polarization),mmp=cross(m[i],mp);
        for(int q=0;q<3;q++) r[i][q]=scale*(-mh[q]-p.alpha*mmh[q]-b*mmp[q]+p.alpha*b*mp[q]);
    }
    return r;
}
static State shifted(const State&m,const State&k,double dt) {
    State u{};for(int i=0;i<3;i++) for(int q=0;q<3;q++) u[i][q]=m[i][q]+dt*k[i][q];return u;
}
static void step(State&m,const Parameters&p,double B,double b,double dt,int mode,double&pre,double&post) {
    bool frozen=mode<2,reference=mode==3;
    State h=fields(m,p,B,reference),k1=rhs(m,h,p,b);
    State u=shifted(m,k1,dt/2);
    State k2=rhs(u,frozen?h:fields(u,p,B,reference),p,b);
    u=shifted(m,k2,dt/2);
    State k3=rhs(u,frozen?h:fields(u,p,B,reference),p,b);
    u=shifted(m,k3,dt);
    State k4=rhs(u,frozen?h:fields(u,p,B,reference),p,b);
    for(int i=0;i<3;i++) {
        for(int q=0;q<3;q++) m[i][q]+=dt/6*(k1[i][q]+2*k2[i][q]+2*k3[i][q]+k4[i][q]);
        double norm=std::sqrt(dot(m[i],m[i]));
        if(!std::isfinite(norm)||norm==0) throw std::runtime_error("nonfinite/zero occupied spin");
        pre=std::max(pre,std::abs(norm-1));
        for(double &v:m[i]) v/=norm;
        post=std::max(post,std::abs(std::sqrt(dot(m[i],m[i]))-1));
    }
}
static void output(std::ofstream&f,const State&m,long completed,double dt) {
    // Little-endian native doubles: physical seconds, then ax ay az bx ... cz.
    double t=completed*dt;f.write(reinterpret_cast<char*>(&t),sizeof(t));
    for(auto &spin:m) f.write(reinterpret_cast<const char*>(spin.data()),3*sizeof(double));
}
int main(int argc,char**argv) {
    if(argc!=13&&argc!=14) return 2;
    Parameters p;std::ifstream input(argv[2]);input>>p.j>>p.d>>p.k>>p.gamma>>p.alpha>>p.b;
    if(!input) return 2;
    int mode=std::stoi(argv[3]);double dt=std::stod(argv[4])*1e-15,seed=std::stod(argv[5]);
    double B=std::stod(argv[6])*1e-4,current=std::stod(argv[7]);
    long count=std::lround(std::stod(argv[8])*1e-9/dt),on=std::lround(std::stod(argv[9])*1e-9/dt);
    long width=std::lround(std::stod(argv[10])*1e-9/dt),sample=std::lround(std::stod(argv[11])*1e-12/dt);
    // Last argument is an explicit output contract tag to prevent accidental legacy invocation.
    if(std::string(argv[12])!="physical"||sample<1||count<1||on<0||on+width>count||mode<0||mode>3) return 2;
    double pi=std::acos(-1.0),aa[3]={pi/6,9*pi/6,5*pi/6};State m{};
    for(int i=0;i<3;i++) {
        axes[i]={0,std::sin(aa[i]),std::cos(aa[i])};
        double a=seed/180*pi+i*120.0/180*pi;
        m[i]={std::sin(a)*std::cos(90.0/180*pi),std::sin(a)*std::sin(90.0/180*pi),std::cos(a)};
    }
    if(argc==14) {
        std::ifstream s(argv[13]);char comma;
        for(auto&v:m) {s>>v[0]>>comma>>v[1]>>comma>>v[2];}
        if(!s) return 2;
    }
    std::ofstream out(argv[1],std::ios::binary);if(!out) return 3;
    double pre=0,post=0;output(out,m,0,dt);
    for(long n=0;n<count;n++) {
        bool driven=n>=on&&n<on+width;
        step(m,p,driven?B:0,driven?p.b*current/5e10:0,dt,mode,pre,post);
        if((n+1)%sample==0||n+1==count) output(out,m,n+1,dt);
    }
    Vec mean{};for(int q=0;q<3;q++) for(auto&v:m) mean[q]+=v[q]/3;
    std::cout<<std::setprecision(17)<<"{\"updates\":"<<count<<",\"final_time_s\":"<<count*dt
        <<",\"max_preprojection_norm_error\":"<<pre<<",\"max_postprojection_norm_error\":"<<post
        <<",\"final_angle_deg\":"<<std::atan2(mean[1],mean[2])*180/pi<<",\"mean\":["
        <<mean[0]<<","<<mean[1]<<","<<mean[2]<<"]}\n";
}
