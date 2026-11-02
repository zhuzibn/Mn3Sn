function [dx,dy,dz]=atom_llgs_rhs(ssx,ssy,ssz,psjSHEx,psjSHEy,psjSHEz,psjSTTx,psjSTTy,psjSTTz,scal,alph,hhx,hhy,hhz,BdSOT,BfSOT,BdSTT,BfSTT)
% Explicit LLGS derivative in 1/s; fields and torque amplitudes are in T.
% Original atomgpurk4 first-stage algebra, including Gilbert cross terms.
% No timestep or normalization here: the lattice integrator owns RK stages.
u1=hhx;u2=hhy;u3=hhz;
v1=ssx;v2=ssy;v3=ssz;
dsdt1x=u2*v3-u3*v2;
dsdt1y=u3*v1-u1*v3;
dsdt1z=u1*v2-u2*v1;
%cross(cross(sss,hh),sss) Beff DLT
u1=-dsdt1x;u2=-dsdt1y;u3=-dsdt1z;
v1=ssx;v2=ssy;v3=ssz;
dsdt2x=u2*v3-u3*v2;
dsdt2y=u3*v1-u1*v3;
dsdt2z=u1*v2-u2*v1;
%cross(-sss,cross(sss,ey)) BdSOT DLT
u1tmp=ssx;u2tmp=ssy;u3tmp=ssz;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=-ssx;u2=-ssy;u3=-ssz;
dbdSOTx=u2*v3-u3*v2;
dbdSOTy=u3*v1-u1*v3;
dbdSOTz=u1*v2-u2*v1;
%cross(sss,ey) BdSOT FLT
u1=ssx;u2=ssy;u3=ssz;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbdSOTx=u2*v3-u3*v2;
fbdSOTy=u3*v1-u1*v3;
fbdSOTz=u1*v2-u2*v1;
%cross(sss,cross(sss,ey)) BfSOT DLT
u1tmp=ssx;u2tmp=ssy;u3tmp=ssz;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=ssx;u2=ssy;u3=ssz;
dbfSOTx=u2*v3-u3*v2;
dbfSOTy=u3*v1-u1*v3;
dbfSOTz=u1*v2-u2*v1;
%cross(sss,ey) BfSOT FLT
u1=ssx;u2=ssy;u3=ssz;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbfSOTx=u2*v3-u3*v2;
fbfSOTy=u3*v1-u1*v3;
fbfSOTz=u1*v2-u2*v1;
%cross(-sss,cross(sss,ey)) BdSTT DLT
u1tmp=ssx;u2tmp=ssy;u3tmp=ssz;
v1tmp=psjSTTx;v2tmp=psjSTTy;v3tmp=psjSTTz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=-ssx;u2=-ssy;u3=-ssz;
dbdSTTx=u2*v3-u3*v2;
dbdSTTy=u3*v1-u1*v3;
dbdSTTz=u1*v2-u2*v1;
%cross(sss,ey) BdSTT FLT
u1=ssx;u2=ssy;u3=ssz;
v1=psjSTTx;v2=psjSTTy;v3=psjSTTz;
fbdSTTx=u2*v3-u3*v2;
fbdSTTy=u3*v1-u1*v3;
fbdSTTz=u1*v2-u2*v1;
%cross(sss,cross(sss,ey)) BfSTT DLT
u1tmp=ssx;u2tmp=ssy;u3tmp=ssz;
v1tmp=psjSTTx;v2tmp=psjSTTy;v3tmp=psjSTTz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=ssx;u2=ssy;u3=ssz;
dbfSTTx=u2*v3-u3*v2;
dbfSTTy=u3*v1-u1*v3;
dbfSTTz=u1*v2-u2*v1;
%cross(sss,ey) BfSTT FLT
u1=ssx;u2=ssy;u3=ssz;
v1=psjSTTx;v2=psjSTTy;v3=psjSTTz;
fbfSTTx=u2*v3-u3*v2;
fbfSTTy=u3*v1-u1*v3;
fbfSTTz=u1*v2-u2*v1;
%
dsdtx=dsdt1x+alph*dsdt2x+BdSOT*dbdSOTx+alph*BdSOT*fbdSOTx+alph*BfSOT*dbfSOTx+BfSOT*fbfSOTx+...
    BdSTT*dbdSTTx+alph*BdSTT*fbdSTTx+alph*BfSTT*dbfSTTx+BfSTT*fbfSTTx;
dsdty=dsdt1y+alph*dsdt2y+BdSOT*dbdSOTy+alph*BdSOT*fbdSOTy+alph*BfSOT*dbfSOTy+BfSOT*fbfSOTy+...
    BdSTT*dbdSTTy+alph*BdSTT*fbdSTTy+alph*BfSTT*dbfSTTy+BfSTT*fbfSTTy;
dsdtz=dsdt1z+alph*dsdt2z+BdSOT*dbdSOTz+alph*BdSOT*fbdSOTz+alph*BfSOT*dbfSOTz+BfSOT*fbfSOTz+...
    BdSTT*dbdSTTz+alph*BdSTT*fbdSTTz+alph*BfSTT*dbfSTTz+BfSTT*fbfSTTz;
%
dx=scal*dsdtx; dy=scal*dsdty; dz=scal*dsdtz;
end
