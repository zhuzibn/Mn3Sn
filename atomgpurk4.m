%LLG solver for gpu calculation using Heun Method
% "clear xx" is not allowed in gpu version of arrayfun
function [sxx,syy,szz]=atomgpurk4(ssx,ssy,ssz,psjSHEx,psjSHEy,psjSHEz,psjSTTx,psjSTTy,psjSTTz,scal,alph,ts,hhx,hhy,hhz,BdSOT,BfSOT,BdSTT,BfSTT)
%cross(u,v)=(u2v3-u3v2)i+(u3v1-u1v3)j+(u1v2-u2v1)k
%------------------kk1--------------------------
%-cross(sss,hh)=cross(hh,sss) Beff FLT
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
kk1x=ts*scal*dsdtx;kk1y=ts*scal*dsdty;kk1z=ts*scal*dsdtz;
%-------------------kk2------------------------
%sss=ss1+ts/2*kk1;%y[i+1]
sxtmp=ssx+1/2*kk1x;
sytmp=ssy+1/2*kk1y;
sztmp=ssz+1/2*kk1z;
%cross(u,v)=(u2v3-u3v2)i+(u3v1-u1v3)j+(u1v2-u2v1)k
%-cross(sss,hh)=cross(hh,sss) Beff FLT
u1=hhx;u2=hhy;u3=hhz;
v1=sxtmp;v2=sytmp;v3=sztmp;
dsdt1x=u2*v3-u3*v2;
dsdt1y=u3*v1-u1*v3;
dsdt1z=u1*v2-u2*v1;
%cross(cross(sss,hh),sss) Beff DLT
u1=-dsdt1x;u2=-dsdt1y;u3=-dsdt1z;
v1=sxtmp;v2=sytmp;v3=sztmp;
dsdt2x=u2*v3-u3*v2;
dsdt2y=u3*v1-u1*v3;
dsdt2z=u1*v2-u2*v1;
%cross(-sss,cross(sss,ey)) BdSOT DLT
u1tmp=sxtmp;u2tmp=sytmp;u3tmp=sztmp;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=-sxtmp;u2=-sytmp;u3=-sztmp;
dbdSOTx=u2*v3-u3*v2;
dbdSOTy=u3*v1-u1*v3;
dbdSOTz=u1*v2-u2*v1;
%cross(sss,ey) BdSOT FLT
u1=sxtmp;u2=sytmp;u3=sztmp;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbdSOTx=u2*v3-u3*v2;
fbdSOTy=u3*v1-u1*v3;
fbdSOTz=u1*v2-u2*v1;
%cross(sss,cross(sss,ey)) BfSOT DLT
u1tmp=sxtmp;u2tmp=sytmp;u3tmp=sztmp;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=sxtmp;u2=sytmp;u3=sztmp;
dbfSOTx=u2*v3-u3*v2;
dbfSOTy=u3*v1-u1*v3;
dbfSOTz=u1*v2-u2*v1;
%cross(sss,ey) BfSOT FLT
u1=sxtmp;u2=sytmp;u3=sztmp;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbfSOTx=u2*v3-u3*v2;
fbfSOTy=u3*v1-u1*v3;
fbfSOTz=u1*v2-u2*v1;
%cross(-sss,cross(sss,ey)) BdSTT DLT
u1tmp=sxtmp;u2tmp=sytmp;u3tmp=sztmp;
v1tmp=psjSTTx;v2tmp=psjSTTy;v3tmp=psjSTTz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=-sxtmp;u2=-sytmp;u3=-sztmp;
dbdSTTx=u2*v3-u3*v2;
dbdSTTy=u3*v1-u1*v3;
dbdSTTz=u1*v2-u2*v1;
%cross(sss,ey) BdSTT FLT
u1=sxtmp;u2=sytmp;u3=sztmp;
v1=psjSTTx;v2=psjSTTy;v3=psjSTTz;
fbdSTTx=u2*v3-u3*v2;
fbdSTTy=u3*v1-u1*v3;
fbdSTTz=u1*v2-u2*v1;
%cross(sss,cross(sss,ey)) BfSTT DLT
u1tmp=sxtmp;u2tmp=sytmp;u3tmp=sztmp;
v1tmp=psjSTTx;v2tmp=psjSTTy;v3tmp=psjSTTz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=sxtmp;u2=sytmp;u3=sztmp;
dbfSTTx=u2*v3-u3*v2;
dbfSTTy=u3*v1-u1*v3;
dbfSTTz=u1*v2-u2*v1;
%cross(sss,ey) BfSTT FLT
u1=sxtmp;u2=sytmp;u3=sztmp;
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
kk2x=ts*scal*dsdtx;kk2y=ts*scal*dsdty;kk2z=ts*scal*dsdtz;
%-------------------kk3------------------------
%sss=ss1+ts/2*kk1;%y[i+1]
sxtmp2=ssx+1/2*kk2x;
sytmp2=ssy+1/2*kk2y;
sztmp2=ssz+1/2*kk2z;

%cross(u,v)=(u2v3-u3v2)i+(u3v1-u1v3)j+(u1v2-u2v1)k
%-cross(sss,hh)=cross(hh,sss) Beff FLT
u1=hhx;u2=hhy;u3=hhz;
v1=sxtmp2;v2=sytmp2;v3=sztmp2;
dsdt1x=u2*v3-u3*v2;
dsdt1y=u3*v1-u1*v3;
dsdt1z=u1*v2-u2*v1;
%cross(cross(sss,hh),sss) Beff DLT
u1=-dsdt1x;u2=-dsdt1y;u3=-dsdt1z;
v1=sxtmp2;v2=sytmp2;v3=sztmp2;
dsdt2x=u2*v3-u3*v2;
dsdt2y=u3*v1-u1*v3;
dsdt2z=u1*v2-u2*v1;
%cross(-sss,cross(sss,ey)) BdSOT DLT
u1tmp=sxtmp2;u2tmp=sytmp2;u3tmp=sztmp2;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=-sxtmp2;u2=-sytmp2;u3=-sztmp2;
dbdSOTx=u2*v3-u3*v2;
dbdSOTy=u3*v1-u1*v3;
dbdSOTz=u1*v2-u2*v1;
%cross(sss,ey) BdSOT FLT
u1=sxtmp2;u2=sytmp2;u3=sztmp2;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbdSOTx=u2*v3-u3*v2;
fbdSOTy=u3*v1-u1*v3;
fbdSOTz=u1*v2-u2*v1;
%cross(sss,cross(sss,ey)) BfSOT DLT
u1tmp=sxtmp2;u2tmp=sytmp2;u3tmp=sztmp2;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=sxtmp2;u2=sytmp2;u3=sztmp2;
dbfSOTx=u2*v3-u3*v2;
dbfSOTy=u3*v1-u1*v3;
dbfSOTz=u1*v2-u2*v1;
%cross(sss,ey) BfSOT FLT
u1=sxtmp2;u2=sytmp2;u3=sztmp2;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbfSOTx=u2*v3-u3*v2;
fbfSOTy=u3*v1-u1*v3;
fbfSOTz=u1*v2-u2*v1;
%cross(-sss,cross(sss,ey)) BdSTT DLT
u1tmp=sxtmp2;u2tmp=sytmp2;u3tmp=sztmp2;
v1tmp=psjSTTx;v2tmp=psjSTTy;v3tmp=psjSTTz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=-sxtmp2;u2=-sytmp2;u3=-sztmp2;
dbdSTTx=u2*v3-u3*v2;
dbdSTTy=u3*v1-u1*v3;
dbdSTTz=u1*v2-u2*v1;
%cross(sss,ey) BdSTT FLT
u1=sxtmp2;u2=sytmp2;u3=sztmp2;
v1=psjSTTx;v2=psjSTTy;v3=psjSTTz;
fbdSTTx=u2*v3-u3*v2;
fbdSTTy=u3*v1-u1*v3;
fbdSTTz=u1*v2-u2*v1;
%cross(sss,cross(sss,ey)) BfSTT DLT
u1tmp=sxtmp2;u2tmp=sytmp2;u3tmp=sztmp2;
v1tmp=psjSTTx;v2tmp=psjSTTy;v3tmp=psjSTTz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=sxtmp2;u2=sytmp2;u3=sztmp2;
dbfSTTx=u2*v3-u3*v2;
dbfSTTy=u3*v1-u1*v3;
dbfSTTz=u1*v2-u2*v1;
%cross(sss,ey) BfSTT FLT
u1=sxtmp2;u2=sytmp2;u3=sztmp2;
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
kk3x=ts*scal*dsdtx;kk3y=ts*scal*dsdty;kk3z=ts*scal*dsdtz;
%-------------------kk4------------------------
%sss=ss1+ts/2*kk1;%y[i+1]
sxtmp3=ssx+kk3x;
sytmp3=ssy+kk3y;
sztmp3=ssz+kk3z;

%cross(u,v)=(u2v3-u3v2)i+(u3v1-u1v3)j+(u1v2-u2v1)k
%-cross(sss,hh)=cross(hh,sss) Beff FLT
u1=hhx;u2=hhy;u3=hhz;
v1=sxtmp3;v2=sytmp3;v3=sztmp3;
dsdt1x=u2*v3-u3*v2;
dsdt1y=u3*v1-u1*v3;
dsdt1z=u1*v2-u2*v1;
%cross(cross(sss,hh),sss) Beff DLT
u1=-dsdt1x;u2=-dsdt1y;u3=-dsdt1z;
v1=sxtmp3;v2=sytmp3;v3=sztmp3;
dsdt2x=u2*v3-u3*v2;
dsdt2y=u3*v1-u1*v3;
dsdt2z=u1*v2-u2*v1;
%cross(-sss,cross(sss,ey)) BdSOT DLT
u1tmp=sxtmp3;u2tmp=sytmp3;u3tmp=sztmp3;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=-sxtmp3;u2=-sytmp3;u3=-sztmp3;
dbdSOTx=u2*v3-u3*v2;
dbdSOTy=u3*v1-u1*v3;
dbdSOTz=u1*v2-u2*v1;
%cross(sss,ey) BdSOT FLT
u1=sxtmp3;u2=sytmp3;u3=sztmp3;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbdSOTx=u2*v3-u3*v2;
fbdSOTy=u3*v1-u1*v3;
fbdSOTz=u1*v2-u2*v1;
%cross(sss,cross(sss,ey)) BfSOT DLT
u1tmp=sxtmp3;u2tmp=sytmp3;u3tmp=sztmp3;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=sxtmp3;u2=sytmp3;u3=sztmp3;
dbfSOTx=u2*v3-u3*v2;
dbfSOTy=u3*v1-u1*v3;
dbfSOTz=u1*v2-u2*v1;
%cross(sss,ey) BfSOT FLT
u1=sxtmp3;u2=sytmp3;u3=sztmp3;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbfSOTx=u2*v3-u3*v2;
fbfSOTy=u3*v1-u1*v3;
fbfSOTz=u1*v2-u2*v1;
%cross(-sss,cross(sss,ey)) BdSTT DLT
u1tmp=sxtmp3;u2tmp=sytmp3;u3tmp=sztmp3;
v1tmp=psjSTTx;v2tmp=psjSTTy;v3tmp=psjSTTz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=-sxtmp3;u2=-sytmp3;u3=-sztmp3;
dbdSTTx=u2*v3-u3*v2;
dbdSTTy=u3*v1-u1*v3;
dbdSTTz=u1*v2-u2*v1;
%cross(sss,ey) BdSTT FLT
u1=sxtmp3;u2=sytmp3;u3=sztmp3;
v1=psjSTTx;v2=psjSTTy;v3=psjSTTz;
fbdSTTx=u2*v3-u3*v2;
fbdSTTy=u3*v1-u1*v3;
fbdSTTz=u1*v2-u2*v1;
%cross(sss,cross(sss,ey)) BfSTT DLT
u1tmp=sxtmp3;u2tmp=sytmp3;u3tmp=sztmp3;
v1tmp=psjSTTx;v2tmp=psjSTTy;v3tmp=psjSTTz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=sxtmp3;u2=sytmp3;u3=sztmp3;
dbfSTTx=u2*v3-u3*v2;
dbfSTTy=u3*v1-u1*v3;
dbfSTTz=u1*v2-u2*v1;
%cross(sss,ey) BfSTT FLT
u1=sxtmp3;u2=sytmp3;u3=sztmp3;
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
kk4x=ts*scal*dsdtx;kk4y=ts*scal*dsdty;kk4z=ts*scal*dsdtz;
%-----------------final---------------------
snx=ssx+1/6*(kk1x+2*kk2x+2*kk3x+kk4x);
sny=ssy+1/6*(kk1y+2*kk2y+2*kk3y+kk4y);
snz=ssz+1/6*(kk1z+2*kk2z+2*kk3z+kk4z);
normsn=sqrt(snx^2+sny^2+snz^2);
snx=snx/normsn;
sny=sny/normsn;
snz=snz/normsn;

sxx=snx;syy=sny;szz=snz;
end