distrib();
mx_init=zeros(natomW,natomL);%initial magnetization
my_init=zeros(natomW,natomL);
mz_init=zeros(natomW,natomL);

if dwcalc
    phi_=0;
    for ctL=1:natomL
        for ctW=1:natomW
            if ctL<round(natomL/2)
                if atomtype_(ctW,ctL)==0%TM
                    thet_=5/180*pi;
                else
                    thet_=(5+180)/180*pi;
                end
                mx_init(ctW,ctL)=sin(thet_)*cos(phi_);
                my_init(ctW,ctL)=sin(thet_)*sin(phi_);
                mz_init(ctW,ctL)=cos(thet_);
            else
                if atomtype_(ctW,ctL)==0%TM
                    thet_=(5+180)/180*pi;
                else
                    thet_=5/180*pi;
                end
                mx_init(ctW,ctL)=sin(thet_)*cos(phi_);
                my_init(ctW,ctL)=sin(thet_)*sin(phi_);
                mz_init(ctW,ctL)=cos(thet_);
            end
        end
    end
else
    for ctL=1:natomL
        for ctW=1:natomW
            phi_=90/180*pi;% keep in y-z plane
            thet_1=30/180*pi;
            if atomtype_(ctW,ctL)==1
                if mod(ctW,2)==1 %奇数行
                    thet_=(thet_1+120/180*pi);
                elseif mod(ctW,4)==2 && mod(ctL,2)==1 %2行奇数列
                    thet_=(thet_1+240/180*pi);
                elseif mod(ctW,4)==0 && mod(ctL,2)==0 %4行偶数列
                    thet_=(thet_1+240/180*pi);
                elseif mod(ctW,4)==2 && mod(ctL,2)==0 %2行偶数列
                    thet_=(thet_1);
                elseif mod(ctW,4)==0 && mod(ctL,2)==1 %4行奇数列
                    thet_=(thet_1);
                end
            else
                phi_=0/180*pi;
                thet_=90/180*pi;
            end
            mx_init(ctW,ctL)=sin(thet_)*cos(phi_);
            my_init(ctW,ctL)=sin(thet_)*sin(phi_);
            mz_init(ctW,ctL)=cos(thet_);
        end
    end
end
clear ctL ctW

if loadstartm
    clear mx_init my_init mz_init
    load(startmname);
    if natomW~=natomWcheck || natomL~=natomLcheck
        error('system not consistent')
    end
    clear natomxcheck natomycheck
    mx_init=mmxstart;
    my_init=mmystart;
    mz_init=mmzstart;
end
