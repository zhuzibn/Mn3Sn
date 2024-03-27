mmx_=zeros(natomW,natomL,gpurun_number*final_m_savestep);
mmy_=zeros(natomW,natomL,gpurun_number*final_m_savestep);
mmz_=zeros(natomW,natomL,gpurun_number*final_m_savestep);

BDSOT=zeros(natomW,natomL,'gpuArray');
BDSTT=zeros(natomW,natomL,'gpuArray');
muigpu=zeros(natomW,natomL,'gpuArray');
scalgpu=zeros(natomW,natomL,'gpuArray');
AsimnextL=zeros(natomW,natomL,'gpuArray');
AsimnextW=zeros(natomW,natomL,'gpuArray');
AsimpreviousL=zeros(natomW,natomL,'gpuArray');
AsimpreviousW=zeros(natomW,natomL,'gpuArray');

AsimpreviousWleft=zeros(natomW,natomL,'gpuArray');
AsimpreviousWright=zeros(natomW,natomL,'gpuArray');
AsimnextWleft=zeros(natomW,natomL,'gpuArray');
AsimnextWright=zeros(natomW,natomL,'gpuArray');
% if mod(ctW,2) %ctW奇数行
for ctL=1:natomL
    for ctW=1:natomW
        if atomtype_(ctW,ctL)==1%RE
            muigpu(ctW,ctL)=musRE;
            scalgpu(ctW,ctL)=gamRE/(1+alp^2);%scale parameter
            BDSOT(ctW,ctL)=BDSOTRE;
            BDSTT(ctW,ctL)=BDSTTRE;
        else

        end
AsimnextL(1:2:end,:)=zeros(L,w);  

        if atomtype_(ctW,ctL)==1%local atom is RE
            if ctL==natomL
                if bc
                    AsimnextL(ctW,ctL)=0;
                    AsimpreviousW(ctW,ctL)=0;
                end
            else
                if mod(ctW,4)==2 || mod(ctW,4)==0 %第二行和第四行
                    AsimnextL(ctW,ctL)=Jintra;
                else
                    AsimnextL(ctW,ctL)=0;
                end
            end

            if ctL==1
                if bc
                    AsimpreviousL(ctW,ctL)=0;
                    AsimpreviousWleft(ctW,ctL)=0;
                end
            else
                if mod(ctW,4)==2 || mod(ctW,4)==0 %第二行和第四行
                    AsimpreviousL(ctW,ctL)=Jintra;
                else
                    AsimpreviousL(ctW,ctL)=0;
                end
            end

            if ctW==natomW
                if bc
                    AsimpreviousW(ctW,ctL)=0;
                    AsimpreviousWleft(ctW,ctL)=0;
                    AsimpreviousWright(ctW,ctL)=0;
                end
            else
                if mod(ctW,4)==1
                    if mod(ctL,2)==1
                        if mod(ctL,2)==1
                            AsimpreviousWright(ctW,ctL)=Jintra;
                            AsimpreviousW(ctW,ctL)=Jintra;
                        end
                    end
                elseif mod(ctW,4)==2
                    if mod(ctL,2)==1
                        AsimpreviousWleft(ctW,ctL)=Jintra;
                    else
                        AsimpreviousW(ctW,ctL)=Jintra;
                    end
                elseif mod(ctW,4)==3
                    if mod(ctL,2)==0
                        AsimpreviousWright(ctW,ctL)=Jintra;
                        AsimpreviousW(ctW,ctL)=Jintra;
                    end
                elseif mod(ctW,4)==0
                    if mod(ctL,2)==1
                    AsimpreviousW(ctW,ctL)=Jintra;
                    else
                    AsimpreviousWleft(ctW,ctL)=Jintra;
                    end
                end
            end

            if ctW==1
                if bc
                    AsimnextW(ctW,ctL)=0;
                    AsimnextWleft(ctW,ctL)=0;
                    AsimnextWright(ctW,ctL)=0;
                end
            else
                if mod(ctW,4)==1
                    if mod(ctL,2)==1
                        AsimnextWright(ctW,ctL)=Jintra;
                        AsimnextW(ctW,ctL)=Jintra;
                    end
                elseif mod(ctW,4)==2                  
                    if mod(ctL,2)==1
                        AsimnextW(ctW,ctL)=Jintra;
                    else
                        AsimnextWleft(ctW,ctL)=Jintra;
                    end
                elseif mod(ctW,4)==3
                    if mod(ctL,2)==0
                        AsimnextWright(ctW,ctL)=Jintra;
                        AsimnextW(ctW,ctL)=Jintra;
                    end  
                elseif mod(ctW,4)==0
                    if mod(ctL,2)==1
                        AsimnextWleft(ctW,ctL)=Jintra;
                    else
                        AsimnextW(ctW,ctL)=Jintra;
                    end
                end
            end

        elseif atomtype_(ctW,ctL)==0%local atom is TM

        end
    end
end
gamatom=scalgpu*(1+alp^2);
clear ctW ctL
BFSOT=chi*BDSOT;
BFSTT=chi*BDSTT;

if dipolemode==3 %group atoms into macrocells
    if ~mod(natomW,natom_mc_W)
        nW_group=natomW/natom_mc_W;%number of macrocells along width direction
    else
        error('natomW should be multiple integer times of natom_mc_W')
    end

    if ~mod(natomL,natom_mc_L)
        nL_group=natomL/natom_mc_L;%number of macrocells along length direction
    else
        error('natomL should be multiple integer times of natom_mc_L')
    end

    n_macrocell=nW_group*nL_group;
    volume_mc=(natom_mc_W*d)*(natom_mc_L*d)*d;

    pmcW_=zeros(nW_group,nL_group);
    pmcL_=zeros(nW_group,nL_group);
    pmcz_=zeros(nW_group,nL_group);

    if(0)%use simpler numbers for debug
        d=1.1;
        muigpu=[1.2,1.4,1.2,1.2;...
            1.2,1.2,1.2,1.2;...
            1.2,1.2,1.2,1.4;...
            1.2,1.2,1.2,1.2];
        mmxtmp=[0.08,-0.08,0.08,0.08;...
            0.08,0.08,0.08,0.08;...
            0.08,0.08,0.08,-0.08;...
            0.08,0.08,0.08,0.08];
        mmytmp=zeros(4,4);
        mmztmp=zeros(4,4);
    end

    [rijL_tmp,rijW_tmp]=meshgrid(1:natomL,1:natomW);
    distL_mui=rijL_tmp.*muigpu;
    distW_mui=rijW_tmp.*muigpu;
    for ctL2=1:nL_group
        for ctW2=1:nW_group
            W_start=natom_mc_W*(ctW2-1)+1;
            W_end=natom_mc_W*ctW2;
            L_start=natom_mc_L*(ctL2-1)+1;
            L_end=natom_mc_L*ctL2;

            sum_distL_mui=sum(sum(distL_mui(W_start:W_end,L_start:L_end)));
            sum_distW_mui=sum(sum(distW_mui(W_start:W_end,L_start:L_end)));
            sum_mui=sum(sum(muigpu(W_start:W_end,L_start:L_end)));
            pmcL_(ctW2,ctL2)=d*sum_distL_mui/sum_mui;
            pmcW_(ctW2,ctL2)=d*sum_distW_mui/sum_mui;
        end
    end
    clear ctL2 ctW2
end

ct3run=round((runtime)/gpusave);
ct3=1;
while ~(ct3>ct3run)

    mmx=zeros(natomW,natomL,gpusteps,'gpuArray');
    mmy=zeros(natomW,natomL,gpusteps,'gpuArray');
    mmz=zeros(natomW,natomL,gpusteps,'gpuArray');

    if ~(ct3==1)
        mmx(:,:,1)=tmp2xn0;mmy(:,:,1)=tmp2yn0;mmz(:,:,1)=tmp2zn0;
    else
        mmx(:,:,1)=mx_init;mmy(:,:,1)=my_init;mmz(:,:,1)=mz_init;
    end
    clear tmpx tmpy tmpz
    ct1=1; %count 1
    while ct1<gpusteps

        mmxtmp=mmx(:,:,ct1);
        mmytmp=mmy(:,:,ct1);
        mmztmp=mmz(:,:,ct1);

        mmxnextL=circshift(mmxtmp,[0,-1]);%length direction
        mmynextL=circshift(mmytmp,[0,-1]);
        mmznextL=circshift(mmztmp,[0,-1]);
        mmxpreviousL=circshift(mmxtmp,[0,1]);
        mmypreviousL=circshift(mmytmp,[0,1]);
        mmzpreviousL=circshift(mmztmp,[0,1]);

        mmxnextW=circshift(mmxtmp,[1,0]);%width direction
        mmynextW=circshift(mmytmp,[1,0]);
        mmznextW=circshift(mmztmp,[1,0]);
        mmxpreviousW=circshift(mmxtmp,[-1,0]);
        mmypreviousW=circshift(mmytmp,[-1,0]);
        mmzpreviousW=circshift(mmztmp,[-1,0]);

        mmxnextWleft=circshift(mmxtmp,[1,1]);
        mmynextWleft=circshift(mmytmp,[1,1]);
        mmznextWleft=circshift(mmztmp,[1,1]);
        mmxnextWright=circshift(mmxtmp,[1,-1]);
        mmynextWright=circshift(mmytmp,[1,-1]);
        mmznextWright=circshift(mmztmp,[1,-1]);

        mmxpreviousWleft=circshift(mmxtmp,[-1,1]);
        mmypreviousWleft=circshift(mmytmp,[-1,1]);
        mmzpreviousWleft=circshift(mmztmp,[-1,1]);
        mmxpreviousWright=circshift(mmxtmp,[-1,-1]);
        mmypreviousWright=circshift(mmytmp,[-1,-1]);
        mmzpreviousWright=circshift(mmztmp,[-1,-1]);

        if bc%not periodic condition
            mmxnextL(:,end)=0;mmynextL(:,end)=0;mmznextL(:,end)=0;
            mmxpreviousL(:,1)=0;mmypreviousL(:,1)=0;mmzpreviousL(:,1)=0;
            mmxnextW(1,:)=0;mmynextW(1,:)=0;mmznextW(1,:)=0;
            mmxpreviousW(end,:)=0;mmypreviousW(end,:)=0;mmzpreviousW(end,:)=0;
            AsimnextW(:,end)=0;AsimpreviousWleft(:,1)=0;AsimpreviousW(:,end)=0;AsimnextWleft(:,1)=0;
        else%periodic condition
            %do nothing
        end

        hex_x=-(AsimnextL.*mmxnextL+AsimpreviousL.*mmxpreviousL+...
            AsimnextW.*mmxnextW+AsimpreviousW.*mmxpreviousW+...
            AsimpreviousWleft.*mmxpreviousWleft+AsimpreviousWright.*mmxpreviousWright+AsimnextWleft.*mmxnextWleft+AsimnextWright.*mmxnextWright)./muigpu;%[T]
        hex_y=-(AsimnextL.*mmynextL+AsimpreviousL.*mmypreviousL+...
            AsimnextW.*mmynextW+AsimpreviousW.*mmypreviousW+...
            AsimpreviousWleft.*mmypreviousWleft+AsimpreviousWright.*mmypreviousWright+AsimnextWleft.*mmynextWleft+AsimnextWright.*mmynextWright)./muigpu;
        hex_z=-(AsimnextL.*mmznextL+AsimpreviousL.*mmzpreviousL+...
            AsimnextW.*mmznextW+AsimpreviousW.*mmzpreviousW+...
            AsimpreviousWleft.*mmzpreviousWleft+AsimpreviousWright.*mmzpreviousWright+AsimnextWleft.*mmznextWleft+AsimnextWright.*mmznextWright)./muigpu;

        % hk_atom=Ksim./muigpu.*dot([coeffs_1,coeffs_2,coeffs_3],[mmxtmp,mmytmp,mmztmp]).*[coeffs_1,coeffs_2,coeffs_3];%hk_Fe in macrospin_model %anisotropy
        hani_x=Ksim./muigpu.*(coeffs_1.*mmxtmp+coeffs_2.*mmytmp+coeffs_3.*mmztmp).*coeffs_1;%hk_Fe in macrospin_model %anisotropy
        hani_y=Ksim./muigpu.*(coeffs_1.*mmxtmp+coeffs_2.*mmytmp+coeffs_3.*mmztmp).*coeffs_2;%hk_Fe in macrospin_model %anisotropy
        hani_z=Ksim./muigpu.*(coeffs_1.*mmxtmp+coeffs_2.*mmytmp+coeffs_3.*mmztmp).*coeffs_3;%hk_Fe in macrospin_model %anisotropy

        % Renjie pointed out the error in calculating hdmi in the previous
        % version of code, which has been corrected below.
        % 2023.09, During the plot of skyrmion, Zhang Xue find the top-down
        % direction is wrong. It is induced by the wrong nextW direction,
        % which has been corrected.
        hdmi_x=Dsim./muigpu.*(-mmznextL+mmzpreviousL);%[T]
        hdmi_y=Dsim./muigpu.*(-mmznextW+mmzpreviousW);
        hdmi_z=Dsim./muigpu.*(mmxnextL-mmxpreviousL+mmynextW-mmypreviousW);
        if thermalenable
            %equation (15) in Atomistic spin model simulations of magnetic nanomaterials
            %calculate once for one time step
            Hthermal1=sqrt(2*kb*T*alp./(muigpu.*gamatom.*tstep));%[T]
            Hthermalx=normrnd(0,Hthermal1);
            Hthermaly=normrnd(0,Hthermal1);
            Hthermalz=normrnd(0,Hthermal1);
        else
            Hthermalx=zeros(natomW,natomL,'gpuArray');
            Hthermaly=zeros(natomW,natomL,'gpuArray');
            Hthermalz=zeros(natomW,natomL,'gpuArray');
        end
        dipole_calc();
        hhx=hex_x+hani_x+hdmi_x+hdipolex_+Hext(1)+Hthermalx;
        hhy=hex_y+hani_y+hdmi_y+hdipoley_+Hext(2)+Hthermaly;
        hhz=hex_z+hani_z+hdmi_z+hdipolez_+Hext(3)+Hthermalz;
         if rk4==2%4th predictor-corrector
            if ct3==1 && ~(ct1>3)
                [sxx,syy,szz]=arrayfun(@atomgpurk4,mmxtmp,mmytmp,mmztmp,scalgpu,alp,...
                    tstep,hhx,hhy,hhz);
            else
                [sxx,syy,szz]=arrayfun(@atomgpupc4,tmpxn0,tmpyn0,tmpzn0,...
                    tmpxn1,tmpyn1,tmpzn1,tmpxn2,tmpyn2,tmpzn2,tmpxn3,tmpyn3,tmpzn3,...
                    scalgpu,alp,tstep,hhx,hhy,hhz);
            end
        elseif rk4==1 %rk4
            [sxx,syy,szz]=arrayfun(@atomgpurk4,mmxtmp,mmytmp,mmztmp,psjSHEx,...
                psjSHEy,psjSHEz,psjSTTx,psjSTTy,psjSTTz,scalgpu,alp,...
                tstep,hhx,hhy,hhz,BDSOT,BFSOT,BDSTT,BFSTT);
        else%heun
            [sxx,syy,szz]=arrayfun(@atomgpu,mmxtmp,mmytmp,mmztmp,scalgpu,alp,...
                tstep,hhx,hhy,hhz);%
        end

        mmx(:,:,ct1+1)=sxx; mmy(:,:,ct1+1)=syy; mmz(:,:,ct1+1)=szz;
        if enablefixedge
            mmx(:,1:fixededgeL,ct1+1)=mxleft;
            mmy(:,1:fixededgeL,ct1+1)=myleft;
            mmz(:,1:fixededgeL,ct1+1)=mzleft;

            mmx(:,natomL-fixededgeL:end,ct1+1)=mxright;
            mmy(:,natomL-fixededgeL:end,ct1+1)=myright;
            mmz(:,natomL-fixededgeL:end,ct1+1)=mzright;
        end
        ct1=ct1+1;
        if ~(ct3==1 && ~(ct1>3)) && ct1>3
            tmpxn0=mmx(:,:,ct1);tmpyn0=mmy(:,:,ct1);tmpzn0=mmz(:,:,ct1);
            tmpxn1=mmx(:,:,ct1-1);tmpyn1=mmy(:,:,ct1-1);tmpzn1=mmz(:,:,ct1-1);
            tmpxn2=mmx(:,:,ct1-2);tmpyn2=mmy(:,:,ct1-2);tmpzn2=mmz(:,:,ct1-2);
            tmpxn3=mmx(:,:,ct1-3);tmpyn3=mmy(:,:,ct1-3);tmpzn3=mmz(:,:,ct1-3);
        elseif ~(ct3==1 && ~(ct1>3)) && ct1==2
            tmpxn0=mmx(:,:,ct1);tmpyn0=mmy(:,:,ct1);tmpzn0=mmz(:,:,ct1);
            tmpxn1=tmp2xn0;tmpyn1=tmp2yn0;tmpzn1=tmp2zn0;
            tmpxn2=tmp2xn1;tmpyn2=tmp2yn1;tmpzn2=tmp2zn1;
            tmpxn3=tmp2xn2;tmpyn3=tmp2yn2;tmpzn3=tmp2zn2;
        elseif ~(ct3==1 && ~(ct1>3)) && ct1==3
            tmpxn0=mmx(:,:,ct1);tmpyn0=mmy(:,:,ct1);tmpzn0=mmz(:,:,ct1);
            tmpxn1=mmx(:,:,ct1-1);tmpyn1=mmy(:,:,ct1-1);tmpzn1=mmz(:,:,ct1-1);
            tmpxn2=tmp2xn0;tmpyn2=tmp2yn0;tmpzn2=tmp2zn0;
            tmpxn3=tmp2xn1;tmpyn3=tmp2yn1;tmpzn3=tmp2zn1;
        end
    end
    tmp2xn0=mmx(:,:,end);tmp2yn0=mmy(:,:,end);tmp2zn0=mmz(:,:,end);
    tmp2xn1=mmx(:,:,end-1);tmp2yn1=mmy(:,:,end-1);tmp2zn1=mmz(:,:,end-1);
    tmp2xn2=mmx(:,:,end-2);tmp2yn2=mmy(:,:,end-2);tmp2zn2=mmz(:,:,end-2);
    mmx_(:,:,(ct3-1)*final_m_savestep+1:ct3*final_m_savestep)=gather(mmx(:,:,1:savetstep:end));
    mmy_(:,:,(ct3-1)*final_m_savestep+1:ct3*final_m_savestep)=gather(mmy(:,:,1:savetstep:end));
    mmz_(:,:,(ct3-1)*final_m_savestep+1:ct3*final_m_savestep)=gather(mmz(:,:,1:savetstep:end));
    ct3=ct3+1;
end

clear mmx mmy mmz tmp2xn0 tmp2yn0 tmp2zn0 tmp2xn1 tmp2yn1 tmp2zn1
clear tmp2xn2 tmp2yn2 tmp2zn2
mmx=mmx_;
mmy=mmy_;
mmz=mmz_;
clear mmx_ mmy_ mmz_
t=t(1:savetstep:end);
