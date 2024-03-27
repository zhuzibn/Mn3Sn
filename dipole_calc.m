switch dipolemode
    case 0
        hdipolex_=zeros(natomW,natomL);
        hdipoley_=zeros(natomW,natomL);
        hdipolez_=zeros(natomW,natomL);
    case 1
        %cpu calculation of dipole field, used for benchmaking
        cpuhdipolex=zeros(natomW,natomL);
        cpuhdipoley=zeros(natomW,natomL);
        cpuhdipolez=zeros(natomW,natomL);
        mmxttmp=gather(mmxtmp);
        mmyttmp=gather(mmytmp);
        mmzttmp=gather(mmztmp);
        for ctL=1:natomL
            for ctW=1:natomW
                for ctL2=1:natomL
                    for ctW2=1:natomW
                        if ctW==ctW2 && ctL==ctL2
                            tmp=[0,0,0];
                        else
                            dist=sqrt(((ctW2-ctW)*d)^2+((ctL2-ctL)*d)^2);
                            rij=[(ctL2-ctL)*d,(ctW2-ctW)*d,0];
                            sj=[mmxttmp(ctW2,ctL2),mmyttmp(ctW2,ctL2),mmzttmp(ctW2,ctL2)];
                            tmp=mu_0/(4*pi)*(3*rij*gather(muigpu(ctW2,ctL2))*dot(sj,rij)/dist^5-...
                                sj*gather(muigpu(ctW2,ctL2))/dist^3);
                        end
                        cpuhdipolex(ctW,ctL)=cpuhdipolex(ctW,ctL)+tmp(1);
                        cpuhdipoley(ctW,ctL)=cpuhdipoley(ctW,ctL)+tmp(2);
                        cpuhdipolez(ctW,ctL)=cpuhdipolez(ctW,ctL)+tmp(3);
                    end
                end
            end
        end
        clear mmxttmp mmyttmp mmzttmp ctW ctL ctW2 ctL2 sj tmp
    case 2%the code might be wrong, don't use for now
        %gpu calculation of dipole field, used for production
        hdipolex_=zeros(natomW,natomL,'gpuArray');
        hdipoley_=zeros(natomW,natomL,'gpuArray');
        hdipolez_=zeros(natomW,natomL,'gpuArray');
        for ctL=1:natomL
            for ctW=1:natomW
                [rijL_tmp,rijW_tmp]=meshgrid(0:natomL-1,0:natomW-1);
                rijL_tmp=rijL_tmp-ctL+1;%note:rijW_tmp corresponds to column
                rijW_tmp=rijW_tmp-ctW+1;%note:rijL_tmp corresponds to row
                dist_=d*sqrt(rijL_tmp.^2+rijW_tmp.^2);
                rijL_=d*rijL_tmp;
                rijW_=d*rijW_tmp;
                dot_mr=muigpu.*(mmxtmp.*rijW_+mmytmp.*rijL_);%problematic, should exchange mmxtmp and mmytmp?
                hdipolex=mu_0/(4*pi)*(3*rijW_.*dot_mr./dist_.^5-muigpu.*mmxtmp./dist_.^3);%[T]
                %problematic, rijW_ should change to rijL_?
                hdipolex(ctW,ctL)=0;
                hdipoley=mu_0/(4*pi)*(3*rijL_.*dot_mr./dist_.^5-muigpu.*mmytmp./dist_.^3);
                hdipoley(ctW,ctL)=0;
                hdipolez=mu_0/(4*pi)*(-muigpu.*mmztmp./dist_.^3);
                hdipolez(ctW,ctL)=0;
                hdipolex_(ctW,ctL)=sum(sum(hdipolex));
                hdipoley_(ctW,ctL)=sum(sum(hdipoley));
                hdipolez_(ctW,ctL)=sum(sum(hdipolez));
            end
        end
        clear ctW ctL
    case 3
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
        % cpu calculation with macrocell method, refer to 2014-Atomistic spin model simulations of-JPCM-R F L Evans-print-1-14
        mui_sx=muigpu.*mmxtmp;
        mui_sy=muigpu.*mmytmp;
        mui_sz=muigpu.*mmytmp;
        
        mmcx_=zeros(nW_group,nL_group);
        mmcy_=zeros(nW_group,nL_group);
        mmcz_=zeros(nW_group,nL_group);
        
        for ctL2=1:nL_group
            for ctW2=1:nW_group
                W_start=natom_mc_W*(ctW2-1)+1;
                W_end=natom_mc_W*ctW2;
                L_start=natom_mc_L*(ctL2-1)+1;
                L_end=natom_mc_L*ctL2;
                
                mmcx_(ctW2,ctL2)=sum(sum(mui_sx(W_start:W_end,L_start:L_end)));
                mmcy_(ctW2,ctL2)=sum(sum(mui_sy(W_start:W_end,L_start:L_end)));
                mmcz_(ctW2,ctL2)=sum(sum(mui_sz(W_start:W_end,L_start:L_end)));
            end
        end
        clear ctW2 ctL2
        
        hdipolex_mc=zeros(nW_group,nL_group);
        hdipoley_mc=zeros(nW_group,nL_group);
        hdipolez_mc=zeros(nW_group,nL_group);
        
        for ctL2=1:nL_group
            for ctW2=1:nW_group
                rijL_=pmcL_-pmcL_(ctW2,ctL2);
                rijW_=pmcW_-pmcW_(ctW2,ctL2);
                rijz_=pmcz_-pmcz_(ctW2,ctL2);
                dist_=sqrt(rijL_.^2+rijW_.^2+rijz_.^2);
                dot_mr=(mmcx_.*rijL_+mmcy_.*rijW_+mmcz_.*rijz_);
                hdipolex=mu_0/(4*pi)*(3*rijL_.*dot_mr./dist_.^5-mmcx_./dist_.^3);%[T]
                %The first component should be divided by dist_.^5 since r is not normalized
                hdipolex(ctW2,ctL2)=0;
                hdipoley=mu_0/(4*pi)*(3*rijW_.*dot_mr./dist_.^5-mmcy_./dist_.^3);
                hdipoley(ctW2,ctL2)=0;
                hdipolez=mu_0/(4*pi)*(3*rijz_.*dot_mr./dist_.^5-mmcz_./dist_.^3);
                hdipolez(ctW2,ctL2)=0;
                hdipolex_mc(ctW2,ctL2)=sum(sum(hdipolex))-mu_0/3*mmcx_(ctW2,ctL2)/volume_mc;
                hdipoley_mc(ctW2,ctL2)=sum(sum(hdipoley))-mu_0/3*mmcy_(ctW2,ctL2)/volume_mc;
                hdipolez_mc(ctW2,ctL2)=sum(sum(hdipolez))-mu_0/3*mmcz_(ctW2,ctL2)/volume_mc;
            end
        end
        clear ctW2 ctL2
        %distribute to all atoms
        hdipolex_=zeros(natomW,natomL);
        hdipoley_=zeros(natomW,natomL);
        hdipolez_=zeros(natomW,natomL);
        for ctL2=1:nL_group
            for ctW2=1:nW_group
                W_start=natom_mc_W*(ctW2-1)+1;
                W_end=natom_mc_W*ctW2;
                L_start=natom_mc_L*(ctL2-1)+1;
                L_end=natom_mc_L*ctL2;
                
                hdipolex_(W_start:W_end,L_start:L_end)=hdipolex_mc(ctW2,ctL2);
                hdipoley_(W_start:W_end,L_start:L_end)=hdipoley_mc(ctW2,ctL2);
                hdipolez_(W_start:W_end,L_start:L_end)=hdipolez_mc(ctW2,ctL2);
            end
        end
        
        if (0)%hand calculation to verify the calculation of dipole field is correct
            m00=[0.176,0,0];
            m01=[0.384,0,0];
            m10=[0.384,0,0];
            m11=[0.176,0,0];
            
            r00=[0,0,0];
            r10=[-0.022,2.222,0];
            r01=[2.178,0.022,0];
            r11=[2.2,2.2,0];
            
            hdipole00=mu_0/(4*pi)*(3*dot(m10,r10)*r10/(norm(r10)^5)-m10/(norm(r10)^3)+...
                3*dot(m01,r01)*r01/(norm(r01)^5)-m01/(norm(r01)^3)+...
                3*dot(m11,r11)*r11/(norm(r11)^5)-m11/(norm(r11)^3))-mu_0/3*m00/volume_mc;
            
            r00=[-2.178,-0.022,0];
            r10=[-2.2,2.2,0];
            r01=[0,0,0];
            r11=[0.022,2.178,0];
            hdipole01=mu_0/(4*pi)*(3*dot(m10,r10)*r10/(norm(r10))^5-m10/(norm(r10))^3+...
                3*dot(m00,r00)*r00/(norm(r00))^5-m00/(norm(r00))^3+...
                3*dot(m11,r11)*r11/(norm(r11))^5-m11/(norm(r11))^3)-mu_0/3*m01/volume_mc;
            
            hdipole10=mu_0/(4*pi)*(3*dot(m00,r00)*r00/(norm(r00))^5-m00/(norm(r00))^3+...
                3*dot(m01,r01)*r01/(norm(r01))^5-m01/(norm(r01))^3+...
                3*dot(m11,r11)*r11/(norm(r11))^5-m11/(norm(r11))^3)-mu_0/3*m10/volume_mc;
            
            hdipole11=mu_0/(4*pi)*(3*dot(m10,r10)*r10/(norm(r10))^5-m10/(norm(r10))^3+...
                3*dot(m01,r01)*r01/(norm(r01))^5-m01/(norm(r01))^3+...
                3*dot(m00,r00)*r00/(norm(r00))^5-m00/(norm(r00))^3)-mu_0/3*m11/volume_mc;
            %}
            %% the below comparisons should be zero
            
            hdipole00(1)-hdipolex_mc(1,1)
            hdipole00(2)-hdipoley_mc(1,1)
            hdipole00(3)-hdipolez_mc(1,1)
            
            hdipole01(1)-hdipolex_mc(1,2)
            hdipole01(2)-hdipoley_mc(1,2)
            hdipole01(3)-hdipolez_mc(1,2)
            %{
                    hdipole10(1)-hdipolex_mc(2,1)
                    hdipole10(2)-hdipoley_mc(2,1)
                    hdipole10(3)-hdipolez_mc(2,1)
                    
                    hdipole11(1)-hdipolex_mc(2,2)
                    hdipole11(2)-hdipoley_mc(2,2)
                    hdipole11(3)-hdipolez_mc(2,2)
            %}
        end
        
        clear ctW ctL ctW2 ctL2
end

