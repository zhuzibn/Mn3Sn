%% distribute atoms
%1 RE,0 TM
% tmp=randperm(natom,round(natom*compositionn));
% atomtype_=10*ones(natomW,natomL);
% for ctL=1:natomL
%     for ctW=1:natomW
%         atomtype_(ctW,ctL)=ismember(ctW+(ctL-1)*natomW,tmp);
%     end
% end
%% Mn3Sn
atomtype_=10*ones(natomW,natomL);
for ctL=1:natomL
    for ctW=1:natomW
        if mod(ctW,4)==1 && mod(ctL,2)==0
            atomtype_(ctW,ctL)=3;
        elseif mod(ctW,4)==3 && mod(ctL,2)==1
            atomtype_(ctW,ctL)=3;
        elseif mod(ctW,4)==1 && ctL==natomL
            atomtype_(ctW,ctL)=3;
        elseif mod(ctW,4)==3 && ctL==natomL
            atomtype_(ctW,ctL)=3;
        else
            atomtype_(ctW,ctL)=1;
        end
    end
end
clear tmp ctW ctL ctz
%

hani_x = zeros(natomW, natomL);
hani_y = zeros(natomW, natomL);
hani_z = zeros(natomW, natomL);
coeffs_nine_2 = zeros(natomW, natomL);
coeffs_nine_3 = zeros(natomW, natomL);
coeffs_five_2 = zeros(natomW, natomL);
coeffs_five_3 = zeros(natomW, natomL);
coeffs_one_2 = zeros(natomW, natomL);
coeffs_one_3 = zeros(natomW, natomL);
hani_constant = ones(natomW, natomL);

% 计算奇数行/偶数行的系数
coefficient_nine = [0, sin(9*pi/6), cos(9*pi/6)];
coefficient_five = [0, sin(5*pi/6), cos(5*pi/6)];
coefficient_one = [0, sin(1*pi/6), cos(1*pi/6)]; % 此处添加 sin(1*pi/6) 项
coeffs_1 = zeros(natomW, natomL);
coeffs_2 = zeros(natomW, natomL);
coeffs_3 = zeros(natomW, natomL);
for ctW = 1:natomW
    for ctL = 1:natomL
        if atomtype_(ctW, ctL)==1
            if mod(ctW, 4) == 1 || mod(ctW, 4) == 3 %1行3行
                coeffs_nine_2(ctW,ctL) = sin(9*pi/6);
                coeffs_nine_3(ctW,ctL) = cos(9*pi/6);
            elseif mod(ctW, 4) == 2 && mod(ctL, 2) == 0 % 2行偶数列
                coeffs_one_2(ctW,ctL) = sin(1*pi/6);
                coeffs_one_3(ctW,ctL) = cos(1*pi/6);
            elseif mod(ctW, 4) == 0 && mod(ctL, 2) == 1 % 4行奇数列
                coeffs_one_2(ctW,ctL) = sin(1*pi/6);
                coeffs_one_3(ctW,ctL) = cos(1*pi/6);
            elseif mod(ctW, 4) == 2 && mod(ctL, 2) == 1 % 2行奇数列
                coeffs_five_2(ctW,ctL) = sin(5*pi/6);
                coeffs_five_3(ctW,ctL) = cos(5*pi/6);
            elseif mod(ctW, 4) == 0 && mod(ctL, 2) == 0 % 4行偶数列
                coeffs_five_2(ctW,ctL) = sin(5*pi/6);
                coeffs_five_3(ctW,ctL) = cos(5*pi/6);
            else
            end
        end
    end
end
coeffs_2 = coeffs_nine_2+coeffs_one_2+coeffs_five_2;
coeffs_3 = coeffs_nine_3+coeffs_one_3+coeffs_five_3;

% if ctW>2 && ctL>2
%     for ctW = 2:natomW-1
%         for ctL = 2:natomL-1
%             hani_constant(ctW,ctL) = 1;
%         end
%     end
% end

if load_fixed_atom_distrib
    clear atomtype_
    load('atomtypee.mat');
elseif save_fixed_atom_distrib%save debug data
    save('atomtypee.mat','atomtype_');%change this to save(ddebugfilename);
    error('distribution mat file has been saved, run the program again by setting load_fixed_atom_distrib=1')
end