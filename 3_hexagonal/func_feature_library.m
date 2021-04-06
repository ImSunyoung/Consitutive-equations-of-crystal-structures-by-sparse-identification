function [ A4, refer_tmodulus, material_coeff ] = func_feature_library(refer, train_loading_set)
nCP =  9;
adCP = [1 2 3 7 8 12 16 19 21];
refer_matF = [];
refer_tmodulus = [];
for i = train_loading_set
    refer_matF = [refer_matF;refer{i,5}];
    nPt = size(refer{i,5},1);
    for j = 1:nPt
        temp = refer{i,4}(j,adCP)';
        refer_tmodulus = [refer_tmodulus;temp];
    end
end
nPoint = length(refer_matF);
nData = length(refer_tmodulus);

%% User define material coefficients 
material_coeff = zeros(1,4);
iter = 0;
for i = 1:9
    iter = iter + 1;
    material_coeff(iter,1) = i;
end
for i = 1:9
    for j = i:9
        iter = iter + 1;
        material_coeff(iter,1:2) = [i j];
    end
end
for i = 1:9
    for j = i:9
        for k = j:9
            iter = iter + 1;
            material_coeff(iter,1:3) = [i j k];
        end
    end
end


% order computation 
od_invariant = [1 1 2 2 3 3 4 5 6];
od_trans = zeros(size(material_coeff,1),3);
od_trans(:,1) = od_invariant(material_coeff(:,1));
od_trans(material_coeff(:,2)>0,2) = od_invariant(material_coeff(material_coeff(:,2)>0,2));
od_trans(material_coeff(:,3)>0,3) = od_invariant(material_coeff(material_coeff(:,3)>0,3));

od_sum = sum(od_trans,2);
material_coeff_new = material_coeff(od_sum<5,:);
%----------------------------------------------------------
ncoeff = length(material_coeff_new);
temp = find(material_coeff_new(:,2) == 0);
ncoeff_1 = length(temp);
temp = find(material_coeff_new(:,3) == 0);
ncoeff_2 = length(temp) - ncoeff_1;
ncoeff_3 = ncoeff - ncoeff_1 - ncoeff_2;
material_coeff = material_coeff_new;

%% The compoenet of coefficient matrix A4 
A = zeros(nData,ncoeff);
for i = 1:nPoint
    coeff = 0;
    temp = refer_matF(i,:);
    matF = reshape(temp,3,[]);
    %----------------------------------------------------------
    matC = matF'*matF;
    matG = matC - eye(3);
    %----------------------------------------------------------
    JJ_temp = get_invariant_hexa4(matG,0);        % scalar      ~ 9
    JJ_C_temp = get_invariant_hexa4(matG,1);      % vector 6*1  ~ 9
    JJ_CC_temp = get_invariant_hexa4(matG,2);     % matrix 6*6  ~ 9
    JJ = zeros(1,9);
    JJ_C = zeros(6,9);
    JJ_CC = zeros(6,6,9);
    for tp = 1: 9
        JJ(tp) = JJ_temp{tp};
        JJ_C(:,tp) = JJ_C_temp{tp};
        JJ_CC(:,:,tp) = JJ_CC_temp{tp};
    end
        %----------------------------------------------------------
    % 1st-order terms
    for kk=1:ncoeff_1
        aa = material_coeff(kk,1);
        
        Ja_CC = JJ_CC(:,:,aa);
        
        Ja_CC(4:6,:) = Ja_CC(4:6,:) / 2;
        Ja_CC(:,4:6) = Ja_CC(:,4:6) / 2;
        
        %%% Full A %%%
        mat = Ja_CC;
        coeff = coeff + 1; 
        A((i-1)*nCP+1:i*nCP, coeff) = [mat(1,1) mat(1,2) mat(1,3)...
                                                mat(2,2) mat(2,3)...
                                                         mat(3,3)...
                                       mat(4,4) mat(5,5) mat(6,6)]';
    end
    %----------------------------------------------------------
    % 2nd-order terms
    for kk=1:ncoeff_2
        ab = kk + ncoeff_1;
        aa = material_coeff(ab,1);
        bb = material_coeff(ab,2);
        
        
        Ja = JJ(aa);
        Ja_C = JJ_C(:,aa);
        Ja_CC = JJ_CC(:,:,aa);
        
        Jb = JJ(bb);
        Jb_C = JJ_C(:,bb);
        Jb_CC = JJ_CC(:,:,bb);
        
        
        JaJb_CC = Ja_CC*Jb + Ja_C*Jb_C' + Jb_C*Ja_C' + Ja*Jb_CC;
        JaJb_CC(4:6,:) = JaJb_CC(4:6,:) / 2;
        JaJb_CC(:,4:6) = JaJb_CC(:,4:6) / 2;
        
        %%% Full A %%%
        mat = JaJb_CC;
        coeff = coeff + 1; 
        A((i-1)*nCP+1:i*nCP, coeff) = [mat(1,1) mat(1,2) mat(1,3)...
                                                mat(2,2) mat(2,3)...
                                                         mat(3,3)...
                                       mat(4,4) mat(5,5) mat(6,6)]';
    end
    %     PE_CC
    %----------------------------------------------------------
    % 3rd-order terms
    for kk=1:ncoeff_3
        abc = kk + ncoeff_1 + ncoeff_2;
        aa = material_coeff(abc,1);
        bb = material_coeff(abc,2);
        cc = material_coeff(abc,3);
        
        
        Ja = JJ(aa);
        Ja_C = JJ_C(:,aa);
        Ja_CC = JJ_CC(:,:,aa);
        
        Jb = JJ(bb);
        Jb_C = JJ_C(:,bb);
        Jb_CC = JJ_CC(:,:,bb);
        
        Jc = JJ(cc);
        Jc_C = JJ_C(:,cc);
        Jc_CC = JJ_CC(:,:,cc);

        JaJbJc_CC = Ja_CC*Jb*Jc + Ja_C*Jb_C'*Jc + Ja_C*Jb*Jc_C' ...
            + Jb_C*Ja_C'*Jc + Ja*Jb_CC*Jc + Ja*Jb_C*Jc_C' ... 
            + Jc_C*Jb*Ja_C' + Ja*Jc_C*Jb_C' + Ja*Jb*Jc_CC;
        
        JaJbJc_CC(4:6,:) = JaJbJc_CC(4:6,:) / 2;
        JaJbJc_CC(:,4:6) = JaJbJc_CC(:,4:6) / 2;
        
        %%% Full A %%%
        mat = JaJbJc_CC;
        coeff = coeff + 1; 
        A((i-1)*nCP+1:i*nCP, coeff) = [mat(1,1) mat(1,2) mat(1,3)...
                                                mat(2,2) mat(2,3)...
                                                         mat(3,3)...
                                       mat(4,4) mat(5,5) mat(6,6)]';
    end
    %----------------------------------------------------------
    A4 = 4 * A;
    %----------------------------------------------------------
    
end
end

