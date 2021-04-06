function [var_out] = func_cubic_predict(varargin)
%----------------------------------------------------------
fid             = varargin{1};
material_coeff	= varargin{2};
matF            = varargin{3};
%----------------------------------------------------------
opt_speed = 1;
%----------------------------------------------------------
ncoeff = length(material_coeff);
temp = find(material_coeff(:,2) == 0);
ncoeff_1 = length(temp);
temp = find(material_coeff(:,3) == 0);
ncoeff_2 = length(temp) - ncoeff_1;
ncoeff_3 = ncoeff - ncoeff_1 - ncoeff_2;

ci = material_coeff(:,4);
%----------------------------------------------------------
matC = matF'*matF;
matG = matC - eye(3);
%----------------------------------------------------------

%     JJ = get_invariant_cubic_func(matG,0);
%     JJ_C = get_invariant_cubic_func(matG,1);
%     JJ_CC = get_invariant_cubic_func(matG,2);

JJ_temp = get_invariant_cubic(matG,0);
JJ_C_temp = get_invariant_cubic(matG,1);
JJ_CC_temp = get_invariant_cubic(matG,2);
JJ = zeros(1,9);
JJ_C = zeros(6,9);
JJ_CC = zeros(6,6,9);
for tp = 1: 9
    JJ(tp) = JJ_temp{tp};
    JJ_C(:,tp) = JJ_C_temp{tp};
    JJ_CC(:,:,tp) = JJ_CC_temp{tp};
end


%----------------------------------------------------------
matL_PK2 = zeros(6, ncoeff);
%----------------------------------------------------------
% 1st-order terms
for kk=1:ncoeff_1
    aa = material_coeff(kk,1);
    
    Ja_C = JJ_C(:,aa);
    
    
    if (opt_speed)
        t1 = 1:3;
        t2 = 4:6;
        matL_PK2(t1,kk) = Ja_C(t1);
        matL_PK2(t2,kk) = Ja_C(t2) / 2;
    else
        for tt=1:6
            matL_PK2(tt,kk) = Ja_C(tt);
            if (tt > 3)
                matL_PK2(tt,kk) = matL_PK2(tt,kk) / 2;
            end
        end
    end
end
%----------------------------------------------------------
% 2nd-order terms
for kk=1:ncoeff_2
    ab = kk + ncoeff_1;
    aa = material_coeff(ab,1);
    bb = material_coeff(ab,2);
    
    
    Ja = JJ(aa);
    Ja_C = JJ_C(:,aa);
    
    Jb = JJ(bb);
    Jb_C = JJ_C(:,bb);
    
    
    if (opt_speed)
        tt = 1:6;
        matL_PK2(tt,ab) = Ja_C(tt)*Jb + Ja*Jb_C(tt);
        tt = 4:6;
        matL_PK2(tt,ab) = matL_PK2(tt,ab) / 2;
    else
        for tt=1:6
            matL_PK2(tt,ab) = Ja_C(tt)*Jb + Ja*Jb_C(tt);
            if (tt > 3)
                matL_PK2(tt,ab) = matL_PK2(tt,ab) / 2;
            end
        end
    end
end
%----------------------------------------------------------
% 3rd-order terms
for kk=1:ncoeff_3
    abc = kk + ncoeff_1 + ncoeff_2;
    aa = material_coeff(abc,1);
    bb = material_coeff(abc,2);
    cc = material_coeff(abc,3);
    
    
    Ja = JJ(aa);
    Ja_C = JJ_C(:,aa);
    
    Jb = JJ(bb);
    Jb_C = JJ_C(:,bb);
    
    Jc = JJ(cc);
    Jc_C = JJ_C(:,cc);
    
    
    if (opt_speed)
        tt = 1:6;
        matL_PK2(tt,abc) = Ja_C(tt)*Jb*Jc + Ja*Jb_C(tt)*Jc + Ja*Jb*Jc_C(tt);
        tt = 4:6;
        matL_PK2(tt,abc) = matL_PK2(tt,abc) / 2;
    else
        for tt=1:6
            matL_PK2(tt,abc) = Ja_C(tt)*Jb*Jc + Ja*Jb_C(tt)*Jc + Ja*Jb*Jc_C(tt);
            if (tt > 3)
                matL_PK2(tt,abc) = matL_PK2(tt,abc) / 2;
            end
        end
    end
end
%----------------------------------------------------------
pk = 2 * matL_PK2 * ci;
PK2 = [
    pk(1), pk(6), pk(5);
    pk(6), pk(2), pk(4);
    pk(5), pk(4), pk(3);
    ];
%----------------------------------------------------------
if (1 == fid)
    var_out = cell(1);
    var_out{1} = PK2;
    return;
else
    %----------------------------------------------------------
    %----------------------------------------------------------
    PE_CC = zeros(6,6);
    %----------------------------------------------------------
    % 1st-order terms
    for kk=1:ncoeff_1
        aa = material_coeff(kk,1);
        
        
        Ja_CC = JJ_CC(:,:,aa);
        
        
        Ja_CC(4:6,:) = Ja_CC(4:6,:) / 2;
        Ja_CC(:,4:6) = Ja_CC(:,4:6) / 2;
        
        PE_CC = PE_CC + ci(kk) * Ja_CC;
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
        
        PE_CC = PE_CC + ci(ab) * JaJb_CC;
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
        
        
        %     JaJbJc_CC = Ja_CC*Jb*Jc + Ja_C*Jb_C'*Jc + Ja_C*Jb*Jc_C' ...
        %               + Jb_C*Ja_C'*Jc + Ja*Jb_CC*Jc + Ja*Jb_C*Jc_C' ...
        %               + Ja_C*Jb*Jc_C' + Ja*Jc_C*Jb_C' + Ja*Jb*Jc_CC;
        JaJbJc_CC = Ja_CC*Jb*Jc + Ja_C*Jb_C'*Jc + Ja_C*Jb*Jc_C' ...
            + Jb_C*Ja_C'*Jc + Ja*Jb_CC*Jc + Ja*Jb_C*Jc_C' ...
            + Jc_C*Jb*Ja_C' + Ja*Jc_C*Jb_C' + Ja*Jb*Jc_CC;
        
        JaJbJc_CC(4:6,:) = JaJbJc_CC(4:6,:) / 2;
        JaJbJc_CC(:,4:6) = JaJbJc_CC(:,4:6) / 2;
        
        PE_CC = PE_CC + ci(abc) * JaJbJc_CC;
    end
    %----------------------------------------------------------
    Cijkl = 4 * PE_CC;
    %----------------------------------------------------------
    var_out = cell(1,2);
    var_out{1} = PK2;
    var_out{2} = Cijkl;
end
%----------------------------------------------------------