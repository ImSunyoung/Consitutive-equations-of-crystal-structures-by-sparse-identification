
%----------------------------------------------------------
% get_invariant
% 
%   try_01 : 9-parameter (1st order term)
%   try_02 : 9-parameter (2nd order term)
%   try_03 : 9-parameter (3rd order term) --> no change (same as try_02)
% 
%----------------------------------------------------------

% ---------------------------------------------------
function [Ji] = get_invariant_tetra2(matG, derivative_order)
% ---------------------------------------------------
    % voigt notation = [11 22 33 23 31 12]
    J1 = [ 
        1 0 0 0 0 0;
        0 1 0 0 0 0;
        ];
    J2 = [
        0 0 1 0 0 0;
        ];
    J3 = [
        0 0 0 2 0 0;
        0 0 0 0 2 0;
        ];
    J4 = [
        0 0 0 0 0 2;
        ];
    J5 = [
        1 1 0 0 0 0;
        ];
    J6 = [
        0 0 0 1 1 1;
        ];
    J7 = [
        1 0 0 2 0 0;
        0 1 0 0 2 0;
        ];
    J8 = [
        0 0 0 2 2 0;
        ];
    % ---------------------------------------------------
    JJ{1} = J1;
    JJ{2} = J2;
    JJ{3} = J3;
    JJ{4} = J4;
    JJ{5} = J5;
    JJ{6} = J6;
    JJ{7} = J7;
    JJ{8} = J8;
    % ---------------------------------------------------
    vecG = [
        matG(1,1), matG(2,2), matG(3,3), matG(2,3), matG(3,1), matG(1,2);
        ];        
    % ---------------------------------------------------
    switch(derivative_order)
        case 1
            Ji = get_Ji_G(JJ,vecG);
        case 2
            Ji = get_Ji_GG(JJ,vecG);
        otherwise
            Ji = get_Ji(JJ,vecG);
    end
end

% ---------------------------------------------------
function [Ji] = get_Ji(JJ,vecG)
% ---------------------------------------------------
    nJJ = length(JJ);
    for aa=1:nJJ
        [nrow, ncol] = size(JJ{aa});
        Jaa = JJ{aa};
        
        JG = 0;
        for rr=1:nrow
            Jrr = Jaa(rr,:);
            vecGJ = vecG .^ Jrr;
            
            f_idx = find(Jrr);
            if (~isempty(f_idx))
                JG = JG + prod(vecGJ(f_idx));
%                 vecGJ
%                 vecG
%                 Jrr
            end
        end
        Ji{aa} = JG;
    end
% ---------------------------------------------------
end

% ---------------------------------------------------
function [Ji_G] = get_Ji_G(JJ,vecG)
% ---------------------------------------------------
    % Polynomial derivative
    nJJ = length(JJ);
    for aa=1:nJJ
        [nrow, ncol] = size(JJ{aa});
        Jaa = JJ{aa};
        
        JG = zeros(6,1);
        for ii=1:6
            JGi = 0;
            for rr=1:nrow
                Jrr = Jaa(rr,:);
                vecGrr = vecG;
                coeff_GJ = 1;
                
%                 check_GJ = Jrr .* vecGrr;
                f_idx = find(Jrr);
                % ---------------------------------------------------
                if (0 == Jrr(ii))
                    continue;
                else
                    coeff_GJ = coeff_GJ * Jrr(ii);
                    Jrr(ii) = Jrr(ii) - 1;
                end
                % ---------------------------------------------------
                vecGJ = vecGrr .^ Jrr;
                if (~isempty(f_idx))
                    JGi = JGi + coeff_GJ * prod(vecGJ(f_idx));
                    if (1 == aa)
                        vecGJ;
                        vecG;
                        Jrr;
                    end
                end
            end
            JG(ii) = JGi;
        end
        Ji_G{aa} = JG;
    end
    % ---------------------------------------------------
end

% ---------------------------------------------------
function [Ji_GG] = get_Ji_GG(JJ,vecG)
% ---------------------------------------------------
    % Polynomial derivative
    nJJ = length(JJ);
    for aa=1:nJJ
        [nrow, ncol] = size(JJ{aa});
        Jaa = JJ{aa};
        
        JG = zeros(6,6);
        for ii=1:6
            for jj=ii:6
                JGi = 0;
                for rr=1:nrow
                    % ---------------------------------------------------
                    Jrr = Jaa(rr,:);
                    vecGrr = vecG;
                    coeff_GJ = 1;

%                     check_GJ = Jrr .* vecGrr;
                    f_idx = find(Jrr);
                    % ---------------------------------------------------
                    if (0 == Jrr(ii) || 0 == Jrr(jj))
                        continue;
                    end
                    if ((ii == jj) && (Jrr(ii) <= 1))
                        continue;
                    end
                    % ---------------------------------------------------
                    if (Jrr(ii) > 0)
                        coeff_GJ = coeff_GJ * Jrr(ii);
                        Jrr(ii) = Jrr(ii) - 1;
                    end
                    if (Jrr(jj) > 0)
                        coeff_GJ = coeff_GJ * Jrr(jj);
                        Jrr(jj) = Jrr(jj) - 1;
                    end
                    % ---------------------------------------------------
                    vecGJ = vecGrr .^ Jrr;
                    if (~isempty(f_idx))
                        JGi = JGi + coeff_GJ * prod(vecGJ(f_idx));
                    end
%                     Jrr
                    % ---------------------------------------------------
                    if (4 == (ii) && 4 == (jj))
                        if (JGi)
%                             ijr = [aa,ii,jj,rr]
%                             JGi
                        end
                    end
                    % ---------------------------------------------------
                end
                JG(ii,jj) = JGi;
                JG(jj,ii) = JGi;
            end
        end
%         aa
%         JG
        Ji_GG{aa} = JG;
    end
    % ---------------------------------------------------
end
% ---------------------------------------------------