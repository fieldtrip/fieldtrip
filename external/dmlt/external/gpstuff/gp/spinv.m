function Z = spinv(A, B)
% SPINV    Evaluate the sparse inverse matrix
%
% z = sinv(A)  returns the elements of inv(A)_ij, for which A_ij
%      is different from zero. 
%
% z = sinv(LD, 1)  returns the elements of inv(A)_ij, for which A_ij
%     is different from zero, and where LD is the LDL cholesky
%     decomposition of A. LD has to be in the form returned by ldlchol
%     in SuiteSparse by Tim Davis. 
%
%   Note! If z = sinv(LD, 1) is used LD must not be modified in Matlab
%   after ldlchol. Matlab destroys the symbolic sparsity structure in
%   the Cholesky decomposition, which is needed in the spinv
%   algorithm. If LD is modified the worst scenario is memory corruption. 
%   
%
%     For details, see:
%     Jarno Vanhatalo and Aki Vehtari (2008). Modelling local and
%     global phenomena with sparse Gaussian processes. Proceedings of
%     the 24th Conference on Uncertainty in Artificial Intelligence
%

% Copyright (c) 2008-2010      Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.
    
    n = size(A,1);

    if nargin == 1
        [LD, p, q] = ldlchol(A);
    else 
        LD = A;
    end
    
    [I,J,ld] = find(LD);
    temp = [I(:) J(:) ; J(:) I(:)];
    temp = sortrows(unique(temp,'rows'),2);
    Iz = temp(:,1); Jz = temp(:,2); 
    
    % Find the column starting points
    a1=zeros(n,1);
    a2 = cumsum(histc(J,1:n));
    a1(1) = 1; a1(2:end) = a2(1:end-1) + 1;
    az1=zeros(n,1);
    az2 = cumsum(histc(Jz,1:n));
    az1(1) = 1; az1(2:end) = az2(1:end-1) + 1;
    
    for j=1:n
        indaz{j} = az1(j):az2(j);
        indIz{j} = Iz(indaz{j})';
    end

    % Evaluate the sparse inverse
    z = zeros(size(Iz));
    z(end) = 1./ld(end);
    % Allocate memory
    cindit=zeros(n,1);
    for jj = n-1:-1:1
        fil = ld(a1(jj)+1:a1(jj+1)-1);
        fi = I(a1(jj)+1:a1(jj+1)-1);
        lfi = length(fi);
        Zt = zeros(lfi,lfi);
        indz = cumsum(histc(indIz{jj},[0 ; fi]));
        indz = az1(jj) + indz(1:end-1);
        
        i4=0;            
        for i1 = 1:lfi
            cind1=indaz{fi(i1)};
            Icind1=indIz{fi(i1)};
            indfi = lfi;
            i2=length(Icind1);
            go = true;
            while go
                if Icind1(i2)==jj  % Find the indeces for the jj'th rows in fi columns
                    i4=i4+1;
                    cindit(i4)=cind1(i2);
                    go = false;
                end
                if indfi >= 1 && fi(indfi) == Icind1(i2) % Find the indeces for the fi'th rows in i2'nd columns
                    Zt(indfi,i1) = z(cind1(i2));
                    indfi = indfi-1;
                end
                i2 = i2-1;
            end
        end
        % remove extras
        cindi=cindit(1:i4);

        zij = -fil'*Zt;
        z(cindi) = zij;
        z(indz) = zij;
        zij = 1./ld(a1(jj)) - fil'*z(indz);
        z(az1(jj)-1+find(indIz{jj}==jj,1)) = zij;
    end
    
    Z = sparse(Iz,Jz,z);
    
    if nargin == 1
        r(q) = 1:n;
        Z = Z(r,r);
    end
end
