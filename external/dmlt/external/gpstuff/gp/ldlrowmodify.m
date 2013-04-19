function L = ldlrowmodify(L, c2, k)
%LDLROWMODIFY  Function to modify the sparse cholesky factorization
%              L*D*L' = C, when a row and column k of C have changed
%
%  Description
%    L = rowmodify(L, c2, k)
%    Function to modify the Cholesky factorization L*D*L' = C
%    (which is stored in L), when a row and column k of C have
%    changed from c to c2. The change in C is assumed to be such
%    that the sparsity structure of C is remained the same. c is
%    the old kth column of C and c2 is the new column.
%   
%    See GPEP_E for usage
%
%  Reference  
%    See Davis and Hager 2005 (Row Modification of a sparse
%    Cholesky factorization) section 4 for details of the
%    algorithm.
    
%  Copyright (c) 2009-2010 Jarno Vanhatalo
    
% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

% First evaluate the l_21 vector. That is the kth row vector of L
    
% Evaluate the old c. NOTE! The C-version (mex-file) of this algorithm
% does not need the old c. However, this m-version is an old
% implementation
    [L2, D] = ldlsplit(L);
    c = L2(k,:)*D(k)*L2(k,:)';
    
    n = size(L,1);
    if k>1
        deltac12 = c2(1:k-1)-c(1:k-1);
        [L2, D] = ldlsplit(L);

        %deltal12 = D(1:k-1,1:k-1) \ (L2(1:k-1,1:k-1)\deltac12);
        
        x = deltac12;
        [I,J,ld] = find(L);
        a2 = cumsum(histc(J,1:n));
        Lp=zeros(n,1);
        Lp(1) = 1; Lp(2:end) = a2(1:end-1) + 1; Lp(end+1) = a2(end);
        Li = I;
        Lx = ld;

        for j=1:k-1
            x(j) = x(j);
            for p = Lp(j)+1:Lp(j+1)-1
                if Li(p) >= k
                    break
                end
                x(Li(p)) = x(Li(p))-Lx(p).*x(j);
            end
            x(j) = x(j)./Lx(Lp(j));
        end
        
        deltal12 = x;
        
        l12 = L(k,1:k-1)';
        l12b = l12 + deltal12;
        L(k,1:k-1) = l12b;
        
        d = D(k,k);
        db = d + c2(k)-c(k) - sum( deltal12.*(l12 + l12b).*diag(D(1:k-1,1:k-1)));
        L(k,k) = db;       
    else
        d = L(k,k);
        db = d + c2(k)-c(k);
        L(k,k) = db;
    end
    
    if k < n
        deltac32 = c2(k+1:n)-c(k+1:n);    
        l32 = L(k+1:n,k);
        if k>1
            l32b = (deltac32 + l32.*d - L(k+1:n,1:k-1)*D(1:k-1,1:k-1)*deltal12)./db;
        else
            l32b = (deltac32 + l32.*d )./db;
        end
        
        w = L(:,k).*sqrt(d);
        L(k+1:n,k) = l32b;
        w2 =L(:,k).*sqrt(db);

        
        % Then update the L_33 part. That is L(k+1:n, k+1:n)
        alpha = 1;
        alpha2 = 1;
        for i = k+1:n    % Li(Lp(k)+1:Lp(k+1)-1)
            
            if w(i) ~= 0
                beta = alpha + w(i).^2 / L(i,i); 
                gamma = w(i) / (beta*L(i,i));
                L(i,i) = (beta/alpha)*L(i,i);
                alpha = beta;
                
                beta2 = alpha2 - w2(i).^2 /L(i,i);
                gamma2 = w2(i) / (beta2*L(i,i));        
                L(i,i) = (beta2/alpha2)*L(i,i);
                alpha2 = beta2;
                
            
                if (i == n) return, end
                
                w(i+1:n) = w(i+1:n) - w(i).*L(i+1:n, i);
                L(i+1:n,i) = L(i+1:n,i) + gamma.*w(i+1:n);
                
                w2(i+1:n) = w2(i+1:n) - w2(i).*L(i+1:n, i);
                L(i+1:n,i) = L(i+1:n,i) - gamma2.*w2(i+1:n);
            end
        end
    end
end