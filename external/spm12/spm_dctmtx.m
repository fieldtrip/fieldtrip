function [C,D] = spm_dctmtx(N,K,n,f)
% Create basis functions for Discrete Cosine Transform
% FORMAT C = spm_dctmtx(N)
% FORMAT C = spm_dctmtx(N,K)
% FORMAT C = spm_dctmtx(N,K,n)
% FORMAT D = spm_dctmtx(N,K,'diff')
% FORMAT D = spm_dctmtx(N,K,n,'diff')
% FORMAT D = spm_dctmtx(N,K,'diff',dx)
%
% N        - dimension
% K        - order
% n        - optional points to sample
%
% C        - DCT matrix or its derivative
%__________________________________________________________________________
%
% spm_dctmtx creates a matrix for the first few basis functions of a one
% dimensional discrete cosine transform.
% With the 'diff' argument, spm_dctmtx produces the derivatives of the DCT.
%
% If N and K are vectors, C is a large prod(N) x prod(K) matrix
% corresponding to the Kronecker tensor product of each N-dimensional
% basis set. This is useful for dealing with vectorised N-arrays. An
% additional argument, dx can be specified to scale the derivatives
%
% Reference:
% Fundamentals of Digital Image Processing (p 150-154). Anil K. Jain, 1989.
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 1996-2022 Wellcome Centre for Human Neuroimaging


% Kroneckor form for n-dimensional DCTs
%--------------------------------------------------------------------------
nN  = numel(N);
if nN > 1
    if nargin > 1
        if numel(K) < nN
            K = repmat(K,1,nN);
        end
    else
        K = N;
    end
    
    if nargin < 3
        
        % Kroneckor form
        %------------------------------------------------------------------
        C = 1;
        for i = 1:nN
            c = spm_dctmtx(N(i),K(i));
            C = kron(c,C);
        end
        
    else
        
        % intervals
        %------------------------------------------------------------------
        if nargin < 4
            f = ones(nN,1);
        end
        
        % derivatives (Kroneckor form)
        %------------------------------------------------------------------
        D     = cell(nN,1);
        for d = 1:nN
            C = 1;
            for i = 1:nN
                if i == d
                    c = spm_dctmtx(N(i),K(i),'diff');
                else
                    c = spm_dctmtx(N(i),K(i));
                end
                C = kron(c,C);
            end
            D{d}  = C/f(d);
        end
        C  = spm_cat(D);
        
    end

    return
end

% DCT matrix
%--------------------------------------------------------------------------
d = 0;

if nargin == 1, K = N; end

if nargin < 3
    n = (0:(N-1))';
elseif nargin == 3
    if strcmp(n,'diff')
        d = 1;
        n = (0:(N-1))';
    elseif strcmp(n,'diff2')
        d = 2;
        n = (0:(N-1))';
    else
        n = n(:);
    end
elseif nargin == 4
    n = n(:);
    if strcmp(f,'diff')
        d = 1;
    elseif strcmp(f,'diff2')
        d = 2;
    else
        error('Incorrect Usage.');
    end
else
    error('Incorrect Usage.');
end

C = zeros(size(n,1),K);

if d == 0
    C(:,1)     = ones(size(n,1),1)/sqrt(N);
    for k=2:K
        C(:,k) = sqrt(2/N)*cos(pi*(2*n+1)*(k-1)/(2*N));
    end
elseif d == 1
    for k=2:K
        C(:,k) = -2^(1/2)*(1/N)^(1/2)*sin(1/2*pi*(2*n*k-2*n+k-1)/N)*pi*(k-1)/N;
    end
elseif d == 2
    for k=2:K
        C(:,k) = -2^(1/2)*(1/N)^(1/2)*cos(1/2*pi*(2*n+1)*(k-1)/N)*pi^2*(k-1)^2/N^2;
    end
else
    error('Incorrect Usage.');
end
