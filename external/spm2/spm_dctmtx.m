function C = spm_dctmtx(N,K,n,f)
% Creates basis functions for Discrete Cosine Transform.
% FORMAT C = spm_dctmtx(N,K,n)
%     OR C = spm_dctmtx(N,K)
%     OR D = spm_dctmtx(N,K,n,'diff')
%     OR D = spm_dctmtx(N,K,'diff')
% N - dimension
% K - order
% n - optional points to sample
%____________________________________________________________________________
% spm_dctmtx creates a matrix for the first few basis functions of a one
% dimensional discrete cosine transform.
% With the 'diff' argument, spm_dctmtx produces the derivatives of the
% DCT.
%
% See:    Fundamentals of Digital Image Processing (p 150-154).
%         Anil K. Jain 1989.
%____________________________________________________________________________
% @(#)spm_dctmtx.m	2.1 John Ashburner MRCCU/FIL 01/08/28

d = 0;

if nargin == 1, K = N; end;

if any(nargin == [1 2]),
	n = (0:(N-1))';
elseif nargin == 3,
	if strcmp(n,'diff'),
		d = 1;
		n = (0:(N-1))';
	elseif strcmp(n,'diff2'),
		d = 2;
		n = (0:(N-1))';
	else
		n = n(:);
	end
elseif nargin == 4,
	n = n(:);
	if strcmp(f,'diff'),
		d = 1;
	elseif strcmp(n,'diff2'),
		d = 2;
	else
		error('Incorrect Usage');
	end
else
	error('Incorrect Usage');
end

C = zeros(size(n,1),K);

if d == 0,
	C(:,1)=ones(size(n,1),1)/sqrt(N);
	for k=2:K
		C(:,k) = sqrt(2/N)*cos(pi*(2*n+1)*(k-1)/(2*N));
	end
elseif d == 1,
	for k=2:K
		C(:,k) = -2^(1/2)*(1/N)^(1/2)*sin(1/2*pi*(2*n*k-2*n+k-1)/N)*pi*(k-1)/N;
	end
elseif d == 2,
	for k=2:K,
		C(:,k) = -2^(1/2)*(1/N)^(1/2)*cos(1/2*pi*(2*n+1)*(k-1)/N)*pi^2*(k-1)^2/N^2;
	end;
else,
	error('Can''t do this');
end


