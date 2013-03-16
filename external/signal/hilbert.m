function [x] = hilbert(xr)
% Computes analytic signal
% FORMAT [x] = hilbert(xr)
%
% Returns analytic signal x = xr + i*xi such that 
% xi is the Hilbert transform of real vector xr.
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging
 
% Will Penny
% $Id$

if ~isreal(xr)
  xr = real(xr);
end

% Work along the first nonsingleton dimension
[xr,nshifts] = shiftdim(xr);

n = size(xr,1);
x = fft(xr,n,1); % n-point FFT over columns.
h  = zeros(n,~isempty(x)); % nx1 for nonempty. 0x0 for empty.
if n > 0 && 2*fix(n/2) == n
  % even and nonempty
  h([1 n/2+1]) = 1;
  h(2:n/2) = 2;
elseif n>0
  % odd and nonempty
  h(1) = 1;
  h(2:(n+1)/2) = 2;
end
x = ifft(x.*h(:,ones(1,size(x,2))));

% Convert back to the original shape.
x = shiftdim(x,-nshifts);
