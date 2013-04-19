function x = mvm_fft(pj, n, fftKcirc, H, B_m, b)
%MVM_FFT Fast matrix vector multiplication using FFT for Logistic-Gaussian
%        Process density model
%
%  See also
%    GPLA_ND_E, LGPDENS, DEMO_LGPDENS
%
% Copyright (c) 2012 Jaakko Riihimäki

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

pjsq=sqrt(pj);
Rb=sqrt(n)*(-pj*(pj'*(b./pjsq))+pjsq.*b);

if size(fftKcirc,1)==1 % 1D case
  
  n1=size(pj,1);
  v=[Rb; zeros(n1,1)];
  q=ifft(fftKcirc.*fft(v'));
  q=q(1:n1)';
  
else % 2D case
  
  n1=size(fftKcirc,2)/2;
  n2=size(fftKcirc,1)/2;
  
  v=zeros(2*n2,2*n1);
  v(1:n2,1:n1)=reshape(Rb,n2,n1);
  q=ifft2(fftKcirc.*fft2(v));
  q=q(1:n2,1:n1);
  q=q(:);
  
end

if ~isempty(H)
  q=q+H'*(B_m*(H*Rb));
end
x=b+sqrt(n)*( pjsq.*q - (pj*(pj'*q))./pjsq );