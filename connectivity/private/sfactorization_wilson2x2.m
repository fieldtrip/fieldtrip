function [H, Z, S, psi] = sfactorization_wilson2x2(S,freq,Niterations,tol,cmbindx,fb,init,checkflag)

% SFACTORIZATION_WILSON2X2 performs pairwise non-parametric spectral factorization on
% cross-spectra, based on Wilson's algorithm.
%
% Usage  : [H, Z, psi] = sfactorization_wilson(S,freq);
%
% Inputs : S (1-sided, 3D-spectral matrix in the form of Channel x Channel x frequency) 
%        : freq (a vector of frequencies) at which S is given. 
%
% Outputs: H (transfer function)
%        : Z (noise covariance)
%        : S (cross-spectral density 1-sided)
%        : psi (left spectral factor)
%
% This function is an implemention of Wilson's algorithm (Eq. 3.1)
% for spectral matrix factorization.
%
% Ref: G.T. Wilson,"The Factorization of Matricial Spectral Densities,"
% SIAM J. Appl. Math.23,420-426(1972).
% Written by M. Dhamala & G. Rangarajan, UF, Aug 3-4, 2006.
% Email addresses: mdhamala@bme.ufl.edu, rangaraj@math.iisc.ernet.in

% Copyright (C) 2009-2017, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if nargin<8, checkflag = true;   end
if nargin<7, init      = 'chol'; end
if nargin<6, fb        = 'none'; end
if nargin<5
  error('FieldTrip:connectivity:sfactorization_wilson2x2', 'when requesting multiple pairwise spectral decomposition, ''cmbindx'' needs to be specified');
end
if nargin<4, tol        = 1e-8;   end
if nargin<3, Niterations = 1000;  end;

dfreq = round(diff(freq)*1e5)./1e5; % allow for some numeric issues
if ~all(dfreq==dfreq(1))
  error('FieldTrip:connectivity:sfactorization_wilson2x2', 'frequency axis is not evenly spaced');
end

if freq(1)~=0
  ft_warning('FieldTrip:connectivity:sfactorization_wilson2x2', 'when performing non-parametric spectral factorization, the frequency axis should ideally start at 0, zero padding the spectral density'); 
  dfreq = mean(dfreq);
  npad  = freq(1)./dfreq;
  
  % update the freq axis and keep track of the frequency bins that are
  % expected in the output
  selfreq  = (1:numel(freq)) + npad;
  freq     = [(0:(npad-1))./dfreq freq];
  S        = cat(3, zeros(size(S,1), size(S,1), npad), S);  
else
  selfreq  = 1:numel(freq);
end

% ensure input S is double (mex-files don't work with single)
S = double(S);

% check whether the last frequency bin is strictly real-valued.
% if that's the case, then it is assumed to be the Nyquist frequency
% and the two-sided spectral density will have an even number of 
% frequency bins. if not, in order to preserve hermitian symmetry,
% the number of frequency bins needs to be odd.
Send = S(:,:,end);
N    = numel(freq);
m    = size(cmbindx,1);
if all(imag(Send(:))<abs(trace(Send)./size(Send,1)*1e-9))
  hasnyq = true;
  N2     = 2*(N-1);
else
  hasnyq = false;
  N2     = 2*(N-1)+1;
end

% preallocate memory for the identity matrix
I      = repmat(eye(2),[1 1 m N2]); % Defining 2 x 2 identity matrix

% %Step 1: Forming 2-sided spectral densities for ifft routine in matlab
% Sarr   = zeros(2,2,m,N2) + 1i.*zeros(2,2,m,N2);
% for c = 1:m
%   Stmp  = S(cmbindx(c,:),cmbindx(c,:),:);
%   
%   % the input cross-spectral density is assumed to be weighted with a
%   % factor of 2 in all non-DC and Nyquist bins, therefore weight the 
%   % DC-bin with a factor of 2 to get a correct two-sided representation
%   Sarr(:,:,c,1) = Stmp(:,:,1).*2;
%   
%   for f_ind = 2:N
%     Sarr(:,:,c,       f_ind) = Stmp(:,:,f_ind);
%     Sarr(:,:,c,(N2+2)-f_ind) = Stmp(:,:,f_ind).';
%   end
% end
% Sarr2 = Sarr;

% preallocate memory for the 2-sided spectral density
Sarr = zeros(2,2,N2,m) + 1i.*zeros(2,2,N2,m);
for c = 1:m
  Sarr(:,:,1:N,c) = S(cmbindx(c,:),cmbindx(c,:),:);
end
if hasnyq
  N1 = N;
else
  N1 = N + 1; % the highest frequency needs to be represented twice, for symmetry purposes
end
Sarr(:,:,N1:N2,:) = flip(Sarr(:,:,2:N,:),3);
Sarr(2,1,N1:N2,:) = conj(Sarr(2,1,N1:N2,:));
Sarr(1,2,N1:N2,:) = conj(Sarr(1,2,N1:N2,:));
Sarr              = permute(Sarr, [1 2 4 3]);
Sarr(:,:,:,1)     = Sarr(:,:,:,1).*2; % weight the DC-bin


% the input cross-spectral density is assumed to be weighted with a
% factor of 2 in all non-DC and Nyquist bins, therefore weight the 
% Nyquist bin with a factor of 2 to get a correct two-sided representation
if hasnyq
  Sarr(:,:,:,N) = Sarr(:,:,:,N).*2;
end

%Step 2: Computing covariance matrices
gam = real(reshape(ifft(reshape(Sarr, [4*m N2]), [], 2),[2 2 m N2]));

%Step 3: Initializing for iterations 
gam0 = gam(:,:,:,1);

h    = complex(zeros(size(gam0)));
for k = 1:m
  switch init
    case 'chol'
      [tmp, dum] = chol(gam0(:,:,k));
      if dum
        warning('initialization with ''chol'' for iterations did not work well, using arbitrary starting condition');
        tmp = rand(2,2); %arbitrary initial condition
        tmp = triu(tmp);
      end
    case 'rand'
      tmp = rand(2,2); %arbitrary initial condition
      tmp = triu(tmp);
    otherwise
      error('initialization method should be eithe ''chol'' or ''rand''');
  end
  h(:,:,k) = tmp;
  
  %h(:,:,k) = chol(gam0(:,:,k));
end
psi  = repmat(h, [1 1 1 N2]);

%Step 4: Iterating to get spectral factors
ft_progress('init', fb, 'computing spectral factorization');
for iter = 1:Niterations
  ft_progress(iter./Niterations, 'computing iteration %d/%d\n', iter, Niterations);
  invpsi = inv2x2(psi);
  g      = sandwich2x2(invpsi, Sarr) + I;
  gp     = PlusOperator2x2(g,m,N); %gp constitutes positive and half of zero lags 
  
  psi_old = psi;
  psi     = mtimes2x2(psi, gp);
  
  if checkflag
    psierr  = abs(psi-psi_old)./abs(psi);
    psierrf = mean(psierr(:));
    if(psierrf<tol) 
      fprintf('reaching convergence at iteration %d\n',iter);
      break; 
    end % checking convergence
  end
end 
ft_progress('close');

%Step 5: Getting covariance matrix from spectral factors
gamtmp = reshape(real(ifft(transpose(reshape(psi, [4*m N2]))))', [2 2 m N2]);

%Step 6: Getting noise covariance & transfer function (see Example pp. 424)
A0    = gamtmp(:,:,:,1); 
A0inv = inv2x2(A0);

Z = zeros(2,2,m);
for k = 1:m
  %Z     = A0*A0.'*fs; %Noise covariance matrix
  Z(:,:,k) = A0(:,:,k)*A0(:,:,k).'; %Noise covariance matrix not multiplied by sampling frequency
  %FIXME check this; at least not multiplying it removes the need to correct later on
  %this also makes it more equivalent to the noisecov estimated by biosig's mvar-function
end

% H = complex(zeros(2,2,m,N));
% S = complex(zeros(2,2,m,N));
% for k = 1:N
%   for kk = 1:m
%     H(:,:,kk,k) = psi(:,:,kk,k)*A0inv(:,:,kk);  % Transfer function
%     S(:,:,kk,k) = psi(:,:,kk,k)*psi(:,:,kk,k)'; % Cross-spectral density
%   end
% end
H = mtimes2x2(psi,A0inv(:,:,:,ones(1,size(psi,4))));
S = mtimes2x2(psi,ctranspose2x2(psi));

siz = [size(H) 1 1];
H   = reshape(H, [4*siz(3) siz(4:end)]);
siz = [size(S) 1 1];
S   = reshape(S, [4*siz(3) siz(4:end)]);
siz = [size(Z) 1 1];
Z   = reshape(Z, [4*siz(3) siz(4:end)]);
siz = [size(psi) 1 1];
psi = reshape(psi, [4*siz(3) siz(4:end)]);

%if numel(selfreq)~=numel(freq)
  % return only the frequency bins that were in the input
  H   =   H(:,selfreq,:,:);
  S   =   S(:,selfreq,:,:);
  psi = psi(:,selfreq,:,:);
%end
  
%---------------------------------------------------------------------
function gp = PlusOperator2x2(g,ncmb,nfreq)

% This function is for [ ]+operation: 
% to take the positive lags & half of the zero lag and reconstitute 
% M. Dhamala, UF, August 2006

g   = transpose(reshape(g, 4*ncmb, []));
gam = ifft(g);

% taking only the positive lags and half of the zero lag
gamp  = gam;
beta0 = 0.5*gam(1,:);

%for k = 1:ncmb
%  gamp(1,(k-1)*4+1:k*4) = reshape(triu(reshape(beta0(1,(k-1)*4+1:k*4),[2 2])),[1 4]);
%end
beta0(2:4:4*ncmb)   = 0;
gamp(1,:)           = beta0;
gamp(nfreq+1:end,:) = 0;

% reconstituting
gp = fft(gamp);
gp = reshape(transpose(gp), [2 2 ncmb numel(gp)/(4*ncmb)]); 
