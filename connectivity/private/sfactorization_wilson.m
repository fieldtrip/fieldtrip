function [H, Z, S, psi] = sfactorization_wilson(S,freq,Niterations,tol,fb,init,checkflag,stabilityfix)

% SFACTORIZATION_WILSON performs multivariate non-parametric spectral factorization on
% cross-spectra, based on Wilson's algorithm.
%
% Usage  : [H, Z, S, psi] = sfactorization_wilson(S,freq);
%
% Inputs : S (1-sided, 3D-spectral matrix in the form of Channel x Channel x frequency) 
%        : freq (a vector of frequencies) at which S is given. 
%
% Outputs: H (transfer function)
%        : Z (noise covariance)
%        : psi (left spectral factor)
%
% This function is an implemention of Wilson's algorithm (Eq. 3.1)
% for spectral matrix factorization
%
% Ref: G.T. Wilson,"The Factorization of Matricial Spectral Densities,"
% SIAM J. Appl. Math.23,420-426(1972).
% Written by M. Dhamala & G. Rangarajan, UF, Aug 3-4, 2006.
% Email addresses: mdhamala@bme.ufl.edu, rangaraj@math.iisc.ernet.in

% Copyright (C) 2009-2013, Jan-Mathijs Schoffelen
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

if nargin<8, stabilityfix = false; end
if nargin<7, checkflag = true;   end
if nargin<6, init      = 'chol'; end
if nargin<5, fb        = 'none'; end
if nargin<4, tol       = 1e-8;   end
if nargin<3, Niterations = 1000; end

dfreq = round(diff(freq)*1e5)./1e5; % allow for some numeric issues
if ~all(dfreq==dfreq(1))
  ft_error('FieldTrip:connectivity:sfactorization_wilson', 'frequency axis is not evenly spaced');
end

if freq(1)~=0
  ft_warning('FieldTrip:connectivity:sfactorization_wilson', 'when performing non-parametric spectral factorization, the frequency axis should ideally start at 0, zero padding the spectral density'); 
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

% check whether the last frequency bin is strictly real-valued.
% if that's the case, then it is assumed to be the Nyquist frequency
% and the two-sided spectral density will have an even number of 
% frequency bins. if not, in order to preserve hermitian symmetry,
% the number of frequency bins needs to be odd.
Send = S(:,:,end);
N    = numel(freq);
m    = size(S,1);
if all(imag(Send(:))<abs(trace(Send)./size(Send,1)*1e-9))
  N2 = 2*(N-1);
else
  N2 = 2*(N-1)+1;
end

% preallocate memory for efficiency
Sarr   = zeros(m,m,N2) + 1i.*zeros(m,m,N2);
gam    = zeros(m,m,N2);
gamtmp = zeros(m,m,N2);
psi    = zeros(m,m,N2);
I      = eye(m); % Defining m x m identity matrix

%Step 1: Forming 2-sided spectral densities for ifft routine in matlab

% the input cross-spectral density is assumed to be weighted with a
% factor of 2 in all non-DC and Nyquist bins, therefore weight the 
% DC-bin with a factor of 2 to get a correct two-sided representation
Sarr(:,:,1) = S(:,:,1).*2;
for f_ind = 2:N
  Sarr(:,:,       f_ind) = S(:,:,f_ind);
  Sarr(:,:,(N2+2)-f_ind) = S(:,:,f_ind).';
end

% the input cross-spectral density is assumed to be weighted with a
% factor of 2 in all non-DC and Nyquist bins, therefore weight the 
% Nyquist bin with a factor of 2 to get a correct two-sided representation
if mod(size(Sarr,3),2)==0
  Sarr(:,:,N) = Sarr(:,:,N).*2;
end

%Step 2: Computing covariance matrices
for k1 = 1:m
  for k2 = 1:m
    %gam(k1,k2,:) = real(ifft(squeeze(Sarr(k1,k2,:)))*fs); %FIXME think about this
    gam(k1,k2,:) = real(ifft(squeeze(Sarr(k1,k2,:))));
  end
end

%Step 3: Initializing for iterations 
gam0 = gam(:,:,1);
switch init
  case 'chol'
    [tmp, dum] = chol(gam0);
    if dum
      ft_warning('initialization with ''chol'' for iterations did not work well, using arbitrary starting condition');
      tmp = randn(m,1000); %arbitrary initial condition
      tmp = (tmp*tmp')./1000;
      %tmp = triu(tmp);
      [tmp, dum] = chol(tmp);
    end
  case 'rand'
    %tmp = randn(m,m); %arbitrary initial condition
    %tmp = triu(tmp);
    tmp = randn(m,1000); %arbitrary initial condition
    tmp = (tmp*tmp')./1000;
    %tmp = triu(tmp);
    [tmp, dum] = chol(tmp);
  otherwise
    ft_error('initialization method should be eithe ''chol'' or ''rand''');
end
h = tmp;

for ind = 1:N2
  psi(:,:,ind) = h; 
end

%Step 4: Iterating to get spectral factors
g = zeros(size(psi));
ft_progress('init', fb, 'computing spectral factorization');
for iter = 1:Niterations
  ft_progress(iter./Niterations, 'computing iteration %d/%d\n', iter, Niterations);
  for ind = 1:N2
    %invpsi     = inv(psi(:,:,ind)); 
    g(:,:,ind) = psi(:,:,ind)\Sarr(:,:,ind)/psi(:,:,ind)'+I;%Eq 3.1
  end
  gp = PlusOperator(g,m,N,stabilityfix); %gp constitutes positive and half of zero lags 

  psi_old = psi;
  for k = 1:N2
    psi(:,:,k) = psi(:,:,k)*gp(:,:,k);
    psierr(k)  = norm(psi(:,:,k)-psi_old(:,:,k),1);
  end

  if checkflag,
    psierrf = mean(psierr);
    if(psierrf<tol), 
      fprintf('reaching convergence at iteration %d\n',iter);
      break; 
    end % checking convergence
  end
end 
ft_progress('close');

%Step 5: Getting covariance matrix from spectral factors
for k1 = 1:m
  for k2 = 1:m
    gamtmp(k1,k2,:) = real(ifft(squeeze(psi(k1,k2,:))));
  end
end

%Step 6: Getting noise covariance & transfer function (see Example pp. 424)
A0    = gamtmp(:,:,1); 
A0inv = inv(A0);

%Z     = A0*A0.'*fs; %Noise covariance matrix
Z     = A0*A0.'; %Noise covariance matrix not multiplied by sampling frequency

%FIXME check this; at least not multiplying it removes the need to correct later on
%this also makes it more equivalent to the noisecov estimated by biosig's mvar-function

H = zeros(m,m,N) + 1i*zeros(m,m,N);
for k = 1:N
  H(:,:,k) = psi(:,:,k)*A0inv;       %Transfer function
  S(:,:,k) = psi(:,:,k)*psi(:,:,k)'; %Updated cross-spectral density
end

if numel(selfreq)~=numel(freq)
  % return only the frequency bins that were in the input
  H   =   H(:,:,selfreq);
  S   =   S(:,:,selfreq);
  psi = psi(:,:,selfreq);
end

%---------------------------------------------------------------------
function gp = PlusOperator(g,nchan,nfreq,stabilityfix)

% This function is for [ ]+operation: 
% to take the positive lags & half of the zero lag and reconstitute 
% M. Dhamala, UF, August 2006

g   = transpose(reshape(g, nchan^2, []));
gam = ifft(g);

% taking only the positive lags and half of the zero lag
gamp  = gam;
beta0 = 0.5*gam(1,:);

gamp(1,          :) = reshape(triu(reshape(beta0, [nchan nchan])),[1 nchan^2]);
gamp(nfreq+1:end,:) = 0;

% smooth with a window, only for the long latency boundary: this is a
% stabilityfix proposed by Martin Vinck
if stabilityfix
  w = tukeywin(nfreq*2, 0.5);
  gamp(1:nfreq,:) = gamp(1:nfreq,:).*repmat(w(nfreq+1:end),[1 nchan^2]);
else
  % nothing to be done here  
end

% reconstituting
gp = fft(gamp);
gp = reshape(transpose(gp), [nchan nchan numel(gp)/(nchan^2)]); 

%------------------------------------------------------
%this is the original code; above is vectorized version
%which is assumed to be faster with many channels present
%for k1 = 1:nchan
%  for k2 = 1:nchan
%    gam(k1,k2,:) = ifft(squeeze(g(k1,k2,:)));
%  end
%end
%
%% taking only the positive lags and half of the zero lag
%gamp  = gam;
%beta0 = 0.5*gam(:,:,1); 
%gamp(:,:,1) = triu(beta0);  %this is Stau
%gamp(:,:,nfreq+1:end) = 0;
%
%% reconstituting
%for k1 = 1:nchan
%  for k2 = 1:nchan
%    gp(k1,k2,:) = fft(squeeze(gamp(k1,k2,:)));
%  end
%end
