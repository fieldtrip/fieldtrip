function [H, Z, S, psi] = sfactorization_wilson(S,freq,Niterations,tol,fb,init)

% Usage  : [H, Z, S, psi] = sfactorization_wilson(S,fs,freq);
% Inputs : S (1-sided, 3D-spectral matrix in the form of Channel x Channel x frequency) 
%        : fs (sampling frequency in Hz)
%        : freq (a vector of frequencies) at which S is given
% Outputs: H (transfer function)
%        : Z (noise covariance)
%        : psi (left spectral factor)
% This function is an implemention of Wilson's algorithm (Eq. 3.1)
% for spectral matrix factorization
% Ref: G.T. Wilson,"The Factorization of Matricial Spectral Densities,"
% SIAM J. Appl. Math.23,420-426(1972).
% Written by M. Dhamala & G. Rangarajan, UF, Aug 3-4, 2006.
% Email addresses: mdhamala@bme.ufl.edu, rangaraj@math.iisc.ernet.in

% number of channels
m   = size(S,1);
N   = length(freq)-1;
N2  = 2*N;

% preallocate memory for efficiency
Sarr   = zeros(m,m,N2) + 1i.*zeros(m,m,N2);
gam    = zeros(m,m,N2);
gamtmp = zeros(m,m,N2);
psi    = zeros(m,m,N2);
I      = eye(m); % Defining m x m identity matrix

%Step 1: Forming 2-sided spectral densities for ifft routine in matlab
f_ind = 0;
for f = freq
  f_ind           = f_ind+1;
  Sarr(:,:,f_ind) = S(:,:,f_ind);
  if(f_ind>1)
    Sarr(:,:,2*N+2-f_ind) = S(:,:,f_ind).';
  end
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
      warning('initialization with ''chol'' for iterations did not work well, using arbitrary starting condition');
      tmp = rand(m,m); %arbitrary initial condition
      tmp = triu(tmp);
    end
  case 'rand'
    tmp = rand(m,m); %arbitrary initial condition
    tmp = triu(tmp);
  otherwise
    error('initialization method should be eithe ''chol'' or ''rand''');
end
h = tmp;

for ind = 1:N2
  psi(:,:,ind) = h; 
end

%Step 4: Iterating to get spectral factors
ft_progress('init', fb, 'computing spectral factorization');
for iter = 1:Niterations
  ft_progress(iter./Niterations, 'computing iteration %d/%d\n', iter, Niterations);
  for ind = 1:N2
    invpsi     = inv(psi(:,:,ind));% + I*eps(psi(:,:,ind))); 
    g(:,:,ind) = invpsi*Sarr(:,:,ind)*invpsi'+I;%Eq 3.1
  end
  gp = PlusOperator(g,m,N+1); %gp constitutes positive and half of zero lags 

  psi_old = psi;
  for k = 1:N2
    psi(:,:,k) = psi(:,:,k)*gp(:,:,k);
    psierr(k)  = norm(psi(:,:,k)-psi_old(:,:,k),1);
  end
  psierrf = mean(psierr);
  if(psierrf<tol), 
    fprintf('reaching convergence at iteration %d\n',iter);
    break; 
  end; % checking convergence
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

H = zeros(m,m,N+1) + 1i*zeros(m,m,N+1);
for k = 1:N+1
  H(:,:,k) = psi(:,:,k)*A0inv;       %Transfer function
  S(:,:,k) = psi(:,:,k)*psi(:,:,k)'; %Updated cross-spectral density
end

%---------------------------------------------------------------------
function gp = PlusOperator(g,nchan,nfreq)

% This function is for [ ]+operation: 
% to take the positive lags & half of the zero lag and reconstitute 
% M. Dhamala, UF, August 2006

g   = transpose(reshape(g, [nchan^2 2*(nfreq-1)]));
gam = ifft(g);

% taking only the positive lags and half of the zero lag
gamp  = gam;
beta0 = 0.5*gam(1,:);

gamp(1,          :) = reshape(triu(reshape(beta0, [nchan nchan])),[1 nchan^2]);
gamp(nfreq+1:end,:) = 0;

% reconstituting
gp = fft(gamp);
gp = reshape(transpose(gp), [nchan nchan 2*(nfreq-1)]); 

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
