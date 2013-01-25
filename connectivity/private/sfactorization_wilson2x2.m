function [H, Z, S, psi] = sfactorization_wilson2x2(S,freq,Niterations,tol,cmbindx,fb,init)

% Usage  : [H, Z, psi] = sfactorization_wilson(S,fs,freq);
% Inputs : S (1-sided, 3D-spectral matrix in the form of Channel x Channel x frequency) 
%        : fs (sampling frequency in Hz)
%        : freq (a vector of frequencies) at which S is given
% Outputs: H (transfer function)
%        : Z (noise covariance)
%        : S (cross-spectral density 1-sided)
%        : psi (left spectral factor)
% This function is an implemention of Wilson's algorithm (Eq. 3.1)
% for spectral matrix factorization
% Ref: G.T. Wilson,"The Factorization of Matricial Spectral Densities,"
% SIAM J. Appl. Math.23,420-426(1972).
% Written by M. Dhamala & G. Rangarajan, UF, Aug 3-4, 2006.
% Email addresses: mdhamala@bme.ufl.edu, rangaraj@math.iisc.ernet.in


m   = size(cmbindx,1);
N   = length(freq)-1;
N2  = 2*N;

% preallocate memory for efficiency
Sarr   = zeros(2,2,m,N2) + 1i.*zeros(2,2,m,N2);
I      = repmat(eye(2),[1 1 m N2]); % Defining 2 x 2 identity matrix

%Step 1: Forming 2-sided spectral densities for ifft routine in matlab
for c = 1:m
  % f_ind = 0;
  Stmp  = S(cmbindx(c,:),cmbindx(c,:),:);
  for f_ind = 1:(N+1)
  % for f = freq
    % f_ind             = f_ind+1;
    Sarr(:,:,c,f_ind) = Stmp(:,:,f_ind);
    if(f_ind>1)
      Sarr(:,:,c,2*N+2-f_ind) = Stmp(:,:,f_ind).';
    end
  end
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
  gp     = PlusOperator2x2(g,m,N+1); %gp constitutes positive and half of zero lags 
  
  psi_old = psi;
  psi     = mtimes2x2(psi, gp);
  %psierr  = sum(sum(abs(psi-psi_old)));
  psierr  = abs(psi-psi_old)./abs(psi);
  
  if 0
    plot(squeeze(psierr(2,1,1,:))); hold on
    plot(squeeze(psierr(1,1,1,:)),'r');drawnow
  end
  psierrf = mean(psierr(:));
  if(psierrf<tol), 
    fprintf('reaching convergence at iteration %d\n',iter);
    break; 
  end; % checking convergence
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

H = complex(zeros(2,2,m,N+1));
S = complex(zeros(2,2,m,N+1));
for k = 1:(N+1)
  for kk = 1:m
    H(:,:,kk,k) = psi(:,:,kk,k)*A0inv(:,:,kk);  % Transfer function
    S(:,:,kk,k) = psi(:,:,kk,k)*psi(:,:,kk,k)'; % Cross-spectral density
  end
end

siz = [size(H) 1 1];
H   = reshape(H, [4*siz(3) siz(4:end)]);
siz = [size(S) 1 1];
S   = reshape(S, [4*siz(3) siz(4:end)]);
siz = [size(Z) 1 1];
Z   = reshape(Z, [4*siz(3) siz(4:end)]);
siz = [size(psi) 1 1];
psi = reshape(psi, [4*siz(3) siz(4:end)]);

%---------------------------------------------------------------------
function gp = PlusOperator2x2(g,ncmb,nfreq)

% This function is for [ ]+operation: 
% to take the positive lags & half of the zero lag and reconstitute 
% M. Dhamala, UF, August 2006

g   = transpose(reshape(g, [4*ncmb 2*(nfreq-1)]));
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
gp = reshape(transpose(gp), [2 2 ncmb 2*(nfreq-1)]); 
