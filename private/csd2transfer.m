function [output] = csd2transfer(freq, varargin)

% CSD2TRANSFER computes the transfer-function from frequency domain data
% using the Wilson-Burg algorithm. The transfer function can be used for
% the computation of directional measures of connectivity, such as granger
% causality, partial directed coherence, or directed transfer functions
%
% Use as
%   [output] = csd2transfer(freq, varargin)
%
% Where freq is a data structure containing frequency domain data containing
% the cross-spectral density computed between all pairs of channels, thus
% containing a 'dimord' of 'chan_chan_freq(_time)'.
%
% Configurable options come in key-value pairs:
%
%   numiteration = scalar value (default: 100) the number of iterations
%   channelcmb   = Nx2 cell-array listing the channel pairs for the spectral
%                    factorization. If not defined or empty (default), a
%                    full multivariate factorization is performed, otherwise
%                    a multiple pairwise factorization is done.
%   tol          = scalar value (default: 1e-18) tolerance limit truncating
%                    the iterations
%
% The code for the Wilson-Burg algorithm has been very generously provided by 
% Dr. Mukesh Dhamala, and Prof. Mingzhou Ding and his group. 
%
% If you use this code for studying directed interactions, please cite from
% the following references:
% -M.Dhamala, R.Rangarajan, M.Ding, Physical Review Letters 100, 018701 (2008)
% -M.Dhamala, R.rangarajan, M.Ding, Neuroimage 41, 354 (2008)

% Undocumented options:
%
%   block
%   blockindx
%   svd
%
% Copyright (C) 2009-2011, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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
% $Id: csd2transfer $

numiteration = ft_getopt(varargin, 'numiteration', 100);
channelcmb   = ft_getopt(varargin, 'channelcmb',   {});
block        = ft_getopt(varargin, 'block',        {});
blockindx    = ft_getopt(varargin, 'blockindx',    cell(0,2));
tol          = ft_getopt(varargin, 'tol',          1e-18);
fb           = ft_getopt(varargin, 'feedback',     'textbar');
sfmethod     = ft_getopt(varargin, 'sfmethod',     'multivariate');
dosvd        = ft_getopt(varargin, 'svd',          'no');

dosvd = istrue(dosvd);

% check whether input data is valid
freq = ft_checkdata(freq, 'datatype', 'freq');
if ~isfield(freq, 'crsspctrm') || ~isfield(freq, 'label')
  error('the input data does not contain cross-spectral density data in the supported format');
end

hasrpt = ~isempty(strfind(freq.dimord, 'rpt'));
if hasrpt,
  nrpt = numel(freq.cumtapcnt);
else
  nrpt = 1;
end

% if bivariate is requested without channelcmb, do all versus all pairwise
if isempty(channelcmb) && strcmp(sfmethod, 'bivariate')
  channelcmb = {'all' 'all'};
end
  
if ~isempty(channelcmb)
  channelcmb = ft_channelcombination(channelcmb, freq.label);
end

if ~isempty(block)
  % sanity check 1
  if ~iscell(block) 
    error('block should be a cell-array containing 3 cells');
  end
  % sanity check 2
  if strcmp(sfmethod, 'bivariate')
    error('when block is specified, it is not OK to do bivariate decomposition');
  end
end

%fsample = findcfg(freq.cfg, 'fsample');
fsample = 0;

if isfield(freq, 'time'),
  ntim = numel(freq.time);
else
  ntim = 1;
end

siz = size(freq.crsspctrm);
if ntim==1,
  siz = [siz 1]; %add dummy dimensionality for time axis
end

if strcmp(sfmethod, 'bivariate')
  fprintf('computing pairwise non-parametric spectral factorization on %d channel pairs\n', size(channelcmb,1) - numel(unique(channelcmb(:))));
elseif strcmp(sfmethod, 'multivariate')
  fprintf('computing multivariate non-parametric spectral factorization on %d channels\n', numel(freq.label));   
else
  error('unknown sfmethod %s', sfmethod);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computations start here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(sfmethod, 'multivariate') && isempty(block) && nrpt>1,
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % multiple repetitions, loop over repetitions
  % multivariate decomposition
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  H = zeros(siz) + 1i.*zeros(siz);
  S = zeros(siz) + 1i.*zeros(siz);
  Z = zeros([siz(1:3) siz(end)]);
  for k = 1:nrpt
    for m = 1:ntim
      tmp = reshape(freq.crsspctrm(k,:,:,:,m), siz(2:end-1));
      [Htmp, Ztmp, Stmp] = sfactorization_wilson(tmp, fsample, freq.freq, ...
                                                   numiteration, tol, fb);
      H(k,:,:,:,m) = Htmp;
      Z(k,:,:,m)   = Ztmp;
      S(k,:,:,:,m) = Stmp;
    end 
  end
elseif strcmp(sfmethod, 'multivariate') && isempty(block),
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % standard code
  % multivariate decomposition
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  H   = zeros(siz) + 1i.*zeros(siz);
  S   = zeros(siz) + 1i.*zeros(siz);
  Z   = zeros([siz(1:2) siz(end)]);
  
  % only do decomposition once
  for m = 1:ntim
    tmp = freq.crsspctrm(:,:,:,m);

    % do SVD to avoid zigzags due to numerical issues
    if dosvd
      dat     = sum(tmp,3);
      [u,s,v] = svd(real(dat));
      for k = 1:size(tmp,3)
        tmp(:,:,k) = u'*tmp(:,:,k)*u;
      end
    end
    
    if any(isnan(tmp(:))),
      Htmp = nan;
      Ztmp = nan;
      Stmp = nan;
    else
      [Htmp, Ztmp, Stmp] = sfactorization_wilson(tmp, fsample, freq.freq, ...
                                                   numiteration, tol, fb);
    end
    
    % undo SVD
    if dosvd
      for k = 1:size(tmp,3)
        Htmp(:,:,k) = u*Htmp(:,:,k)*u';
        Stmp(:,:,k) = u*Stmp(:,:,k)*u';
      end
      Ztmp = u*Ztmp*u';
    end
    
    H(:,:,:,m) = Htmp;
    Z(:,:,m)   = Ztmp;
    S(:,:,:,m) = Stmp;
  end
elseif strcmp(sfmethod, 'bivariate') && nrpt>1,
  % error 
  error('single trial estimates and linear combination indexing is not implemented');
elseif nrpt>1 && ~isempty(block),
  % error 
  error('single trial estimates and blockwise factorisation is not yet implemented');
elseif strcmp(sfmethod, 'bivariate')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % pairwise factorization resulting in linearly indexed transfer functions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %convert list of channel labels into indices
  cmbindx     = zeros(size(channelcmb));
  ok          = true(size(cmbindx,1), 1);
  for k = 1:size(cmbindx,1)
    [tmp, cmbindx(k,:)] = match_str(channelcmb(k,:)', freq.label);
    if ~isempty(intersect(cmbindx(1:(k-1),:), cmbindx(k,:), 'rows'))
      ok(k) = false;
    elseif cmbindx(k,1)==cmbindx(k,2)
      ok(k) = false;
    end
  end
  
  %remove auto-combinations and double occurrences
  cmbindx    = cmbindx(ok,:);
  channelcmb = channelcmb(ok,:);
  
 
  %do multiple 2x2 factorization efficiently
  if ntim>1,
    for kk = 1:ntim
      [Htmp, Ztmp, Stmp] = sfactorization_wilson2x2(freq.crsspctrm(:,:,:,kk), fsample, ...
                               freq.freq, numiteration, tol, cmbindx, fb);
      if kk==1,
        H   = Htmp;
        Z   = Ztmp;
        S   = Stmp;
        H(:,:,2:ntim) = nan;
        Z(:,:,2:ntim) = nan;
        S(:,:,2:ntim) = nan;
      else
        H(:,:,kk) = Htmp;
        Z(:,:,kk) = Ztmp;
        S(:,:,kk) = Stmp;
      end
    end
  else
    % if the number of pairs becomes too big, it seems to slow down quite a
    % bit. try to chunk
    if size(cmbindx,1)>1000
      begchunk = 1:1000:size(cmbindx,1);
      endchunk = [1000:1000:size(cmbindx,1) size(cmbindx,1)];
      H = zeros(4*size(cmbindx,1), numel(freq.freq));
      S = zeros(4*size(cmbindx,1), numel(freq.freq));
      Z = zeros(4*size(cmbindx,1), 1);
      for k = 1:numel(begchunk)
        fprintf('computing factorization of chunck %d/%d\n', k, numel(begchunk));
        [Htmp, Ztmp, Stmp] = sfactorization_wilson2x2(freq.crsspctrm, fsample, freq.freq, ...
                                             numiteration, tol, cmbindx(begchunk(k):endchunk(k),:), fb);
                                           
        begix = (k-1)*4000+1;
        endix = min(k*4000, size(cmbindx,1)*4);
        H(begix:endix, :) = Htmp;
        S(begix:endix, :) = Stmp;
        Z(begix:endix, :) = Ztmp;
                                           
      end
    else
      [H, Z, S] = sfactorization_wilson2x2(freq.crsspctrm, fsample, freq.freq, ...
                                             numiteration, tol, cmbindx, fb);
    end
  end
  
  %convert crsspctrm accordingly
  siz          = [size(H) 1];
  tmpcrsspctrm = complex(zeros([2 2 siz(1)/4 siz(2:end)]));
  for k = 1:size(cmbindx,1)
    tmpcrsspctrm(:,:,k,:,:) = freq.crsspctrm(cmbindx(k,:),cmbindx(k,:),:,:);   
  end
  freq.crsspctrm = reshape(tmpcrsspctrm, siz);
  
  labelcmb = cell(size(cmbindx,1)*4, 2);
  for k = 1:size(cmbindx,1)
    indx = (k-1)*4 + (1:4);
    labelcmb{indx(1),1} = [channelcmb{k,1},'[',channelcmb{k,1},channelcmb{k,2},']'];
    labelcmb{indx(1),2} = [channelcmb{k,1},'[',channelcmb{k,1},channelcmb{k,2},']'];
    labelcmb{indx(2),1} = [channelcmb{k,2},'[',channelcmb{k,1},channelcmb{k,2},']'];
    labelcmb{indx(2),2} = [channelcmb{k,1},'[',channelcmb{k,1},channelcmb{k,2},']'];
    labelcmb{indx(3),1} = [channelcmb{k,1},'[',channelcmb{k,1},channelcmb{k,2},']'];
    labelcmb{indx(3),2} = [channelcmb{k,2},'[',channelcmb{k,1},channelcmb{k,2},']'];
    labelcmb{indx(4),1} = [channelcmb{k,2},'[',channelcmb{k,1},channelcmb{k,2},']'];
    labelcmb{indx(4),2} = [channelcmb{k,2},'[',channelcmb{k,1},channelcmb{k,2},']'];
  end
elseif ~isempty(block)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % blockwise multivariate stuff
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if ntim>1,
    error('blockwise factorization of tfrs is not yet possible');
  end
  
  %reorder the channel order such that the blocks are ordered
  nblocks = unique(blockindx{2});
  for k = 1:numel(nblocks)
    %b{k} = cfg.blockindx{2}(find(cfg.blockindx{2}==nblocks(k)));
    b{k} = find(blockindx{2}==nblocks(k));
  end
  indx = cat(1,b{:});
  freq.crsspctrm = freq.crsspctrm(indx, indx, :);
  freq.label     = freq.label(indx);
  bindx          = blockindx{2}(indx);

  for k = 1:numel(block)
    sel  = find(ismember(bindx, block{k}));
    Stmp = freq.crsspctrm(sel,sel,:);
    
    % do PCA to avoid zigzags due to numerical issues
    dopca = 1;
    if dopca
      dat     = sum(Stmp,3);
      [u,s,v] = svd(real(dat));
      for m = 1:size(Stmp,3)
        Stmp(:,:,m) = u'*Stmp(:,:,m)*u;
      end
    end
    
    [Htmp, Ztmp, Stmp] = sfactorization_wilson(Stmp, fsample, freq.freq, ...
                                                 numiteration, tol, fb);  
    
    % undo PCA
    if dopca
      for m = 1:size(Stmp,3)
        Htmp(:,:,m) = u*Htmp(:,:,m)*u';
        Stmp(:,:,m) = u*Stmp(:,:,m)*u';
      end
      Ztmp = u*Ztmp*u';
    end
                                               
    siz  = [size(Htmp) 1]; siz = [siz(1)*siz(2) siz(3:end)];
    Htmp = reshape(Htmp, siz);
    siz  = [size(Ztmp) 1]; siz = [siz(1)*siz(2) siz(3:end)];
    Ztmp = reshape(Ztmp, siz);
    siz  = [size(Stmp) 1]; siz = [siz(1)*siz(2) siz(3:end)];
    Stmp = reshape(Stmp, siz);
    
    tmpindx = [];
    cmbtmp  = cell(siz(1), 2);
    [tmpindx(:,1), tmpindx(:,2)] = ind2sub(sqrt(siz(1))*[1 1],1:siz(1));
    for kk = 1:size(cmbtmp,1)
      cmbtmp{kk,1} = [freq.label{sel(tmpindx(kk,1))},'[',cat(2,freq.label{sel}),']'];
      cmbtmp{kk,2} = [freq.label{sel(tmpindx(kk,2))},'[',cat(2,freq.label{sel}),']'];
    end
    
    %concatenate
    if k == 1,
      H = Htmp;
      Z = Ztmp;
      S = Stmp;
      labelcmb = cmbtmp;
    else
      H = cat(1,H,Htmp);
      Z = cat(1,Z,Ztmp);
      S = cat(1,S,Stmp);
      labelcmb = cat(1,labelcmb,cmbtmp);
    end 
  end
  
  if strcmp(freq.dimord(1:9), 'chan_chan'),
    freq.dimord = ['chancmb_',freq.dimord(10:end)];
  end
end

% create output
output           = [];
if strcmp(sfmethod, 'multivariate')
  output.dimord    = freq.dimord;
else
  ix = strfind(freq.dimord, 'chan_chan');
  newdimord = [freq.dimord(1:(ix+3)),'cmb',freq.dimord(ix+9:end)];

  %dimtok = tokenize(freq.dimord, '_');
  %chdim  = 0;
  %newdimord = '';
  %for k = 1:numel(dimtok)
  %  if strcmp(dimtok{k}, 'chan'), chdim = chdim+1; end
  %  if chdim==1,
  %    newdimord = [newdimord,'_chancmb'];
  %  elseif chdim
  %    newdimord = [newdimord, '_', dimtok{k}];
  %  end
  %end
  output.dimord    = newdimord;
end
output.label     = freq.label;
output.freq      = freq.freq;
output.crsspctrm = S;
output.transfer  = H;
output.noisecov  = Z;
try
  output.time = freq.time;
catch
end
try
  output.cumtapcnt = freq.cumtapcnt;
catch
end
try 
  output.cumsumcnt = freq.cumsumcnt;
catch
end
if exist('labelcmb', 'var') && ~isempty(labelcmb)
  output.labelcmb = labelcmb;
  output          = rmfield(output, 'label');
end

%-------------------------------------------------------------------
function [H, Z, S, psi] = sfactorization_wilson(S,fs,freq,Niterations,tol,fb)

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
[h, dum] = chol(gam0);
if dum
  warning('initialization for iterations did not work well, using arbitrary starting condition');
  h = rand(m,m); h = triu(h); %arbitrary initial condition
end

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

%------------------------------------------------------------------------------
function [H, Z, S, psi] = sfactorization_wilson2x2(S,fs,freq,Niterations,tol,cmbindx,fb)

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
  f_ind = 0;
  Stmp  = S(cmbindx(c,:),cmbindx(c,:),:);
  for f = freq
    f_ind             = f_ind+1;
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
  [tmp, dum] = chol(gam0(:,:,k));
  if dum
    warning('initialization for iterations did not work well, using arbitrary starting condition');
    tmp = rand(2,2); h(:,:,k) = triu(tmp); %arbitrary initial condition
  else
    h(:,:,k) = tmp;
  end
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
  psierr  = max(sum(abs(psi-psi_old)));
  
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

for k = 1:ncmb
  gamp(1,(k-1)*4+1:k*4) = reshape(triu(reshape(beta0(1,(k-1)*4+1:k*4),[2 2])),[1 4]);
end
gamp(nfreq+1:end,:) = 0;

% reconstituting
gp = fft(gamp);
gp = reshape(transpose(gp), [2 2 ncmb 2*(nfreq-1)]); 

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subfunctions for the 2x2 functionality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%determinant of a 2x2 matrix
function [d] = det2x2(x)

%computes determinant of matrix x, using explicit analytic definition if
%size(x,1) < 4, otherwise use matlab det-function

siz = size(x);
if all(siz(1:2)==2),
  d = x(1,1,:,:).*x(2,2,:,:) - x(1,2,:,:).*x(2,1,:,:);
elseif all(siz(1:2)==3),
  d = x(1,1,:,:).*x(2,2,:,:).*x(3,3,:,:) - ...
      x(1,1,:,:).*x(2,3,:,:).*x(3,2,:,:) - ...
      x(1,2,:,:).*x(2,1,:,:).*x(3,3,:,:) + ...
      x(1,2,:,:).*x(2,3,:,:).*x(3,1,:,:) + ...
      x(1,3,:,:).*x(2,1,:,:).*x(3,2,:,:) - ...
      x(1,3,:,:).*x(2,2,:,:).*x(3,1,:,:);
elseif numel(siz)==2,
  d = det(x);
else
  %error   
  %write for loop
  %for
  %end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% inverse of a 2x2 matrix
function [d] = inv2x2(x)

%computes inverse of matrix x, using explicit analytic definition if
%size(x,1) < 4, otherwise use matlab inv-function

siz = size(x);
if all(siz(1:2)==2),
  adjx  = [x(2,2,:,:) -x(1,2,:,:); -x(2,1,:,:) x(1,1,:,:)];
  denom = det2x2(x);
  d     = adjx./denom([1 1],[1 1],:,:);
elseif all(siz(1:2)==3),
  adjx = [ det2x2(x([2 3],[2 3],:,:)) -det2x2(x([1 3],[2 3],:,:))  det2x2(x([1 2],[2 3],:,:)); ...
          -det2x2(x([2 3],[1 3],:,:))  det2x2(x([1 3],[1 3],:,:)) -det2x2(x([1 2],[1 3],:,:)); ...
	   det2x2(x([2 3],[1 2],:,:)) -det2x2(x([1 3],[1 2],:,:))  det2x2(x([1 2],[1 2],:,:))];
  denom = det2x2(x);
  d     = adjx./denom([1 1 1],[1 1 1],:,:);
elseif numel(siz)==2,
  d = inv(x);
else
  error('cannot compute slicewise inverse');
  %write for loop for the higher dimensions, using normal inv
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix multiplication 2x2
function [z] = mtimes2x2(x, y)

% compute x*y 
% and dimensionatity is 2x2xN

z     = complex(zeros(size(x)));
xconj = conj(x);

z(1,1,:,:) = x(1,1,:,:).*y(1,1,:,:) + x(1,2,:,:).*y(2,1,:,:);
z(1,2,:,:) = x(1,1,:,:).*y(1,2,:,:) + x(1,2,:,:).*y(2,2,:,:);
z(2,1,:,:) = x(2,1,:,:).*y(1,1,:,:) + x(2,2,:,:).*y(2,1,:,:);
z(2,2,:,:) = x(2,1,:,:).*y(1,2,:,:) + x(2,2,:,:).*y(2,2,:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quadratic multiplication 2x2 matrix,  sandwiched matrix assumed hermitian
function [z] = sandwich2x2(x, y)

% compute x*y*x' provided y = hermitian
% and dimensionatity is 2x2xN

%FIXME build in check for hermitianity
z     = complex(zeros(size(x)));
xconj = conj(x);
xabs2 = abs(x).^2;

z(1,1,:,:) = xabs2(1,1,:,:) .* y(1,1,:,:) + ...
             xabs2(1,2,:,:) .* y(2,2,:,:) + ...
           2.*real(y(2,1,:,:).*xconj(1,1,:,:).*x(1,2,:,:));
z(2,1,:,:) = y(1,1,:,:).*xconj(1,1,:,:).*x(2,1,:,:) + ...
           y(2,1,:,:).*xconj(1,1,:,:).*x(2,2,:,:) + ...
           y(1,2,:,:).*xconj(1,2,:,:).*x(2,1,:,:) + ...
           y(2,2,:,:).*xconj(1,2,:,:).*x(2,2,:,:);
z(1,2,:,:) = conj(z(2,1,:,:));
z(2,2,:,:) = xabs2(2,1,:,:) .* y(1,1,:,:) + ...
             xabs2(2,2,:,:) .* y(2,2,:,:) + ...
           2.*real(y(1,2,:,:).*xconj(2,2,:,:).*x(2,1,:,:));

%b1 b2     a1 a2'   b1' b3'
%b3 b4     a2 a3    b2' b4'
%
%b1*a1+b2*a2  b1*a2'+b2*a3  b1' b3'
%b3*a1+b4*a2  b3*a2'+b4*a3  b2' b4'
%
%b1*a1*b1'+b2*a2*b1'+b1*a2'*b2'+b2*a3*b2' b1*a1*b3'+b2*a2*b3'+b1*a2'*b4'+b2*a3*b4'
%b3*a1*b1'+b4*a2*b1'+b3*a2'*b2'+b4*a3*b2' b3*a1*b3'+b4*a2*b3'+b3*a2'*b4'+b4*a3*b4'
%
%a1*abs(b1)^2 + a2*(b1'*b2) + a2'*(b1*b2') + a3*abs(b2)^2    a1*b1*b3'    + a2*b2*b3'   + a2'*b1*b4'   + a3*b2*b4'
%a1*b1'*b3    + a2*b1'*b4   + a2'*b2'*b3   + a3*b2'*b4       a1*abs(b3)^2 + a2*(b3'*b4) + a2'*(b3*b4') + a3*abs(b4)^2
