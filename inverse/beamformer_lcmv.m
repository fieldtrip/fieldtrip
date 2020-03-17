function [dipout] = beamformer_lcmv(dip, grad, headmodel, dat, C, varargin)

% BEAMFORMER_LCMV scans on pre-defined dipole locations with a single dipole
% and returns the beamformer spatial filter output for a dipole on every
% location. Dipole locations that are outside the head will return a
% NaN value.
%
% Use as
%   [dipout] = beamformer_lcmv(dipin, grad, headmodel, dat, cov, varargin)
% where
%   dipin       is the input dipole model
%   grad        is the gradiometer definition
%   headmodel   is the volume conductor definition
%   dat         is the data matrix with the ERP or ERF
%   cov         is the data covariance or cross-spectral density matrix
% and
%   dipout      is the resulting dipole model with all details
%
% The input dipole model consists of
%   dipin.pos   positions for dipole, e.g. regular grid, Npositions x 3
%   dipin.mom   dipole orientation (optional), 3 x Npositions,
% and can additionally contain things like a precomputed filter.
%
% Additional options should be specified in key-value pairs and can be
%  'lambda'           = regularisation parameter
%  'kappa'            = parameter for covariance matrix inversion
%  'tol'              = parameter for covariance matrix inversion
%  'powmethod'        = can be 'trace' or 'lambda1'
%  'feedback'         = give ft_progress indication, can be 'text', 'gui' or 'none' (default)
%  'fixedori'         = use fixed or free orientation,                   can be 'yes' or 'no'
%  'projectnoise'     = project noise estimate through filter,           can be 'yes' or 'no'
%  'projectmom'       = project the dipole moment timecourse on the direction of maximal power, can be 'yes' or 'no'
%  'keepfilter'       = remember the beamformer filter,                  can be 'yes' or 'no'
%  'keepleadfield'    = remember the forward computation,                can be 'yes' or 'no'
%  'keepmom'          = remember the estimated dipole moment timeseries, can be 'yes' or 'no'
%  'keepcov'          = remember the estimated dipole covariance,        can be 'yes' or 'no'
%  'kurtosis'         = compute the kurtosis of the dipole timeseries,   can be 'yes' or 'no'
%  'weightnorm'       = normalize the beamformer weights,                can be 'no', 'unitnoisegain' or 'nai'
%
% These options influence the forward computation of the leadfield
%  'reducerank'       = reduce the leadfield rank, can be 'no' or a number (e.g. 2)
%  'normalize'        = normalize the leadfield
%  'normalizeparam'   = parameter for depth normalization (default = 0.5)
%
% If the dipole definition only specifies the dipole location, a rotating
% dipole (regional source) is assumed on each location. If a dipole moment
% is specified, its orientation will be used and only the strength will
% be fitted to the data.

% Copyright (C) 2003-2020, Robert Oostenveld and Jan-Mathijs Schoffelen
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

if mod(nargin-5,2)
  % the first 5 arguments are fixed, the other arguments should come in pairs
  ft_error('invalid number of optional arguments');
end

% these optional settings do not have defaults
powmethod      = ft_getopt(varargin, 'powmethod'); % the default for this is set below
subspace       = ft_getopt(varargin, 'subspace'); % used to implement an "eigenspace beamformer" as described in Sekihara et al. 2002 in HBM
% these settings pertain to the forward model, the defaults are set in compute_leadfield
reducerank     = ft_getopt(varargin, 'reducerank');
normalize      = ft_getopt(varargin, 'normalize');
normalizeparam = ft_getopt(varargin, 'normalizeparam');
% these optional settings have defaults
feedback       = ft_getopt(varargin, 'feedback',      'text');
keepfilter     = ft_getopt(varargin, 'keepfilter',    'no');
keepleadfield  = ft_getopt(varargin, 'keepleadfield', 'no');
keepcov        = ft_getopt(varargin, 'keepcov',       'no');
keepmom        = ft_getopt(varargin, 'keepmom',       'yes');
lambda         = ft_getopt(varargin, 'lambda',        0);
kappa          = ft_getopt(varargin, 'kappa',         []);
tol            = ft_getopt(varargin, 'tol',           []);
invmethod      = ft_getopt(varargin, 'invmethod',     []);
projectnoise   = ft_getopt(varargin, 'projectnoise',  'yes');
projectmom     = ft_getopt(varargin, 'projectmom',    'no');
fixedori       = ft_getopt(varargin, 'fixedori',      'no');
computekurt    = ft_getopt(varargin, 'kurtosis',      'no');
weightnorm     = ft_getopt(varargin, 'weightnorm',    'no');

% convert the yes/no arguments to the corresponding logical values
keepfilter     = istrue(keepfilter);
keepleadfield  = istrue(keepleadfield);
keepcov        = istrue(keepcov);
keepmom        = istrue(keepmom);
projectnoise   = istrue(projectnoise);
projectmom     = istrue(projectmom);
fixedori       = istrue(fixedori);
computekurt    = istrue(computekurt);
dofeedback     = ~strcmp(feedback, 'none');

% default is to use the trace of the covariance matrix, see Van Veen 1997
if isempty(powmethod)
  powmethod = 'trace';
end

% use these two logical flags instead of doing the string comparisons each time again
powtrace   = strcmp(powmethod, 'trace');
powlambda1 = strcmp(powmethod, 'lambda1');

if isfield(dip, 'mom') && fixedori
  ft_error('you cannot specify a dipole orientation and fixedmom simultaneously');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the dipole positions that are inside/outside the brain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(dip, 'inside')
  dip.inside = ft_inside_headmodel(dip.pos, headmodel);
end

if any(dip.inside>1)
  % convert to logical representation
  tmp = false(size(dip.pos,1),1);
  tmp(dip.inside) = true;
  dip.inside = tmp;
end

% keep the original details on inside and outside positions
originside = dip.inside;
origpos    = dip.pos;

% flags to avoid calling isfield repeatedly in the loop over grid positions (saves a lot of time)
hasmom       = false;
hasleadfield = false;
hasfilter    = false;
hassubspace  = false;

% select only the dipole positions inside the brain for scanning
dip.pos    = dip.pos(originside,:);
dip.inside = true(size(dip.pos,1),1);
if isfield(dip, 'mom')
  hasmom  = 1;
  dip.mom = dip.mom(:, originside);
end
if isfield(dip, 'leadfield')
  hasleadfield = 1;
  if dofeeback
    fprintf('using precomputed leadfields\n');
  end
  dip.leadfield = dip.leadfield(originside);
end
if isfield(dip, 'filter')
  hasfilter = 1;
  if dofeedback
    fprintf('using precomputed filters\n');
  end
  dip.filter = dip.filter(originside);
end
if isfield(dip, 'subspace')
  hassubspace = 1;
  if dofeedback
    fprintf('using subspace projection\n');
  end
  dip.subspace = dip.subspace(originside);
end

rankC = rank(C);

% it is difficult to give a quantitative estimate of lambda, therefore also
% support relative (percentage) measure that can be specified as string (e.g. '10%')
if ~isempty(lambda) && ischar(lambda) && lambda(end)=='%'
  ratio = sscanf(lambda, '%f%%');
  ratio = ratio/100;
  tmplambda = ratio * trace(C)/size(C,1);
elseif ~isempty(lambda)
  tmplambda = lambda;
else
  tmplambda = 0;
end

if projectnoise || strcmp(weightnorm, 'nai')
    % estimate the noise level in the covariance matrix by the smallest singular (non-zero) value
    % always needed for the NAI weight normalization case
    noise = svd(C);
    noise = noise(rankC);
    % estimated noise floor is equal to or higher than lambda
    noise = max(noise, tmplambda);
end

% the inverse of the covariance matrix only has to be computed once for all dipoles
if hassubspace
  if dofeedback
    fprintf('using source-specific subspace projection\n');
  end
  % remember the original data prior to the voxel dependent subspace projection
  dat_pre_subspace = dat;
  C_pre_subspace  = C;
elseif ~isempty(subspace)
  if dofeedback
    fprintf('using data-specific subspace projection\n');
  end
  
  % TODO implement an "eigenspace beamformer" as described in Sekihara et al. 2002 in HBM
  if numel(subspace)==1
    % interpret this as a truncation of the eigenvalue-spectrum
    % if <1 it is a fraction of the largest eigenvalue
    % if >=1 it is the number of largest eigenvalues
    dat_pre_subspace = dat;
    C_pre_subspace  = C;
    [u, s, v] = svd(real(C));
    if subspace<1
      subspace = find(diag(s)./s(1,1) > subspace, 1, 'last');
    end
    
    C       = s(1:subspace,1:subspace);
    % this is equivalent to subspace*C*subspace' but behaves well numerically by construction.
    invC     = diag(1./diag(C + lambda * eye(size(C))));
    subspace = u(:,1:subspace)';
    dat      = subspace*dat;
  else
    dat_pre_subspace = dat;
    C_pre_subspace   = C;
    C                = subspace*C*subspace'; 
    % here the subspace can be different from the singular vectors of C, so we
    % have to do the sandwiching as opposed to line 223
    invC = ft_inv(C, 'lambda', lambda, 'kappa', kappa, 'tolerance', tol, 'method', invmethod);
    dat   = subspace*dat;
  end
else
  invC = ft_inv(C, 'lambda', lambda, 'kappa', kappa, 'tolerance', tol, 'method', invmethod);
end

% compute the square of invC, which might be needed for unitnoisegain or NAI constraint
invC_squared = invC^2;

% start the scanning with the proper metric
ft_progress('init', feedback, 'scanning grid');
for i=1:size(dip.pos,1)
  if hasleadfield && hasmom && size(dip.mom, 1)==size(dip.leadfield{i}, 2)
    % reuse the leadfield that was previously computed and project
    lf = dip.leadfield{i} * dip.mom(:,i);
  elseif  hasleadfield &&  hasmom
    % reuse the leadfield that was previously computed but don't project
    lf = dip.leadfield{i};
  elseif  hasleadfield && ~hasmom
    % reuse the leadfield that was previously computed
    lf = dip.leadfield{i};    
  elseif  ~hasleadfield && hasmom
    % compute the leadfield for a fixed dipole orientation
    lf = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam) * dip.mom(:,i);
  else
    % compute the leadfield
    lf = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam);
  end
  
  if hassubspace
    % do subspace projection of the forward model
    lf    = dip.subspace{i} * lf;
    % the data and the covariance become voxel dependent due to the projection
    dat   = dip.subspace{i} * dat_pre_subspace;
    C    = dip.subspace{i} * C_pre_subspace * dip.subspace{i}';
    invC = ft_inv(dip.subspace{i} * C_pre_subspace * dip.subspace{i}', 'lambda', lambda, 'kappa', kappa, 'tolerance', tol, 'method', invmethod);
  elseif ~isempty(subspace)
    % do subspace projection of the forward model only
    lforig = lf;
    lf     = subspace * lf;
    
    % according to Kensuke's paper, the eigenspace bf boils down to projecting
    % the 'traditional' filter onto the subspace spanned by the first k eigenvectors 
    % [u,s,v] = svd(C); filt = ESES*filt; ESES = u(:,1:k)*u(:,1:k)';
    % however, even though it seems that the shape of the filter is identical to
    % the shape it is obtained with the following code, the w*lf=I does not hold.
  end

  if fixedori
    switch(weightnorm)
      case {'unitnoisegain','nai'}
        % optimal orientation calculation for unit-noise gain beamformer,
        % (also applies to similar NAI), based on equation 4.47 from Sekihara & Nagarajan (2008)
        [vv, dd]     = eig(pinv(lf' * invC_squared *lf)*(lf' * invC *lf));
        [dum,maxeig] = max(diag(dd));
        unitnoiseori = vv(:,maxeig);
        lf  = lf * unitnoiseori;
        if hassubspace, lforig = lforig * unitnoiseori; end
        dipout.ori{i} = unitnoiseori;
      otherwise
        % compute the leadfield for the optimal dipole orientation
        % subsequently the leadfield for only that dipole orientation will be used for the final filter computation
        % filt = pinv(lf' * invC * lf) * lf' * invC;
        % [u, s, v] = svd(real(filt * C * ctranspose(filt)));
        % in this step the filter computation is not necessary, use the quick way to compute the voxel level covariance (cf. van Veen 1997)
        [u, s, v] = svd(real(pinv(lf' * invC *lf)));
        maxpowori = u(:,1);
        eta = s(1,1)./s(2,2); % ratio between the first and second singular valu

        % and compute the leadfield for that orientation
        lf  = lf * maxpowori;
        dipout.ori{i} = maxpowori;
        dipout.eta{i} = eta;
        if ~isempty(subspace), lforig = lforig * maxpowori; end
    end
  end
  
  if hasfilter
    % use the provided filter
    filt = dip.filter{i};
  elseif strcmp(weightnorm,'nai')
    % Van Veen's Neural Activity Index
    % below equation is equivalent to following:  
    % filt = pinv(lf' * invC * lf) * lf' * invC; 
    % filt = filt/sqrt(noise*filt*filt');
    filt = pinv(sqrt(noise * lf' * invC_squared * lf)) * lf' *invC; % based on Sekihara & Nagarajan 2008 eqn. 4.15
  elseif strcmp(weightnorm,'unitnoisegain')
    % Unit-noise gain minimum variance (aka Borgiotti-Kaplan) beamformer
    % below equation is equivalent to following:  
    % filt = pinv(lf' * invC * lf) * lf' * invC; 
    % filt = filt/sqrt(filt*filt');
    filt = pinv(sqrtm(lf' * invC_squared * lf)) * lf' *invC;     % Sekihara & Nagarajan 2008 eqn. 4.15
  else
    % construct the spatial filter
    filt = pinv(lf' * invC * lf) * lf' * invC;              % van Veen eqn. 23, use PINV/SVD to cover rank deficient leadfield
  end

  cfilt = filt * C * ctranpose(filt);
  if projectmom
    [u, s, v] = svd(cfilt);
    mom = u(:,1); % dominant dipole direction
    filt = (mom') * filt;
  end
  if powlambda1
    % dipout.pow(i) = lambda1(pinv(lf' * invC * lf)); % this is more efficient if the filters are not present
    dipout.pow(i,1) = lambda1(cfilt);                 % this is more efficient if the filters are present
  elseif powtrace
    % dipout.pow(i) = trace(pinv(lf' * invC * lf));   % this is more efficient if the filters are not present, van Veen eqn. 24
    dipout.pow(i,1) = trace(cfilt);                   % this is more efficient if the filters are present
  end
  if keepcov
    % compute the source covariance matrix
    dipout.cov{i,1} = cfilt;
  end
  if keepmom && ~isempty(dat)
    % estimate the instantaneous dipole moment at the current position
    dipout.mom{i,1} = filt * dat;
  end
  if computekurt && ~isempty(dat)
    % compute the kurtosis of the dipole time series
    dipout.kurtosis(i,:) = kurtosis((filt*dat)');
  end    
  if projectnoise
    % estimate the power of the noise that is projected through the filter
    if powlambda1
      dipout.noise(i,1) = noise * lambda1(filt * ctranspose(filt));
    elseif powtrace
      dipout.noise(i,1) = noise * trace(filt * ctranspose(filt));
    end
    if keepcov
      dipout.noisecov{i,1} = noise * filt * ctranspose(filt);
    end
  end
  if keepfilter
    if ~isempty(subspace)
      dipout.filter{i,1} = filt*subspace;
      %dipout.filter{i} = filt*pinv(subspace);
    else
      dipout.filter{i,1} = filt;
    end
  end
  if keepleadfield
    if ~isempty(subspace)
      dipout.leadfield{i,1} = lforig;
    else
      dipout.leadfield{i,1} = lf;
    end
  end
  ft_progress(i/size(dip.pos,1), 'scanning grid %d/%d\n', i, size(dip.pos,1));
end

ft_progress('close');

% reassign the scan values over the inside and outside grid positions
dipout.pos     = origpos;
dipout.inside  = originside;
if isfield(dipout, 'leadfield')
  dipout.leadfield( originside) = dipout.leadfield;
  dipout.leadfield(~originside) = {[]};
end
if isfield(dipout, 'filter')
  dipout.filter( originside) = dipout.filter;
  dipout.filter(~originside) = {[]};
end
if isfield(dipout, 'mom')
  dipout.mom( originside) = dipout.mom;
  dipout.mom(~originside) = {[]};
end
if isfield(dipout, 'ori')
  dipout.ori( originside) = dipout.ori;
  dipout.ori(~originside) = {[]};
end
if isfield(dipout, 'cov')
  dipout.cov( originside) = dipout.cov;
  dipout.cov(~originside) = {[]};
end
if isfield(dipout, 'noisecov')
  dipout.noisecov( originside) = dipout.noisecov;
  dipout.noisecov(~originside) = {[]};
end
if isfield(dipout, 'pow')
  dipout.pow( originside) = dipout.pow;
  dipout.pow(~originside) = nan;
end
if isfield(dipout, 'noise')
  dipout.noise( originside) = dipout.noise;
  dipout.noise(~originside) = nan;
end
if isfield(dipout, 'kurtosis')
  dipout.kurtosis( originside) = dipout.kurtosis;
  dipout.kurtosis(~originside) = nan;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to obtain the largest singular value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = lambda1(x)
% determine the largest singular value, which corresponds to the power along the dominant direction
s = svd(x);
s = s(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to compute the pseudo inverse. This is the same as the
% standard MATLAB function, except that the default tolerance is twice as
% high.
%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision$  $Date: 2009/03/23 21:14:42 $
%   default tolerance increased by factor 2 (Robert Oostenveld, 7 Feb 2004)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = pinv(A,varargin)
[m,n] = size(A);
if n > m
  X = pinv(A',varargin{:})';
else
  [U,S,V] = svd(A,0);
  if m > 1, s = diag(S);
  elseif m == 1, s = S(1);
  else s = 0;
  end
  if nargin == 2
    tol = varargin{1};
  else
    tol = 10 * max(m,n) * max(s) * eps;
  end
  r = sum(s > tol);
  if (r == 0)
    X = zeros(size(A'),class(A));
  else
    s = diag(ones(r,1)./s(1:r));
    X = V(:,1:r)*s*U(:,1:r)';
  end
end

