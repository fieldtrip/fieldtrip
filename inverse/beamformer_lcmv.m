function [dipout] = beamformer_lcmv(dip, grad, vol, dat, Cy, varargin)

% BEAMFORMER_LCMV scans on pre-defined dipole locations with a single dipole
% and returns the beamformer spatial filter output for a dipole on every
% location.  Dipole locations that are outside the head will return a
% NaN value.
%
% Use as
%   [dipout] = beamformer_lcmv(dipin, grad, vol, dat, cov, varargin)
% where
%   dipin       is the input dipole model
%   grad        is the gradiometer definition
%   vol         is the volume conductor definition
%   dat         is the data matrix with the ERP or ERF
%   cov         is the data covariance or cross-spectral density matrix
% and
%   dipout      is the resulting dipole model with all details
%
% The input dipole model consists of
%   dipin.pos   positions for dipole, e.g. regular grid, Npositions x 3
%   dipin.mom   dipole orientation (optional), 3 x Npositions
%
% Additional options should be specified in key-value pairs and can be
%  'lambda'           = regularisation parameter
%  'powmethod'        = can be 'trace' or 'lambda1'
%  'feedback'         = give ft_progress indication, can be 'text', 'gui' or 'none' (default)
%  'fixedori'         = use fixed or free orientation,                 can be 'yes' or 'no'
%  'projectnoise'     = project noise estimate through filter,         can be 'yes' or 'no'
%  'projectmom'       = project the dipole moment timecourse on the direction of maximal power, can be 'yes' or 'no'
%  'keepfilter'       = remember the beamformer filter,                can be 'yes' or 'no'
%  'keepleadfield'    = remember the forward computation,              can be 'yes' or 'no'
%  'keepmom'          = remember the estimated dipole moment,          can be 'yes' or 'no'
%  'keepcov'          = remember the estimated dipole covariance,      can be 'yes' or 'no'
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

% Copyright (C) 2003-2008, Robert Oostenveld
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
% $Id$

if mod(nargin-5,2)
  % the first 5 arguments are fixed, the other arguments should come in pairs
  error('invalid number of optional arguments');
end

% these optional settings do not have defaults
powmethod      = keyval('powmethod',     varargin); % the default for this is set below
subspace       = keyval('subspace',      varargin); % used to implement an "eigenspace beamformer" as described in Sekihara et al. 2002 in HBM
% these settings pertain to the forward model, the defaults are set in compute_leadfield
reducerank     = keyval('reducerank',     varargin);
normalize      = keyval('normalize',      varargin);
normalizeparam = keyval('normalizeparam', varargin);
% these optional settings have defaults
feedback       = keyval('feedback',      varargin); if isempty(feedback),      feedback = 'text';            end
keepfilter     = keyval('keepfilter',    varargin); if isempty(keepfilter),    keepfilter = 'no';            end
keepleadfield  = keyval('keepleadfield', varargin); if isempty(keepleadfield), keepleadfield = 'no';         end
keepcov        = keyval('keepcov',       varargin); if isempty(keepcov),       keepcov = 'no';               end
keepmom        = keyval('keepmom',       varargin); if isempty(keepmom),       keepmom = 'yes';              end
lambda         = keyval('lambda',        varargin); if isempty(lambda  ),      lambda = 0;                   end
projectnoise   = keyval('projectnoise',  varargin); if isempty(projectnoise),  projectnoise = 'yes';         end
projectmom     = keyval('projectmom',    varargin); if isempty(projectmom),    projectmom = 'no';            end
fixedori       = keyval('fixedori',      varargin); if isempty(fixedori),      fixedori = 'no';              end

% convert the yes/no arguments to the corresponding logical values
keepfilter     = strcmp(keepfilter,    'yes');
keepleadfield  = strcmp(keepleadfield, 'yes');
keepcov        = strcmp(keepcov,       'yes');
keepmom        = strcmp(keepmom,       'yes');
projectnoise   = strcmp(projectnoise,  'yes');
projectmom     = strcmp(projectmom,    'yes');
fixedori       = strcmp(fixedori,      'yes');

% default is to use the trace of the covariance matrix, see Van Veen 1997
if isempty(powmethod)
  powmethod = 'trace';
end

% use these two logical flags instead of doing the string comparisons each time again
powtrace   = strcmp(powmethod, 'trace');
powlambda1 = strcmp(powmethod, 'lambda1');

if isfield(dip, 'mom') && fixedori
  error('you cannot specify a dipole orientation and fixedmom simultaneously');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the dipole positions that are inside/outside the brain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(dip, 'inside') && ~isfield(dip, 'outside');
  insideLogical = ft_inside_vol(dip.pos, vol);
  dip.inside = find(insideLogical);
  dip.outside = find(~dip.inside);
elseif isfield(dip, 'inside') && ~isfield(dip, 'outside');
  dip.outside    = setdiff(1:size(dip.pos,1), dip.inside);
elseif ~isfield(dip, 'inside') && isfield(dip, 'outside');
  dip.inside     = setdiff(1:size(dip.pos,1), dip.outside);
end

% select only the dipole positions inside the brain for scanning
dip.origpos     = dip.pos;
dip.originside  = dip.inside;
dip.origoutside = dip.outside;
if isfield(dip, 'mom')
  dip.mom = dip.mom(:, dip.inside);
end
if isfield(dip, 'leadfield')
  fprintf('using precomputed leadfields\n');
  dip.leadfield = dip.leadfield(dip.inside);
end
if isfield(dip, 'filter')
  fprintf('using precomputed filters\n');
  dip.filter = dip.filter(dip.inside);
end
if isfield(dip, 'subspace')
  fprintf('using subspace projection\n');
  dip.subspace = dip.subspace(dip.inside);
end
dip.pos     = dip.pos(dip.inside, :);
dip.inside  = 1:size(dip.pos,1);
dip.outside = [];

isrankdeficient = (rank(Cy)<size(Cy,1));

% it is difficult to give a quantitative estimate of lambda, therefore also
% support relative (percentage) measure that can be specified as string (e.g. '10%')
if ~isempty(lambda) && ischar(lambda) && lambda(end)=='%'
  ratio = sscanf(lambda, '%f%%');
  ratio = ratio/100;
  lambda = ratio * trace(Cy)/size(Cy,1);
end

if projectnoise
  % estimate the noise power, which is further assumed to be equal and uncorrelated over channels
  if isrankdeficient
    % estimated noise floor is equal to or higher than lambda
    noise = lambda;
  else
    % estimate the noise level in the covariance matrix by the smallest singular value
    noise = svd(Cy);
    noise = noise(end);
    % estimated noise floor is equal to or higher than lambda
    noise = max(noise, lambda);
  end
end

% the inverse only has to be computed once for all dipoles
invCy = pinv(Cy + lambda * eye(size(Cy)));
if isfield(dip, 'subspace')
  fprintf('using source-specific subspace projection\n');
  % remember the original data prior to the voxel dependent subspace projection
  dat_pre_subspace = dat;
  Cy_pre_subspace  = Cy;
elseif ~isempty(subspace)
  fprintf('using data-specific subspace projection\n');
  % TODO implement an "eigenspace beamformer" as described in Sekihara et al. 2002 in HBM
  if numel(subspace)==1,
    % interpret this as a truncation of the eigenvalue-spectrum
    % if <1 it is a fraction of the largest eigenvalue
    % if >=1 it is the number of largest eigenvalues
    dat_pre_subspace = dat;
    Cy_pre_subspace  = Cy;
    [u, s, v] = svd(real(Cy));
    if subspace<1,
      sel      = find(diag(s)./s(1,1) > subspace);
      subspace = max(sel);
    else
      Cy       = s(1:subspace,1:subspace);
      % this is equivalent to subspace*Cy*subspace' but behaves well numerically
      % by construction.
      invCy    = diag(1./diag(Cy));
      subspace = u(:,1:subspace)';
      dat      = subspace*dat;
    end
  else
    dat_pre_subspace = dat;
    Cy_pre_subspace  = Cy;
    Cy    = subspace*Cy*subspace'; % here the subspace can be different from
    % the singular vectors of Cy, so we have to do the sandwiching as opposed
    % to line 216
    invCy = pinv(Cy);
    dat   = subspace*dat;
  end
end

% start the scanning with the proper metric
ft_progress('init', feedback, 'scanning grid');

for i=1:size(dip.pos,1)
  if isfield(dip, 'leadfield') && isfield(dip, 'mom') && size(dip.mom, 1)==size(dip.leadfield{i}, 2)
    % reuse the leadfield that was previously computed and project
    lf = dip.leadfield{i} * dip.mom(:,i);
  elseif  isfield(dip, 'leadfield') &&  isfield(dip, 'mom')
    % reuse the leadfield that was previously computed but don't project
    lf = dip.leadfield{i};
  elseif  isfield(dip, 'leadfield') && ~isfield(dip, 'mom')
    % reuse the leadfield that was previously computed
    lf = dip.leadfield{i};    
  elseif  ~isfield(dip, 'leadfield') && isfield(dip, 'mom')
    % compute the leadfield for a fixed dipole orientation
    lf = ft_compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam) * dip.mom(:,i);
  else
    % compute the leadfield
    lf = ft_compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam);
  end
  
  if isfield(dip, 'subspace')
    % do subspace projection of the forward model
    lf = dip.subspace{i} * lf;
    % the data and the covariance become voxel dependent due to the projection
    dat   = dip.subspace{i} * dat_pre_subspace;
    Cy    = dip.subspace{i} * (Cy_pre_subspace + lambda * eye(size(Cy_pre_subspace))) * dip.subspace{i}';
    invCy = pinv(dip.subspace{i} * (Cy_pre_subspace + lambda * eye(size(Cy_pre_subspace))) * dip.subspace{i}');
  elseif ~isempty(subspace)
    % do subspace projection of the forward model only
    lforig = lf;
    lf     = subspace * lf;
    
    % according to Kensuke's paper, the eigenspace bf boils down to projecting
    % the 'traditional' filter onto the subspace
    % spanned by the first k eigenvectors [u,s,v] = svd(Cy); filt = ESES*filt; 
    % ESES = u(:,1:k)*u(:,1:k)';
    % however, even though it seems that the shape of the filter is identical to
    % the shape it is obtained with the following code, the w*lf=I does not hold.
  end
  
  if fixedori
    % compute the leadfield for the optimal dipole orientation
    % subsequently the leadfield for only that dipole orientation will be used for the final filter computation
    % filt = pinv(lf' * invCy * lf) * lf' * invCy;
    % [u, s, v] = svd(real(filt * Cy * ctranspose(filt)));
    % in this step the filter computation is not necessary, use the quick way to compute the voxel level covariance (cf. van Veen 1997)
    [u, s, v] = svd(real(pinv(lf' * invCy *lf)));
    eta = u(:,1);
    lf  = lf * eta;
    if ~isempty(subspace), lforig = lforig * eta; end
    dipout.ori{i} = eta;
  end
  
  if isfield(dip, 'filter')
    % use the provided filter
    filt = dip.filter{i};
  else
    % construct the spatial filter
    filt = pinv(lf' * invCy * lf) * lf' * invCy;              % van Veen eqn. 23, use PINV/SVD to cover rank deficient leadfield
  end
  if projectmom
    [u, s, v] = svd(filt * Cy * ctranspose(filt));
    mom = u(:,1);
    filt = (mom') * filt;
  end
  if powlambda1
    % dipout.pow(i) = lambda1(pinv(lf' * invCy * lf));        % this is more efficient if the filters are not present
    dipout.pow(i) = lambda1(filt * Cy * ctranspose(filt));    % this is more efficient if the filters are present
  elseif powtrace
    % dipout.pow(i) = trace(pinv(lf' * invCy * lf));          % this is more efficient if the filters are not present, van Veen eqn. 24
    dipout.pow(i) = trace(filt * Cy * ctranspose(filt));      % this is more efficient if the filters are present
  end
  if keepcov
    % compute the source covariance matrix
    dipout.cov{i} = filt * Cy * ctranspose(filt);
  end
  if keepmom && ~isempty(dat)
    % estimate the instantaneous dipole moment at the current position
    dipout.mom{i} = filt * dat;
  end
  if projectnoise
    % estimate the power of the noise that is projected through the filter
    if powlambda1
      dipout.noise(i) = noise * lambda1(filt * ctranspose(filt));
    elseif powtrace
      dipout.noise(i) = noise * trace(filt * ctranspose(filt));
    end
    if keepcov
      dipout.noisecov{i} = noise * filt * ctranspose(filt);
    end
  end
  if keepfilter
    if ~isempty(subspace)
      dipout.filter{i} = filt*subspace;
      %dipout.filter{i} = filt*pinv(subspace);
    else
      dipout.filter{i} = filt;
    end
  end
  if keepleadfield
    if ~isempty(subspace)
      dipout.leadfield{i} = lforig;
    else
      dipout.leadfield{i} = lf;
    end
  end
  ft_progress(i/size(dip.pos,1), 'scanning grid %d/%d\n', i, size(dip.pos,1));
end

ft_progress('close');

dipout.inside  = dip.originside;
dipout.outside = dip.origoutside;
dipout.pos     = dip.origpos;

% reassign the scan values over the inside and outside grid positions
if isfield(dipout, 'leadfield')
  dipout.leadfield(dipout.inside)  = dipout.leadfield;
  dipout.leadfield(dipout.outside) = {[]};
end
if isfield(dipout, 'filter')
  dipout.filter(dipout.inside)  = dipout.filter;
  dipout.filter(dipout.outside) = {[]};
end
if isfield(dipout, 'mom')
  dipout.mom(dipout.inside)  = dipout.mom;
  dipout.mom(dipout.outside) = {[]};
end
if isfield(dipout, 'ori')
  dipout.ori(dipout.inside)  = dipout.ori;
  dipout.ori(dipout.outside) = {[]};
end
if isfield(dipout, 'cov')
  dipout.cov(dipout.inside)  = dipout.cov;
  dipout.cov(dipout.outside) = {[]};
end
if isfield(dipout, 'noisecov')
  dipout.noisecov(dipout.inside)  = dipout.noisecov;
  dipout.noisecov(dipout.outside) = {[]};
end
if isfield(dipout, 'pow')
  dipout.pow(dipout.inside)  = dipout.pow;
  dipout.pow(dipout.outside) = nan;
end
if isfield(dipout, 'noise')
  dipout.noise(dipout.inside)  = dipout.noise;
  dipout.noise(dipout.outside) = nan;
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
% standard Matlab function, except that the default tolerance is twice as
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

