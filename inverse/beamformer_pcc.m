function [dipout] = beamformer_pcc(dip, grad, headmodel, dat, C, varargin)

% BEAMFORMER_PCC implements a linearly-constrained miminum variance  beamformer
% that allows for post-hoc computation of canonical or partial coherence or
% correlation. Moreover, if cortico-cortical interactions are computed, the
% spatial filters are computed with a paired dipole as sourcemodel, thus
% suppressing the distortive effect of correlated activity from the seed region.
% Dipole locations that are outside the head will return a NaN value.
%
% Use as
%   [dipout] = beamformer_pcc(dipin, grad, headmodel, dat, cov, ...)
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
%   dipin.mom   dipole orientation (optional), 3 x Npositions
% and can additionally contain things like a precomputed filter.
%
% Additional options should be specified in key-value pairs and can be
%   refchan
%   refdip
%   supchan
%   supdip
%   reducerank
%   normalize
%   normalizeparam
%   feedback
%   keepcsd
%   keepfilter
%   keepleadfield
%   keepmom
%   lambda
%   projectnoise
%   realfilter
%   fixedori

% Copyright (C) 2005-2020, Robert Oostenveld & Jan-Mathijs Schoffelen
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
refchan        = ft_getopt(varargin, 'refchan',       []);
refdip         = ft_getopt(varargin, 'refdip',        []);
supchan        = ft_getopt(varargin, 'supchan',       []);
supdip         = ft_getopt(varargin, 'supdip',        []);

% these optional settings have defaults
feedback       = ft_getopt(varargin, 'feedback',      'text');
keepcsd        = ft_getopt(varargin, 'keepcsd',       'yes');
keepfilter     = ft_getopt(varargin, 'keepfilter',    'no');
keepleadfield  = ft_getopt(varargin, 'keepleadfield', 'no');
keepmom        = ft_getopt(varargin, 'keepmom',       'yes');
projectnoise   = ft_getopt(varargin, 'projectnoise',  'yes');
realfilter     = ft_getopt(varargin, 'realfilter',    'yes');
fixedori       = ft_getopt(varargin, 'fixedori',      'no');
weightnorm     = ft_getopt(varargin, 'weightnorm',    'no');

% construct the low-level options for the covariance matrix inversion as key-value pairs, these are passed to FT_INV
invopt = {};
invopt = ft_setopt(invopt, 'lambda',    ft_getopt(varargin, 'lambda', 0));
invopt = ft_setopt(invopt, 'kappa',     ft_getopt(varargin, 'kappa'));
invopt = ft_setopt(invopt, 'tolerance', ft_getopt(varargin, 'tol'));
invopt = ft_setopt(invopt, 'method',    ft_getopt(varargin, 'invmethod'));

% construct the low-level options for the leadfield computation as key-value pairs, these are passed to FT_COMPUTE_LEADFIELD
leadfieldopt = {};
leadfieldopt = ft_setopt(leadfieldopt, 'reducerank',     ft_getopt(varargin, 'reducerank'));
leadfieldopt = ft_setopt(leadfieldopt, 'backproject',    ft_getopt(varargin, 'backproject'));
leadfieldopt = ft_setopt(leadfieldopt, 'normalize',      ft_getopt(varargin, 'normalize'));
leadfieldopt = ft_setopt(leadfieldopt, 'normalizeparam', ft_getopt(varargin, 'normalizeparam'));
leadfieldopt = ft_setopt(leadfieldopt, 'weight',         ft_getopt(varargin, 'weight'));

% convert the yes/no arguments to the corresponding logical values
fixedori       = istrue(fixedori);
keepcsd        = istrue(keepcsd);  % see below
keepfilter     = istrue(keepfilter);
keepleadfield  = istrue(keepleadfield);
keepmom        = istrue(keepmom);
projectnoise   = istrue(projectnoise);
realfilter     = istrue(realfilter);

if ~isequal(weightnorm, 'no')
  ft_error('weight normalization is not supported');
end

if ~keepcsd
  ft_warning('setting keepcsd to ''yes'', overruling specified option');
  keepcsd = true; % this needs to be true by definition
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

% flags to avoid calling isfield repeatedly in the loop over grid positions (saves a lot of time)
hasmom        = isfield(dip, 'mom');
hasleadfield  = isfield(dip, 'leadfield');
hasfilter     = isfield(dip, 'filter');

% keep the original details on inside and outside positions
originside = dip.inside;
origpos    = dip.pos;

% select only the dipole positions inside the brain for scanning
dip.pos    = dip.pos(originside,:);
dip.inside = true(size(dip.pos,1),1);

if hasmom
  dip.mom = dip.mom(:, originside);
end
if hasleadfield
  ft_info('using precomputed leadfields\n');
  dip.leadfield = dip.leadfield(originside);
end
if hasfilter
  ft_info('using precomputed filters\n');
  dip.filter = dip.filter(originside);
end

if ~isempty(refdip)
  rf = ft_compute_leadfield(refdip, grad, headmodel, leadfieldopt{:});
else
  rf = [];
end

if ~isempty(supdip)
  sf = ft_compute_leadfield(supdip, grad, headmodel, leadfieldopt{:});
else
  sf = [];
end

% sanity check
if (~isempty(rf) || ~isempty(sf)) && hasfilter
  ft_error('precomputed filters cannot be used in combination with a refdip or supdip')
end

megchan  = setdiff(1:size(C,1), [refchan supchan]);
Nrefchan = length(refchan);
Nsupchan = length(supchan);
Nmegchan = length(megchan);
Cmeg     = C(megchan,megchan);   %  the filter uses the csd between all MEG channels

% it is difficult to give a quantitative estimate of lambda, therefore also
% support relative (percentage) measure that can be specified as string (e.g. '10%')
% the converted value needs to be passed on to ft_inv
lambda = ft_getopt(invopt, 'lambda');
if ~isempty(lambda) && ischar(lambda) && lambda(end)=='%'
  ratio  = sscanf(lambda, '%f%%');
  ratio  = ratio/100;
  lambda = ratio * trace(C)/size(C,1);
  invopt = ft_setopt(invopt, 'lambda', lambda);
end

if projectnoise
  % estimate the noise level in the covariance matrix by the smallest singular (non-zero) value
  noise = svd(Cmeg);
  noise = noise(rank(Cmeg));
  % estimated noise floor is equal to or higher than a numeric lambda
  noise = max(noise, lambda);
end

if realfilter
  % construct the filter only on the real part of the CSD matrix, i.e. filter is real
  invCmeg = ft_inv(real(Cmeg), invopt{:});
else
  % construct the filter on the complex CSD matrix, i.e. filter contains imaginary component as well
  % this results in a phase rotation of the channel data if the filter is applied to the data
  invCmeg = ft_inv(Cmeg, invopt{:});
end

% start the scanning with the proper metric
ft_progress('init', feedback, 'beaming sources');

for i=1:size(dip.pos,1)
  if hasfilter
    % precomputed filter is provided, the leadfield is not needed
  elseif hasleadfield && isfield(dip, 'mom') && size(dip.mom, 1)==size(dip.leadfield{i}, 2)
    % reuse the leadfield that was previously computed and project
    lf = dip.leadfield{i} * dip.mom(:,i);
  elseif  hasleadfield &&  isfield(dip, 'mom')
    % reuse the leadfield that was previously computed but don't project
    lf = dip.leadfield{i};
  elseif hasleadfield && ~isfield(dip, 'mom')
    % reuse the leadfield that was previously computed
    lf = dip.leadfield{i};
  elseif ~hasleadfield && isfield(dip, 'mom')
    % compute the leadfield for a fixed dipole orientation
    lf = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, leadfieldopt{:}) * dip.mom(:,i);
  else
    % compute the leadfield
    lf = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, leadfieldopt{:});
  end
  
  % determine the dimensionality of the source (usually 3, or 6 for a dipole pair)
  if hasfilter
    sourcedim = size(dip.filter{i},1);
  else
    % concatenate scandip, refdip and supdip
    lfa = [lf rf sf];
    sourcedim = size(lfa,2);
  end
  
  if fixedori
    if isempty(refdip) && isempty(supdip) && isempty(refchan) && isempty(supchan)
      % compute the leadfield for the optimal dipole orientation
      % subsequently the leadfield for only that dipole orientation will
      % be used for the final filter computation
      if hasfilter && size(dip.filter{i},1)==1
        % nothing to do
        ft_warning('Ignoring ''fixedori''. The fixedori option is supported only if there is ONE dipole for location.')
      else
        if hasfilter && size(dip.filter{i},1)~=1
          filt = dip.filter{i};
        else
          filt = pinv(lfa' * invCmeg * lfa) * lfa' * invCmeg;
        end
        [u, s, v] = svd(real(filt * Cmeg * ctranspose(filt)));
        maxpowori = u(:,1);
        if numel(s)>1
          eta = s(1,1)./s(2,2);
        else
          eta = nan;
        end
        lfa  = lfa * maxpowori;
        dipout.ori{i} = maxpowori;
        dipout.eta(i) = eta;
        % update the number of dipole components
        sourcedim = size(lfa,2);
      end
    else
      ft_warning('Ignoring ''fixedori''. The fixedori option is supported only if there is ONE dipole for location.')
    end
  end
  
  if hasfilter
    % use the provided filter
    filt = dip.filter{i};
  else
    % construct the spatial filter
    filt = pinv(lfa' * invCmeg * lfa) * lfa' * invCmeg;                % use PINV/SVD to cover rank deficient leadfield
  end
  
  % concatenate the source filters with the channel filters
  filtn = zeros(sourcedim+Nrefchan+Nsupchan, Nmegchan+Nrefchan+Nsupchan);
  % this part of the filter relates to the sources
  filtn(1:sourcedim,megchan) = filt;
  % this part of the filter relates to the channels
  filtn((sourcedim+1):end,setdiff(1:(Nmegchan+Nrefchan+Nsupchan), megchan)) = eye(Nrefchan+Nsupchan);
  filt = filtn;
  clear filtn
  
  if keepcsd
    dipout.csd{i,1} = filt * C * ctranspose(filt);
  end
  if projectnoise
    dipout.noisecsd{i,1} = noise * (filt * ctranspose(filt));
  end
  if keepmom && ~isempty(dat)
    dipout.mom{i,1} = filt * dat;
  end
  if keepfilter
    dipout.filter{i,1} = filt;
  end
  if keepleadfield && needleadfield
    dipout.leadfield{i,1} = lf;
  end
  
  ft_progress(i/size(dip.pos,1), 'beaming source %d from %d\n', i, size(dip.pos,1));
  
  % remember how all components in the output csd should be interpreted
  %scandiplabel = repmat({'scandip'}, 1, size(lf, 2));    % based on last leadfield
  scandiplabel = repmat({'scandip'}, 1, size(filt, 1)-size(rf, 2)-size(sf, 2)-Nrefchan-Nsupchan); % robust if lf does not exist
  refdiplabel  = repmat({'refdip'},  1, size(rf, 2));
  supdiplabel  = repmat({'supdip'},  1, size(sf, 2));
  refchanlabel = repmat({'refchan'}, 1, Nrefchan);
  supchanlabel = repmat({'supchan'}, 1, Nsupchan);
  % concatenate all the labels
  dipout.csdlabel{i,1} = [scandiplabel refdiplabel supdiplabel refchanlabel supchanlabel];
end % for all dipoles

ft_progress('close');

% wrap it all up, prepare the complete output
dipout.inside  = originside;
dipout.pos     = origpos;

fnames_cell   = {'leadfield' 'filter' 'ori' 'csd' 'noisecsd' 'subspace' 'mom' 'csdlabel'};
for k = 1:numel(fnames_cell)
  if isfield(dipout, fnames_cell{k})
    dipout.(fnames_cell{k})( originside) = dipout.(fnames_cell{k});
    dipout.(fnames_cell{k})(~originside) = {[]};
  end
end

fnames_scalar = {'pow' 'noise' 'eta' 'coh'};
for k = 1:numel(fnames_scalar)
  if isfield(dipout, fnames_scalar{k})
    dipout.(fnames_scalar{k})( originside) = dipout.(fnames_scalar{k});
    dipout.(fnames_scalar{k})(~originside) = nan;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to compute the pseudo inverse. This is the same as the
% standard MATLAB function, except that the default tolerance is twice as
% high.
%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision$  $Date: 2009/01/07 13:12:03 $
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
