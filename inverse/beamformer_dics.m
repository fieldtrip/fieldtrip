function [dipout] = beamformer_dics(dip, grad, headmodel, dat, Cf, varargin)

% BEAMFORMER_DICS scans on pre-defined dipole locations with a single dipole
% and returns the beamformer spatial filter output for a dipole on every
% location. Dipole locations that are outside the head will return a
% NaN value.
%
% Use as
%   [dipout] = beamformer_dics(dipin, grad, headmodel, dat, cov, varargin)
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
%  'Pr'               = power of the external reference channel
%  'Cr'               = cross spectral density between all data channels and the external reference channel
%  'refdip'           = location of dipole with which coherence is computed
%  'lambda'           = regularisation parameter
%  'powmethod'        = can be 'trace' or 'lambda1'
%  'feedback'         = give ft_progress indication, can be 'text', 'gui' or 'none'
%  'fixedori'         = use fixed or free orientation,                 can be 'yes' or 'no'
%  'projectnoise'     = project noise estimate through filter,         can be 'yes' or 'no'
%  'realfilter'       = construct a real-valued filter,                can be 'yes' or 'no'
%  'keepfilter'       = remember the beamformer filter,                can be 'yes' or 'no'
%  'keepleadfield'    = remember the forward computation,              can be 'yes' or 'no'
%  'keepcsd'          = remember the estimated cross-spectral density, can be 'yes' or 'no'
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
Pr             = ft_getopt(varargin, 'Pr');
Cr             = ft_getopt(varargin, 'Cr');
refdip         = ft_getopt(varargin, 'refdip');
powmethod      = ft_getopt(varargin, 'powmethod');  % the default for this is set below
realfilter     = ft_getopt(varargin, 'realfilter'); % the default for this is set below
subspace       = ft_getopt(varargin, 'subspace');
% these settings pertain to the forward model, the defaults are set in compute_leadfield
reducerank     = ft_getopt(varargin, 'reducerank');
normalize      = ft_getopt(varargin, 'normalize');
normalizeparam = ft_getopt(varargin, 'normalizeparam');
% these optional settings have defaults
feedback       = ft_getopt(varargin, 'feedback', 'text');
keepcsd        = ft_getopt(varargin, 'keepcsd', 'no');
keepfilter     = ft_getopt(varargin, 'keepfilter', 'no');
keepleadfield  = ft_getopt(varargin, 'keepleadfield', 'no');
lambda         = ft_getopt(varargin, 'lambda', 0);
projectnoise   = ft_getopt(varargin, 'projectnoise', 'yes');
fixedori       = ft_getopt(varargin, 'fixedori', 'no');

% convert the yes/no arguments to the corresponding logical values
keepcsd        = strcmp(keepcsd,       'yes');
keepfilter     = strcmp(keepfilter,    'yes');
keepleadfield  = strcmp(keepleadfield, 'yes');
projectnoise   = strcmp(projectnoise,  'yes');
fixedori       = strcmp(fixedori,      'yes');
dofeedback     = ~strcmp(feedback,     'none');
% FIXME besides regular/complex lambda1, also implement a real version

% default is to use the largest singular value of the csd matrix, see Gross 2001
if isempty(powmethod)
  powmethod = 'lambda1';
end

% default is to be consistent with the original description of DICS in Gross 2001
if isempty(realfilter)
  realfilter = 'no';
end

% use these two logical flags instead of doing the string comparisons each time again
powtrace   = strcmp(powmethod, 'trace');
powlambda1 = strcmp(powmethod, 'lambda1');

if ~isempty(Cr)
  % ensure that the cross-spectral density with the reference signal is a column matrix
  Cr = Cr(:);
end

if isfield(dip, 'mom') && fixedori
  ft_error('you cannot specify a dipole orientation and fixedmom simultaneously');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the dipole positions that are inside/outside the brain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(dip, 'inside')
  dip.inside = ft_inside_vol(dip.pos, headmodel);
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
hasmom        = false;
hasleadfield  = false;
hasfilter     = false;
hassubspace   = false;

% select only the dipole positions inside the brain for scanning
dip.pos    = dip.pos(originside,:);
dip.inside = true(size(dip.pos,1),1);
if isfield(dip, 'mom')
  hasmom = 1;
  dip.mom = dip.mom(:,originside);
end
if isfield(dip, 'leadfield')
  hasleadfield = 1;
  if dofeedback
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

% dics has the following sub-methods, which depend on the function input arguments
% power only, cortico-muscular coherence and cortico-cortical coherence
if ~isempty(Cr) && ~isempty(Pr) && isempty(refdip)
  % compute cortico-muscular coherence, using reference cross spectral density
  submethod = 'dics_refchan';
elseif isempty(Cr) && isempty(Pr) && ~isempty(refdip)
  % compute cortico-cortical coherence with a dipole at the reference position
  submethod = 'dics_refdip';
elseif isempty(Cr) && isempty(Pr) && isempty(refdip)
  % only compute power of a dipole at the grid positions
  submethod = 'dics_power';
else
  ft_error('invalid combination of input arguments for dics');
end

isrankdeficient = (rank(Cf)<size(Cf,1));
rankCf = rank(Cf);

% it is difficult to give a quantitative estimate of lambda, therefore also
% support relative (percentage) measure that can be specified as string (e.g. '10%')
if ~isempty(lambda) && ischar(lambda) && lambda(end)=='%'
  ratio = sscanf(lambda, '%f%%');
  ratio = ratio/100;
  if ~isempty(subspace) && numel(subspace)>1,
    lambda = ratio * trace(subspace*Cf*subspace')./size(subspace,1);
  else
    lambda = ratio * trace(Cf)/size(Cf,1);
  end
end

if projectnoise
  % estimate the noise level in the covariance matrix by the smallest (non-zero) singular value
  noise = svd(Cf);
  noise = noise(rankCf);
  % estimated noise floor is equal to or higher than lambda
  noise = max(noise, lambda);
end

% the inverse only has to be computed once for all dipoles
if strcmp(realfilter, 'yes')
  % the filter is computed using only the leadfield and the inverse covariance or CSD matrix
  % therefore using the real-valued part of the CSD matrix here ensures a real-valued filter
  invCf = pinv(real(Cf) + lambda * eye(size(Cf)));
else
  invCf = pinv(Cf + lambda * eye(size(Cf)));
end

if hassubspace
  if dofeedback
    fprintf('using source-specific subspace projection\n');
  end
  % remember the original data prior to the voxel dependent subspace projection
  dat_pre_subspace = dat;
  Cf_pre_subspace = Cf;
  if strcmp(submethod, 'dics_refchan')
    Cr_pre_subspace = Cr;
    Pr_pre_subspace = Pr;
  end
elseif ~isempty(subspace)
  if dofeedback
    fprintf('using data-specific subspace projection\n');
  end
  % TODO implement an "eigenspace beamformer" as described in Sekihara et al. 2002 in HBM
  if numel(subspace)==1,
    % interpret this as a truncation of the eigenvalue-spectrum
    % if <1 it is a fraction of the largest eigenvalue
    % if >=1 it is the number of largest eigenvalues
    dat_pre_subspace = dat;
    Cf_pre_subspace  = Cf;
    [u, s, v] = svd(real(Cf));
    if subspace<1,
      sel      = find(diag(s)./s(1,1) > subspace);
      subspace = max(sel);
    end
    
    Cf       = s(1:subspace,1:subspace);
    % this is equivalent to subspace*Cf*subspace' but behaves well numerically
    % by construction.
    invCf    = diag(1./diag(Cf + lambda * eye(size(Cf))));
    subspace = u(:,1:subspace)';
    if ~isempty(dat), dat = subspace*dat; end
    
    if strcmp(submethod, 'dics_refchan')
      Cr = subspace*Cr;
    end
    
  else
    Cf_pre_subspace  = Cf;
    Cf    = subspace*Cf*subspace'; % here the subspace can be different from
    % the singular vectors of Cy, so we have to do the sandwiching as opposed
    % to line 216
    if strcmp(realfilter, 'yes')
      invCf = pinv(real(Cf) + lambda * eye(size(Cf)));
    else
      invCf = pinv(Cf + lambda * eye(size(Cf)));
    end
    
    if strcmp(submethod, 'dics_refchan')
      Cr = subspace*Cr;
    end
  end
end

% start the scanning with the proper metric
ft_progress('init', feedback, 'scanning grid');
switch submethod
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % dics_power
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'dics_power'
    % only compute power of a dipole at the grid positions
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
      elseif ~hasleadfield && hasmom
        % compute the leadfield for a fixed dipole orientation
        lf = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam) * dip.mom(:,i);
      else
        % compute the leadfield
        lf = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam);
      end
      if hassubspace
        % do subspace projection of the forward model
        lf = dip.subspace{i} * lf;
        % the cross-spectral density becomes voxel dependent due to the projection
        Cf    = dip.subspace{i} * Cf_pre_subspace * dip.subspace{i}';
        if strcmp(realfilter, 'yes')
          invCf = pinv(dip.subspace{i} * (real(Cf_pre_subspace) + lambda * eye(size(Cf_pre_subspace))) * dip.subspace{i}');
        else
          invCf = pinv(dip.subspace{i} * (Cf_pre_subspace + lambda * eye(size(Cf_pre_subspace))) * dip.subspace{i}');
        end
      elseif ~isempty(subspace)
        % do subspace projection of the forward model only
        lforig = lf;
        lf     = subspace * lf;
        
        % according to Kensuke's paper, the eigenspace bf boils down to projecting
        % the 'traditional' filter onto the subspace
        % spanned by the first k eigenvectors [u,s,v] = svd(Cy); filt = ESES*filt;
        % ESES = u(:,1:k)*u(:,1:k)';
        % however, even though it seems that the shape of the filter is identical to
        % the shape it is obtained with the following code, the w*lf=I does not
        % hold.
      end
      
      if hasfilter
        % use precomputed filter
        filt = dip.filter{i};
      else
        % compute filter
        filt = pinv(lf' * invCf * lf) * lf' * invCf;              % Gross eqn. 3, use PINV/SVD to cover rank deficient leadfield
      end
      
      if fixedori
        % use single dipole orientation
        if hasfilter && size(filt,1) == 1
          % provided precomputed filter already projects to one
          % orientation, nothing to be done here
        else
          % find out the optimal dipole orientation
          [u, s, v] = svd(real(filt * Cf * ctranspose(filt)));
          maxpowori = u(:,1);
          eta = s(1,1)./s(2,2);
          
          % and compute the leadfield for that orientation
          lf  = lf * maxpowori;
          dipout.ori{i} = maxpowori;
          dipout.eta(i) = eta;
          if ~isempty(subspace), lforig = lforig * maxpowori; end
          
          % recompute the filter to only use that orientation
          filt = pinv(lf' * invCf * lf) * lf' * invCf;
        end
      elseif hasfilter && size(filt,1) == 1
        ft_error('the precomputed filter you provided projects to a single dipole orientation, but you request fixedori=''no''; this is invalid. Either provide a filter with the three orientations retained, or specify fixedori=''yes''.');
      end
      
      csd = filt * Cf * ctranspose(filt);                         % Gross eqn. 4 and 5
      if powlambda1
        if size(csd,1) == 1
          % only 1 orientation, no need to do svd
          dipout.pow(i,1) = real(csd);
        else
          dipout.pow(i,1) = lambda1(csd);                           % compute the power at the dipole location, Gross eqn. 8
        end
      elseif powtrace
        dipout.pow(i,1) = real(trace(csd));                         % compute the power at the dipole location
      end
      if keepcsd
        dipout.csd{i,1} = csd;
      end
      if projectnoise
        if powlambda1
          dipout.noise(i,1) = noise * lambda1(filt * ctranspose(filt));
        elseif powtrace
          dipout.noise(i,1) = noise * real(trace(filt * ctranspose(filt)));
        end
        if keepcsd
          dipout.noisecsd{i,1} = noise * filt * ctranspose(filt);
        end
      end
      if keepfilter
        if ~isempty(subspace)
          dipout.filter{i,1} = filt*subspace; %FIXME should this be subspace, or pinv(subspace)?
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dics_refchan
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'dics_refchan'
    % compute cortico-muscular coherence, using reference cross spectral density
    for i=1:size(dip.pos,1)
      if hasleadfield
        % reuse the leadfield that was previously computed
        lf = dip.leadfield{i};
      elseif hasmom
        % compute the leadfield for a fixed dipole orientation
        lf = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, 'reducerank', reducerank, 'normalize', normalize) .* dip.mom(i,:)';
      else
        % compute the leadfield
        lf = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, 'reducerank', reducerank, 'normalize', normalize);
      end
      if hassubspace
        % do subspace projection of the forward model
        lforig = lf;
        lf = dip.subspace{i} * lf;
        % the cross-spectral density becomes voxel dependent due to the projection
        Cf    = dip.subspace{i} * Cf_pre_subspace * dip.subspace{i}';
        invCf = pinv(dip.subspace{i} * (Cf_pre_subspace + lambda * eye(size(Cf))) * dip.subspace{i}');
      elseif ~isempty(subspace)
        % do subspace projection of the forward model only
        lforig = lf;
        lf     = subspace * lf;
        
        % according to Kensuke's paper, the eigenspace bf boils down to projecting
        % the 'traditional' filter onto the subspace
        % spanned by the first k eigenvectors [u,s,v] = svd(Cy); filt = ESES*filt;
        % ESES = u(:,1:k)*u(:,1:k)';
        % however, even though it seems that the shape of the filter is identical to
        % the shape it is obtained with the following code, the w*lf=I does not
        % hold.
      end
      
      if hasfilter
        % use precomputed filter
        filt = dip.filter{i};
      else
        % compute filter
        filt = pinv(lf' * invCf * lf) * lf' * invCf;              % Gross eqn. 3, use PINV/SVD to cover rank deficient leadfield
      end
      
      if fixedori
        % use single dipole orientation
        if hasfilter && size(filt,1) == 1
          % provided precomputed filter already projects to one
          % orientation, nothing to be done here
        else
          % find out the optimal dipole orientation
          [u, s, v] = svd(real(filt * Cf * ctranspose(filt)));
          maxpowori = u(:,1);
          
          % compute the leadfield for that orientation
          lf  = lf * maxpowori;
          dipout.ori{i,1} = maxpowori;
          
          % recompute the filter to only use that orientation
          filt = pinv(lf' * invCf * lf) * lf' * invCf;
        end
      elseif hasfilter && size(filt,1) == 1
        ft_error('the precomputed filter you provided projects to a single dipole orientation, but you request fixedori=''no''; this is invalid. Either provide a filter with the three orientations retained, or specify fixedori=''yes''.');
      end
      
      if powlambda1
        [pow, ori] = lambda1(filt * Cf * ctranspose(filt));            % compute the power and orientation at the dipole location, Gross eqn. 4, 5 and 8
      elseif powtrace
        pow = real(trace(filt * Cf * ctranspose(filt)));               % compute the power at the dipole location
      end
      csd = filt*Cr;                                                   % Gross eqn. 6
      if powlambda1
        % FIXME this should use the dipole orientation with maximum power
        coh = lambda1(csd)^2 / (pow * Pr);                             % Gross eqn. 9
      elseif powtrace
        coh = norm(csd)^2 / (pow * Pr);
      end
      dipout.pow(i,1) = pow;
      dipout.coh(i,1) = coh;
      if keepcsd
        dipout.csd{i,1} = csd;
      end
      if projectnoise
        if powlambda1
          dipout.noise(i,1) = noise * lambda1(filt * ctranspose(filt));
        elseif powtrace
          dipout.noise(i,1) = noise * real(trace(filt * ctranspose(filt)));
        end
        if keepcsd
          dipout.noisecsd{i,1} = noise * filt * ctranspose(filt);
        end
      end
      if keepfilter
        dipout.filter{i,1} = filt;
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dics_refdip
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'dics_refdip'
    if hassubspace || ~isempty(subspace)
      ft_error('subspace projections are not supported for beaming cortico-cortical coherence');
    end
    if fixedori
      ft_error('fixed orientations are not supported for beaming cortico-cortical coherence');
    end
    if isstruct(refdip) && isfield(refdip, 'filter') % check if precomputed filter is present
      assert(iscell(refdip.filter) && numel(refdip.filter)==1);
      filt1 = refdip.filter{1};
    elseif isstruct(refdip) && isfield(refdip, 'leadfield') % check if precomputed leadfield is present
      assert(iscell(refdip.leadfield) && numel(refdip.leadfield)==1);
      lf1 = refdip.leadfield{1};
      filt1 = pinv(lf1' * invCf * lf1) * lf1' * invCf;       % use PINV/SVD to cover rank deficient leadfield
    elseif isstruct(refdip) && isfield(refdip, 'pos') % check if only position of refdip is present
      assert(isnumeric(refdip.pos) && numel(refdip.pos)==3);
      lf1 = ft_compute_leadfield(refdip.pos, grad, headmodel, 'reducerank', reducerank, 'normalize', normalize);
      if isfield(refdip,'mom'); % check for fixed orientation
        lf1 = lf1.*refdip.mom(:); 
      end 
      filt1 = pinv(lf1' * invCf * lf1) * lf1' * invCf;       % use PINV/SVD to cover rank deficient leadfield
    else % backwards compatible with previous implementation - only position of refdip is present
      % compute cortio-cortical coherence with a dipole at the reference position
      lf1 = ft_compute_leadfield(refdip, grad, headmodel, 'reducerank', reducerank, 'normalize', normalize);
      % construct the spatial filter for the first (reference) dipole location
      filt1 = pinv(lf1' * invCf * lf1) * lf1' * invCf;       % use PINV/SVD to cover rank deficient leadfield
    end
    if powlambda1
      Pref = lambda1(filt1 * Cf * ctranspose(filt1));      % compute the power at the first dipole location, Gross eqn. 8
    elseif powtrace
      Pref = real(trace(filt1 * Cf * ctranspose(filt1)));  % compute the power at the first dipole location
    end
    for i=1:size(dip.pos,1)
      if hasleadfield
        % reuse the leadfield that was previously computed
        lf2 = dip.leadfield{i};
      elseif hasmom
        % compute the leadfield for a fixed dipole orientation
        lf2 = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, 'reducerank', reducerank, 'normalize', normalize) .* dip.mom(i,:)';
      else
        % compute the leadfield
        lf2 = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, 'reducerank', reducerank, 'normalize', normalize);
      end
      if hasfilter
        % use the provided filter
        filt2 = dip.filter{i};
      else
        % construct the spatial filter for the second dipole location
        filt2 = pinv(lf2' * invCf * lf2) * lf2' * invCf;   %  use PINV/SVD to cover rank deficient leadfield
      end
      csd = filt1 * Cf * ctranspose(filt2);                % compute the cross spectral density between the two dipoles, Gross eqn. 4
      if powlambda1
        pow = lambda1(filt2 * Cf * ctranspose(filt2));     % compute the power at the second dipole location, Gross eqn. 8
      elseif powtrace
        pow = real(trace(filt2 * Cf * ctranspose(filt2))); % compute the power at the second dipole location
      end
      if powlambda1
        coh = lambda1(csd)^2 / (pow * Pref);               % compute the coherence between the first and second dipole
      elseif powtrace
        coh = real(trace((csd)))^2 / (pow * Pref);         % compute the coherence between the first and second dipole
      end
      dipout.pow(i,1) = pow;
      dipout.coh(i,1) = coh;
      if keepcsd
        dipout.csd{i,1} = csd;
      end
      if projectnoise
        if powlambda1
          dipout.noise(i,1) = noise * lambda1(filt2 * ctranspose(filt2));
        elseif powtrace
          dipout.noise(i,1) = noise * real(trace(filt2 * ctranspose(filt2)));
        end
        if keepcsd
          dipout.noisecsd{i,1} = noise * filt2 * ctranspose(filt2);
        end
      end
      if keepleadfield
        dipout.leadfield{i,1} = lf2;
      end
      ft_progress(i/size(dip.pos,1), 'scanning grid %d/%d\n', i, size(dip.pos,1));
    end
    
end % switch submethod

ft_progress('close');

% wrap it all up, prepare the complete output
dipout.inside  = originside;
dipout.pos     = origpos;

% reassign the scan values over the inside and outside grid positions
if isfield(dipout, 'leadfield')
  dipout.leadfield( originside) = dipout.leadfield;
  dipout.leadfield(~originside) = {[]};
end
if isfield(dipout, 'filter')
  dipout.filter( originside) = dipout.filter;
  dipout.filter(~originside) = {[]};
end
if isfield(dipout, 'ori')
  dipout.ori( originside) = dipout.ori;
  dipout.ori(~originside) = {[]};
end
if isfield(dipout, 'eta')
  dipout.eta( originside) = dipout.eta;
  dipout.eta(~originside) = nan;
end
if isfield(dipout, 'pow')
  dipout.pow( originside) = dipout.pow;
  dipout.pow(~originside) = nan;
end
if isfield(dipout, 'noise')
  dipout.noise( originside) = dipout.noise;
  dipout.noise(~originside) = nan;
end
if isfield(dipout, 'coh')
  dipout.coh( originside) = dipout.coh;
  dipout.coh(~originside) = nan;
end
if isfield(dipout, 'csd')
  dipout.csd( originside) = dipout.csd;
  dipout.csd(~originside) = {[]};
end
if isfield(dipout, 'noisecsd')
  dipout.noisecsd( originside) = dipout.noisecsd;
  dipout.noisecsd(~originside) = {[]};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to obtain the largest singular value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s, ori] = lambda1(x)
% determine the largest singular value, which corresponds to the power along the dominant direction
[u, s, v] = svd(x);
s   = s(1);
ori = u(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to compute the pseudo inverse. This is the same as the
% standard MATLAB function, except that the default tolerance is twice as
% high.
%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision$  $Date: 2009/06/17 13:40:37 $
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

