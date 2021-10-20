function [estimate] = ft_inverse_dics(sourcemodel, sens, headmodel, dat, C, varargin)

% FT_INVERSE_DICS scans on pre-defined dipole locations with a single dipole
% and returns the beamformer spatial filter output for a dipole on every
% location.
%
% Use as
%   [estimate] = ft_inverse_dics(sourcemodel, sens, headmodel, dat, cov, ...)
% where
%   sourcemodel is the input source model, see FT_PREPARE_SOURCEMODEL
%   sens        is the gradiometer or electrode definition, see FT_DATATYPE_SENS
%   headmodel   is the volume conductor definition, see FT_PREPARE_HEADMODEL
%   dat         is the data matrix with the ERP or ERF
%   cov         is the data covariance or cross-spectral density matrix
% and
%   estimate    contains the estimated source parameters
%
% Additional input arguments should be specified as key-value pairs and can include
%   'Pr'               = power of the external reference channel
%   'Cr'               = cross spectral density between all data channels and the external reference channel
%   'refdip'           = location of dipole with which coherence is computed
%   'powmethod'        = can be 'trace' or 'lambda1'
%   'feedback'         = can be 'none', 'gui', 'dial', 'textbar', 'text', 'textcr', 'textnl' (default = 'text')
%   'fixedori'         = use fixed or free orientation,                 can be 'yes' or 'no'
%   'projectnoise'     = project noise estimate through filter,         can be 'yes' or 'no'
%   'realfilter'       = construct a real-valued filter,                can be 'yes' or 'no'
%   'keepfilter'       = remember the beamformer filter,                can be 'yes' or 'no'
%   'keepleadfield'    = remember the forward computation,              can be 'yes' or 'no'
%   'keepcsd'          = remember the estimated cross-spectral density, can be 'yes' or 'no'
%   'weightnorm'       = normalize the beamformer weights,              can be 'no', 'unitnoisegain' or 'nai'
%
% These options influence the forward computation of the leadfield
%   'reducerank'      = 'no' or number (default = 3 for EEG, 2 for MEG)
%   'backproject'     = 'yes' or 'no', in the case of a rank reduction this parameter determines whether the result will be backprojected onto the original subspace (default = 'yes')
%   'normalize'       = 'no', 'yes' or 'column' (default = 'no')
%   'normalizeparam'  = parameter for depth normalization (default = 0.5)
%   'weight'          = number or Nx1 vector, weight for each dipole position to compensate for the size of the corresponding patch (default = 1)
%
% These options influence the mathematical inversion of the cross-spectral density matrix
%  'lambda'           = regularisation parameter
%  'kappa'            = parameter for covariance matrix inversion
%  'tol'              = parameter for covariance matrix inversion
%
% If the dipole definition only specifies the dipole location, a rotating
% dipole (regional source) is assumed on each location. If a dipole moment
% is specified, its orientation will be used and only the strength will
% be fitted to the data.
%
% See also FT_SOURCEANALYSIS, FT_PREPARE_HEADMODEL, FT_PREPARE_SOURCEMODEL

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

% get the optional input arguments, or use defaults
Pr             = ft_getopt(varargin, 'Pr');
Cr             = ft_getopt(varargin, 'Cr');
refdip         = ft_getopt(varargin, 'refdip');
powmethod      = ft_getopt(varargin, 'powmethod');  % the default for this is set below
realfilter     = ft_getopt(varargin, 'realfilter'); % the default for this is set below
subspace       = ft_getopt(varargin, 'subspace');
feedback       = ft_getopt(varargin, 'feedback',      'text');
keepcsd        = ft_getopt(varargin, 'keepcsd',       'no');
keepfilter     = ft_getopt(varargin, 'keepfilter',    'no');
keepleadfield  = ft_getopt(varargin, 'keepleadfield', 'no');
projectnoise   = ft_getopt(varargin, 'projectnoise',  'yes');
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
keepfilter     = istrue(keepfilter);
keepleadfield  = istrue(keepleadfield);
keepcsd        = istrue(keepcsd);
projectnoise   = istrue(projectnoise);
fixedori       = istrue(fixedori);

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

if isfield(sourcemodel, 'mom') && fixedori
  ft_error('you cannot specify a dipole orientation and fixedori simultaneously');
end

% flags to avoid calling isfield repeatedly in the loop over grid positions (saves a lot of time)
hasmom        = isfield(sourcemodel, 'mom');
hasleadfield  = isfield(sourcemodel, 'leadfield');
hasfilter     = isfield(sourcemodel, 'filter');
hassubspace   = isfield(sourcemodel, 'subspace');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the dipole positions that are inside/outside the brain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(sourcemodel, 'inside')
  if hasfilter
    sourcemodel.inside = ~cellfun(@isempty, sourcemodel.filter);
  elseif hasleadfield
    sourcemodel.inside = ~cellfun(@isempty, sourcemodel.leadfield);
  else
    sourcemodel.inside = ft_inside_headmodel(sourcemodel.pos, headmodel);
  end
end

% convert to logical representation
sourcemodel = fixinside(sourcemodel);

if hasfilter && (fixedori || ~isequal(weightnorm, 'no'))
  ft_warning('with precomputed spatial filters a fixed orientation constraint or weight normalisation options are not applied');
end

% keep the original details on inside and outside positions
originside = sourcemodel.inside;
origpos    = sourcemodel.pos;

% select only the dipole positions inside the brain for scanning
sourcemodel.pos    = sourcemodel.pos(originside,:);
sourcemodel.inside = true(size(sourcemodel.pos,1),1);

if hasmom
  sourcemodel.mom = sourcemodel.mom(:,originside);
end

if hasfilter
  ft_info('using precomputed filters\n');
  sourcemodel.filter = sourcemodel.filter(originside);
elseif hasleadfield
  ft_info('using precomputed leadfields\n');
  sourcemodel.leadfield = sourcemodel.leadfield(originside);
else
  ft_info('computing forward model on the fly\n');
end

if hassubspace
  ft_info('using subspace projection\n');
  sourcemodel.subspace = sourcemodel.subspace(originside);
end

% dics has the following sub-methods, which depend on the function input arguments
% power only, cortico-muscular coherence and cortico-cortical coherence
if ~isempty(Cr) && ~isempty(Pr) && isempty(refdip)
  % compute cortico-muscular coherence, using reference cross spectral density
  submethod = 'dics_refchan';
elseif isempty(Cr) && isempty(Pr) && ~isempty(refdip)
  % compute cortico-cortical coherence with a dipole at the reference position
  submethod = 'dics_refdip';
  % check for forbidden options in combination with this submethod
  if hassubspace || ~isempty(subspace)
    ft_error('subspace projections are not supported for beaming cortico-cortical coherence');
  end
  if fixedori
    ft_error('fixed orientations are not supported for beaming cortico-cortical coherence');
  end
elseif isempty(Cr) && isempty(Pr) && isempty(refdip)
  % only compute power of a dipole at the grid positions
  submethod = 'dics_power';
else
  ft_error('invalid combination of input arguments for dics');
end

% it is difficult to give a quantitative estimate of lambda, therefore also
% support relative (percentage) measure that can be specified as string (e.g. '10%')
% the converted value needs to be passed on to ft_inv
lambda = ft_getopt(invopt, 'lambda');
if ~isempty(lambda) && ischar(lambda) && lambda(end)=='%'
  ratio  = sscanf(lambda, '%f%%');
  ratio  = ratio/100;
  lambda = ratio * trace(C)/size(C,1);
  if ~hassubspace
    invopt = ft_setopt(invopt, 'lambda', lambda);
  end
end

if projectnoise || strcmp(weightnorm, 'nai')
  % estimate the noise level in the covariance matrix by the smallest (non-zero) singular value
  % always needed for the NAI weight normalization case
  noise = svd(C);
  noise = noise(rank(C));
  % estimated noise floor is equal to or higher than lambda
  noise = max(noise, lambda);
end

% the inverse of the cross-spectral density matrix only has to be computed once for all dipoles
if hassubspace
  ft_notice('using source-specific subspace projection\n');
  % remember the original data prior to the voxel dependent subspace projection
  dat_pre_subspace = dat;
  C_pre_subspace = C;
  if strcmp(submethod, 'dics_refchan')
    Cr_pre_subspace = Cr;
    Pr_pre_subspace = Pr;
  end
elseif ~isempty(subspace)
  ft_notice('using data-specific subspace projection\n');
  % TODO implement an "eigenspace beamformer" as described in Sekihara et al. 2002 in HBM
  if numel(subspace)==1
    % interpret this as a truncation of the eigenvalue-spectrum
    % if <1 it is a fraction of the largest eigenvalue
    % if >=1 it is the number of largest eigenvalues
    dat_pre_subspace = dat;
    C_pre_subspace  = C;
    [u, s, v] = svd(real(C));
    if subspace<1
      subspace  = find(diag(s)./s(1,1) > subspace, 1, 'last');
    end
    
    C       = s(1:subspace,1:subspace);
    % this is equivalent to subspace*C*subspace' but behaves well numerically
    % by construction.
    invC     = diag(1./diag(C + lambda * eye(size(C))));
    subspace = u(:,1:subspace)';
    if ~isempty(dat), dat = subspace*dat; end
    
    if strcmp(submethod, 'dics_refchan')
      Cr = subspace*Cr;
    end
    
  else
    C_pre_subspace  = C;
    C               = subspace*C*subspace';
    % here the subspace can be different from the singular vectors of C, so we
    % have to do the sandwiching as opposed to line 256
    if strcmp(realfilter, 'yes')
      invC = ft_inv(real(C), invopt{:});
    else
      invC = ft_inv(C, invopt{:});
    end
    
    if strcmp(submethod, 'dics_refchan')
      Cr = subspace*Cr;
    end
  end
else
  if strcmp(realfilter, 'yes')
    % the filter is computed using only the leadfield and the inverse covariance or CSD matrix
    % therefore using the real-valued part of the CSD matrix here ensures a real-valued filter
    invC = ft_inv(real(C), invopt{:});
  else
    invC = ft_inv(C, invopt{:});
  end
end

if ~hassubspace
  % compute the square of invC, which might be needed for unitnoisegain or NAI constraint
  invC_squared = invC^2;
end

if strcmp(submethod, 'dics_refdip')
  % handle the reference dipole, this needs to be done only once, since in this implementation
  % no fancy pairwise dipole source model is assumed
  if isstruct(refdip) && isfield(refdip, 'filter') % check if precomputed filter is present
    assert(iscell(refdip.filter) && numel(refdip.filter)==1);
    filt1 = refdip.filter{1};
  else
    if ~isstruct(refdip)
      % this is very old behavior: refdip is assumed to be a 1x3 position vector
      refpos     = refdip; clear refdip;
      refdip.pos = refpos;
    end
    % now refdip is always a struct
    if isfield(refdip, 'leadfield') % check if precomputed leadfield is present
      assert(iscell(refdip.leadfield) && numel(refdip.leadfield)==1);
      lf1 = refdip.leadfield{1};
    elseif isfield(refdip, 'pos')
      assert(isnumeric(refdip.pos) && numel(refdip.pos)==3);
      lf1 = ft_compute_leadfield(refdip.pos, sens, headmodel, leadfieldopt{:});
    end
    if isfield(refdip, 'mom') % check for fixed orientation
      lf1 = lf1 * refdip.mom(:);
    end
    
    % construct the spatial filter
    switch weightnorm
      case 'nai'
        % Van Veen's Neural Activity Index
        ft_error('vector version of nai weight normalization is not implemented');
      case 'unitnoisegain'
        % filt*filt' = I
        % Unit-noise gain minimum variance (aka Borgiotti-Kaplan) beamformer
        denom = pinv(lf1' * invC * lf1);
        gamma = denom * (lf1' * invC_squared * lf1) * denom;
        
        % compute the spatial filter, as per eqn. 4.85
        filt1 = diag(1./sqrt(diag(gamma))) * denom * lf1' * invC;
      case 'arraygain'
        % filt*lf = ||lf||, applies to scalar leadfield, and to one of the possibilities of the vector version, eqn. 4.75
        lfn   = lf1./norm(lf1);
        filt1 = pinv(lfn' * invC * lfn) * lfn' * invC; % S&N eqn. 4.09 (scalar version), and eqn. 4.75 (vector version)
      case {'unitgain' 'no'}
        % this is the 'standard' unit gain constraint spatial filter: filt*lf=I, applies both to vector and scalar leadfields
        filt1 = pinv(lf1' * invC * lf1) * lf1' * invC; % Gross eqn. 3 & van Veen eqn. 23, use PINV/SVD to cover rank deficient leadfield
      otherwise
        ft_error('unsupported option for weightnorm');
    end
  end
  
  % compute the power of the reference dipole location
  if powlambda1
    Pr = lambda1(filt1 * C * ctranspose(filt1));      % compute the power at the first dipole location, Gross eqn. 8
  elseif powtrace
    Pr = real(trace(filt1 * C * ctranspose(filt1)));  % compute the power at the first dipole location
  end
end % if dics_refdip

% start the scanning with the proper metric
ft_progress('init', feedback, 'scanning grid');
for i=1:size(sourcemodel.pos,1)
  ft_progress(i/size(sourcemodel.pos,1), 'scanning grid %d/%d\n', i, size(sourcemodel.pos,1));

  if hasfilter
    % precomputed filter is provided, the leadfield is not needed
    filt = sourcemodel.filter{i};

  else
    if hasleadfield && hasmom && size(sourcemodel.mom, 1)==size(sourcemodel.leadfield{i}, 2)
      % reuse the leadfield that was previously computed and project
      lf = sourcemodel.leadfield{i} * sourcemodel.mom(:,i);
    elseif  hasleadfield &&  hasmom
      % reuse the leadfield that was previously computed but don't project
      lf = sourcemodel.leadfield{i};
    elseif  hasleadfield && ~hasmom
      % reuse the leadfield that was previously computed
      lf = sourcemodel.leadfield{i};
    elseif ~hasleadfield &&  hasmom
      % compute the leadfield for a fixed dipole orientation
      lf = ft_compute_leadfield(sourcemodel.pos(i,:), sens, headmodel, leadfieldopt{:}) * sourcemodel.mom(:,i);
    else
      % compute the leadfield
      lf = ft_compute_leadfield(sourcemodel.pos(i,:), sens, headmodel, leadfieldopt{:});
    end
    
    if hassubspace
      % do subspace projection of the forward model
      lforig = lf;
      lf = sourcemodel.subspace{i} * lf;
      % the cross-spectral density becomes voxel dependent due to the projection
      C    = sourcemodel.subspace{i} * C_pre_subspace * sourcemodel.subspace{i}';
      if strcmp(realfilter, 'yes')
        invC = ft_inv(sourcemodel.subspace{i} * real(C_pre_subspace) * sourcemodel.subspace{i}', invopt{:});
      else
        invC = ft_inv(sourcemodel.subspace{i} * C_pre_subspace * sourcemodel.subspace{i}', invopt{:});
      end
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
          % (also applies nai weightnorm constraint), based on equation 4.47 from Sekihara & Nagarajan (2008)
          % the following is a reformulation of the generalized eigenproblem
          [v, d]        = eig(pinv(lf' * invC_squared *lf)*(lf' * invC *lf));
          [d, iv]       = sort(diag(d), 'descend');
          unitnoiseori  = v(:,iv(1));
          lf            = lf * unitnoiseori;
          estimate.ori{i} = unitnoiseori;
          estimate.eta(i) = d(1)./d(2); % ratio between largest and second largest eigenvalues
          if hassubspace, lforig = lforig * unitnoiseori; end
          
        case 'arraygain'
          % optimal orientation calculation for array-gain beamformer, Sekihara & Nagarajan eqn. 4.44
          [v, d]        = eig(pinv(lf' * invC *lf)*(lf' * lf));
          [d, iv]       = sort(diag(d), 'descend');
          arraygainori  = v(:,iv(1));
          lf            = lf * arraygainori;
          estimate.ori{i} = arraygainori;
          estimate.eta(i) = d(1)./d(2); % ratio between largest and second largest eigenvalues
          if hassubspace, lforig = lforig * arraygainori; end
          
        otherwise
          % compute the leadfield for the optimal dipole orientation that maximizes spatial filter output
          % subsequently the leadfield for only that dipole orientation will be used for the final filter computation
          % filt = pinv(lf' * invC * lf) * lf' * invC;
          % [u, s, v] = svd(real(filt * C * ctranspose(filt)));
          % in this step the filter computation is not necessary, use the quick way to compute the voxel level covariance (cf. van Veen 1997)
          [u, s, v]     = svd(real(pinv(lf' * invC *lf)));
          maxpowori     = u(:,1);
          lf            = lf * maxpowori;
          estimate.ori{i} = maxpowori;
          estimate.eta(i) = s(1,1)./s(2,2); % ratio between the first and second singular values
          if hassubspace, lforig = lforig * maxpowori; end
      end
    end
    
    % construct the spatial filter
    switch weightnorm
      case 'nai'
        % Van Veen's Neural Activity Index
        % below equation is equivalent to following:
        % filt = pinv(lf' * invC * lf) * lf' * invC;
        % filt = filt/sqrt(noise*filt*filt');
        % the scaling term in the denominator is sqrt of projected noise, as per eqn. 2.67 of Sekihara & Nagarajan 2008 (S&N)
        if fixedori
          filt = pinv(sqrt(noise * lf' * invC_squared * lf)) * lf' *invC; % based on S&N eqn. 4.08
        else
          ft_error('vector version of nai weight normalization is not implemented');
        end
      case 'unitnoisegain'
        % filt*filt' = I
        % Unit-noise gain minimum variance (aka Borgiotti-Kaplan) beamformer
        % below equation is equivalent to following:
        % filt = pinv(lf' * invC * lf) * lf' * invC;
        % filt = filt/sqrt(filt*filt');
        if fixedori
          filt = pinv(sqrt(lf' * invC_squared * lf)) * lf' *invC; % S&N eqn. 4.15
        else
          % compute the matrix that is used for scaling of the filter's rows, as per eqn. 4.83
          denom = pinv(lf' * invC * lf);
          gamma = denom * (lf' * invC_squared * lf) * denom;
          
          % compute the spatial filter, as per eqn. 4.85
          filt = diag(1./sqrt(diag(gamma))) * denom * lf' * invC;
        end
      case 'arraygain'
        % filt*lf = ||lf||, applies to scalar leadfield, and to one of the possibilities of the vector version, eqn. 4.75
        lfn  = lf./norm(lf);
        filt = pinv(lfn' * invC * lfn) * lfn' * invC; % S&N eqn. 4.09 (scalar version), and eqn. 4.75 (vector version)
        
      case {'unitgain' 'no'}
        % this is the 'standard' unit gain constraint spatial filter: filt*lf=I, applies both to vector and scalar leadfields
        filt = pinv(lf' * invC * lf) * lf' * invC; % Gross eqn. 3 & van Veen eqn. 23, use PINV/SVD to cover rank deficient leadfield
        
      otherwise
    end
  end
  
  cfilt = filt * C * ctranspose(filt); % Gross eqn. 4 and 5
  if powlambda1
    if size(cfilt,1) == 1
      % only 1 orientation, no need to do svd
      estimate.pow(i,1) = real(cfilt);
    else
      estimate.pow(i,1) = lambda1(cfilt);   % compute the power at the dipole location, Gross eqn. 8
    end
  elseif powtrace
    estimate.pow(i,1) = real(trace(cfilt)); % compute the power at the dipole location
  end
  
  if strcmp(submethod, 'dics_refchan')
    % compute the cross-spectrum between the dipole and the reference channel
    cfilt = filt*Cr; % Gross eqn. 6, replaces cfilt computed above
    if powlambda1
      coh = lambda1(cfilt)^2 / (estimate.pow(i,1) * Pr); % Gross eqn. 9
    elseif powtrace
      coh = norm(cfilt)^2 / (estimate.pow(i,1) * Pr);
    end
    estimate.coh(i,1) = coh;
  elseif strcmp(submethod, 'dics_refdip')
    % compute the cross-spectrum between the reference dipole and the dipole
    cfilt = filt1 * C * ctranspose(filt); % Gross eqn. 6, replaces cfilt computed above
    if powlambda1
      coh = lambda1(cfilt)^2 / (estimate.pow(i,1) * Pr); % Gross eqn. 9
    elseif powtrace
      coh = norm(cfilt)^2 / (estimate.pow(i,1) * Pr);
    end
    estimate.coh(i,1) = coh;
  end
  
  if keepcsd
    % keep the dipole's csd or the dipole-ref cross-spectrum in the output
    estimate.csd{i,1} = cfilt;
  end
  
  if projectnoise
    if powlambda1
      estimate.noise(i,1) = noise * lambda1(filt * ctranspose(filt));
    elseif powtrace
      estimate.noise(i,1) = noise * real(trace(filt * ctranspose(filt)));
    end
    if keepcsd
      estimate.noisecsd{i,1} = noise * filt * ctranspose(filt);
    end
  end
  
  if keepfilter
    if hassubspace
      estimate.filter{i,1} = filt*sourcemodel.subspace{i};
    elseif ~isempty(subspace)
      estimate.filter{i,1} = filt*subspace; % FIXME should this be subspace, or pinv(subspace)?
    else
      estimate.filter{i,1} = filt;
    end
  end
  
  if keepleadfield
    if ~isempty(subspace)
      estimate.leadfield{i,1} = lforig;
    else
      estimate.leadfield{i,1} = lf;
    end
  end
  
end % for each dipole position
ft_progress('close');

% reassign the estimated values over the inside and outside grid positions
estimate.inside  = originside;
estimate.pos     = origpos;

fnames_cell   = {'leadfield' 'filter' 'ori' 'csd' 'noisecsd' 'subspace'};
for k = 1:numel(fnames_cell)
  if isfield(estimate, fnames_cell{k})
    estimate.(fnames_cell{k})( originside) = estimate.(fnames_cell{k});
    estimate.(fnames_cell{k})(~originside) = {[]};
  end
end

fnames_scalar = {'pow' 'noise' 'eta' 'coh'};
for k = 1:numel(fnames_scalar)
  if isfield(estimate, fnames_scalar{k})
    estimate.(fnames_scalar{k})( originside) = estimate.(fnames_scalar{k});
    estimate.(fnames_scalar{k})(~originside) = nan;
  end
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
