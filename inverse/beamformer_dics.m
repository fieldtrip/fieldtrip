function [dipout] = beamformer_dics(dip, grad, headmodel, dat, C, varargin)

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
%   dipin.mom   dipole orientation (optional), 3 x Npositions,
% and can additionally contain things like a precomputed filter.
%
% Additional options should be specified in key-value pairs and can be
%  'Pr'               = power of the external reference channel
%  'Cr'               = cross spectral density between all data channels and the external reference channel
%  'refdip'           = location of dipole with which coherence is computed
%  'powmethod'        = can be 'trace' or 'lambda1'
%  'feedback'         = give ft_progress indication, can be 'text', 'gui' or 'none'
%  'fixedori'         = use fixed or free orientation,                 can be 'yes' or 'no'
%  'projectnoise'     = project noise estimate through filter,         can be 'yes' or 'no'
%  'realfilter'       = construct a real-valued filter,                can be 'yes' or 'no'
%  'keepfilter'       = remember the beamformer filter,                can be 'yes' or 'no'
%  'keepleadfield'    = remember the forward computation,              can be 'yes' or 'no'
%  'keepcsd'          = remember the estimated cross-spectral density, can be 'yes' or 'no'
%  'weightnorm'       = normalize the beamformer weights,              can be 'no', 'unitnoisegain' or 'nai'
%
% These options influence the forward computation of the leadfield
%   'reducerank'      = 'no', or number (default = 3 for EEG, 2 for MEG)
%   'backproject'     = 'yes' or 'no',  determines when reducerank is applied whether the 
%                        lower rank leadfield is projected back onto the original linear 
%                        subspace, or not (default = 'yes')
%   'normalize'       = 'yes' or 'no' (default = 'no')
%   'normalizeparam'  = depth normalization parameter (default = 0.5)
%   'weight'          = number or Nx1 vector, weight for each dipole position to compensate 
%                        for the size of the corresponding patch (default = 1)
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
Pr             = ft_getopt(varargin, 'Pr');
Cr             = ft_getopt(varargin, 'Cr');
refdip         = ft_getopt(varargin, 'refdip');
powmethod      = ft_getopt(varargin, 'powmethod');  % the default for this is set below
realfilter     = ft_getopt(varargin, 'realfilter'); % the default for this is set below
subspace       = ft_getopt(varargin, 'subspace');

% these optional settings have defaults
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

% flags to avoid calling isfield repeatedly in the loop over grid positions (saves a lot of time)
hasmom        = isfield(dip, 'mom');
hasleadfield  = isfield(dip, 'leadfield');
hasfilter     = isfield(dip, 'filter');
hassubspace   = isfield(dip, 'subspace');

if hasfilter && (fixedori || ~isequal(weightnorm, 'no')
  ft_warning('with precomputed spatial filters a fixed orientation constraint or weight normalisation options are not applied');
end

% keep the original details on inside and outside positions
originside = dip.inside;
origpos    = dip.pos;

% select only the dipole positions inside the brain for scanning
dip.pos    = dip.pos(originside,:);
dip.inside = true(size(dip.pos,1),1);
if hasmom
  dip.mom = dip.mom(:,originside);
end
if hasleadfield
  ft_info('using precomputed leadfields\n');
  dip.leadfield = dip.leadfield(originside);
end
if hasfilter
  ft_info('using precomputed filters\n');
  dip.filter = dip.filter(originside);
end
if hassubspace
  ft_info('using subspace projection\n');
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
  % estimate the noise level in the covariance matrix by the smallest (non-zero) singular value
  % always needed for the NAI weight normalization case
  noise = svd(C);
  noise = noise(rankC);
  % estimated noise floor is equal to or higher than lambda
  noise = max(noise, tmplambda);
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

% compute the square of invC, which might be needed for unitnoisegain or NAI constraint
invC_squared = invC^2;

% start the scanning with the proper metric
ft_progress('init', feedback, 'scanning grid');
switch submethod
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % dics_power
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'dics_power'
    % only compute power of a dipole at the grid positions
    for i=1:size(dip.pos,1)
      if hasfilter
        % precomputed filter is provided, the leadfield is not needed
        filt = dip.filter{i};

      else

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
          lf = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, leadfieldopt{:}) * dip.mom(:,i);
        else
          % compute the leadfield
          lf = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, leadifeldopt{:});
        end

        if hassubspace
          % do subspace projection of the forward model
          lf = dip.subspace{i} * lf;
          % the cross-spectral density becomes voxel dependent due to the projection
          C    = dip.subspace{i} * C_pre_subspace * dip.subspace{i}';
          if strcmp(realfilter, 'yes')
            invC = ft_inv(dip.subspace{i} * real(C_pre_subspace) * dip.subspace{i}', invopt{:});
          else
            invC = ft_inv(dip.subspace{i} * C_pre_subspace * dip.subspace{i}', invopt{:});
          end
        elseif ~isempty(subspace)
          % do subspace projection of the forward model only
          lforig = lf;
          lf     = subspace * lf;
        
          % according to Kensuke's paper, the eigenspace bf boils down to projecting
          % the 'traditional' filter onto the subspace spanned by the first k eigenvectors
          % [u,s,v] = svd(Cy); filt = ESES*filt; ESES = u(:,1:k)*u(:,1:k)';
          % however, even though it seems that the shape of the filter is identical to
          % the shape it is obtained with the following code, the w*lf=I does not hold.
        end
      
        if fixedori
          % use single dipole orientation
          if hasfilter && size(filt,1) == 1
            % provided precomputed filter already projects onto one
            % orientation, nothing to be done here
          else
            % find out the optimal dipole orientation
            [u, s, v] = svd(real(filt * C * ctranspose(filt)));
            maxpowori = u(:,1);
            eta = s(1,1)./s(2,2);
            
            % and compute the leadfield for that orientation
            lf  = lf * maxpowori;
            dipout.ori{i} = maxpowori;
            dipout.eta(i) = eta;
            if ~isempty(subspace), lforig = lforig * maxpowori; end
            
            % recompute the filter to only use that orientation
            filt = pinv(lf' * invC * lf) * lf' * invC;
          end
        end
        % compute filter
        filt = pinv(lf' * invC * lf) * lf' * invC; % Gross eqn. 3, use pinv/SVD to cover rank deficient leadfield
      end
      
      csd = filt * C * ctranspose(filt);    % Gross eqn. 4 and 5
      if powlambda1
        if size(csd,1) == 1
          % only 1 orientation, no need to do svd
          dipout.pow(i,1) = real(csd);
        else
          dipout.pow(i,1) = lambda1(csd);   % compute the power at the dipole location, Gross eqn. 8
        end
      elseif powtrace
        dipout.pow(i,1) = real(trace(csd)); % compute the power at the dipole location
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
      if hasfilter
        % precomputed filter is provided, the leadfield is not needed
      elseif hasleadfield
        % reuse the leadfield that was previously computed
        lf = dip.leadfield{i};
      elseif hasmom
        % compute the leadfield for a fixed dipole orientation
        lf = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, leadfieldopt{:}) .* dip.mom(i,:)';
      else
        % compute the leadfield
        lf = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, leadfieldopt{:});
      end
      if hassubspace
        % do subspace projection of the forward model
        lforig = lf;
        lf = dip.subspace{i} * lf;
        % the cross-spectral density becomes voxel dependent due to the projection
        C    = dip.subspace{i} * C_pre_subspace * dip.subspace{i}';
        invC = ft_inv(dip.subspace{i} * (C_pre_subspace + lambda * eye(size(C))) * dip.subspace{i}', invopt{:});
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
        filt = pinv(lf' * invC * lf) * lf' * invC;              % Gross eqn. 3, use pinv/SVD to cover rank deficient leadfield
      end
      
      if fixedori
        % use single dipole orientation
        if hasfilter && size(filt,1) == 1
          % provided precomputed filter already projects to one
          % orientation, nothing to be done here
        else
          % find out the optimal dipole orientation
          [u, s, v] = svd(real(filt * C * ctranspose(filt)));
          maxpowori = u(:,1);
          
          % compute the leadfield for that orientation
          lf  = lf * maxpowori;
          dipout.ori{i,1} = maxpowori;
          
          % recompute the filter to only use that orientation
          filt = pinv(lf' * invC * lf) * lf' * invC;
        end
      elseif hasfilter && size(filt,1) == 1
        ft_error('the precomputed filter you provided projects to a single dipole orientation, but you request fixedori=''no''; this is invalid. Either provide a filter with the three orientations retained, or specify fixedori=''yes''.');
      end
      
      if powlambda1
        [pow, ori] = lambda1(filt * C * ctranspose(filt));            % compute the power and orientation at the dipole location, Gross eqn. 4, 5 and 8
      elseif powtrace
        pow = real(trace(filt * C * ctranspose(filt)));               % compute the power at the dipole location
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
      filt1 = pinv(lf1' * invC * lf1) * lf1' * invC;       % use pinv/SVD to cover rank deficient leadfield
    elseif isstruct(refdip) && isfield(refdip, 'pos') % check if only position of refdip is present
      assert(isnumeric(refdip.pos) && numel(refdip.pos)==3);
      lf1 = ft_compute_leadfield(refdip.pos, grad, headmodel, leadfieldopt{:});
      if isfield(refdip,'mom') % check for fixed orientation
        lf1 = lf1.*refdip.mom(:); 
      end 
      filt1 = pinv(lf1' * invC * lf1) * lf1' * invC;       % use pinv/SVD to cover rank deficient leadfield
    else % backwards compatible with previous implementation - only position of refdip is present
      % compute cortio-cortical coherence with a dipole at the reference position
      lf1 = ft_compute_leadfield(refdip, grad, headmodel, leadfieldopt{:});
      % construct the spatial filter for the first (reference) dipole location
      filt1 = pinv(lf1' * invC * lf1) * lf1' * invC;       % use pinv/SVD to cover rank deficient leadfield
    end
    if powlambda1
      Pref = lambda1(filt1 * C * ctranspose(filt1));      % compute the power at the first dipole location, Gross eqn. 8
    elseif powtrace
      Pref = real(trace(filt1 * C * ctranspose(filt1)));  % compute the power at the first dipole location
    end
    for i=1:size(dip.pos,1)
      if hasfilter
        % precomputed filter is provided, the leadfield is not needed
      elseif hasleadfield
        % reuse the leadfield that was previously computed
        lf2 = dip.leadfield{i};
      elseif hasmom
        % compute the leadfield for a fixed dipole orientation
        lf2 = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, leadfieldopt{:}) .* dip.mom(i,:)';
      else
        % compute the leadfield
        lf2 = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, leadfieldopt{:});
      end
      if hasfilter
        % use the provided filter
        filt2 = dip.filter{i};
      else
        % construct the spatial filter for the second dipole location
        filt2 = pinv(lf2' * invC * lf2) * lf2' * invC;   %  use pinv/SVD to cover rank deficient leadfield
      end
      csd = filt1 * C * ctranspose(filt2);                % compute the cross spectral density between the two dipoles, Gross eqn. 4
      if powlambda1
        pow = lambda1(filt2 * C * ctranspose(filt2));     % compute the power at the second dipole location, Gross eqn. 8
      elseif powtrace
        pow = real(trace(filt2 * C * ctranspose(filt2))); % compute the power at the second dipole location
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

% reassign the scan values over the inside and outside grid positions
dipout.inside  = originside;
dipout.pos     = origpos;

fnames_cell   = {'leadfield' 'filter' 'ori' 'csd' 'noisecsd' 'subspace'};
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

