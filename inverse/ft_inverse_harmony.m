function [estimate] = ft_inverse_harmony(sourcemodel, sens, headmodel, dat, varargin)

% FT_INVERSE_HARMONY computes a linear estimate of the current in a distributed
% source model using a mesh harmonic based low-pass filter.
%
% Use as
%   [estimate] = ft_inverse_harmony(sourcemodel, sens, headmodel, dat, ...)
% where
%   sourcemodel is the input source model, see FT_PREPARE_SOURCEMODEL
%   sens        is the gradiometer or electrode definition, see FT_DATATYPE_SENS
%   headmodel   is the volume conductor definition, see FT_PREPARE_HEADMODEL
%   dat         is the data matrix with the ERP or ERF
% and
%   estimate    contains the estimated source parameters
%
% Additional input arguments should be specified as key-value pairs and can include
%   'noisecov'             = Nchan x Nchan matrix with noise covariance
%   'noiselambda'          = scalar value, regularisation parameter for the noise covariance matrix (default=0)
%   'filter_order'         = scalar, order of the mesh Butterwirth filter
%   'filter_bs'            = scalar, stop-band of the mesh Butterworth filter
%   'number_harmonics'     = Integer, number of mesh harmonics used (can be empty, the default will then be identity)
%   'lambda'               = scalar, regularisation parameter (can be empty, it will then be estimated from snr)
%   'snr'                  = scalar, signal to noise ratio
%   'scalesourcecov'       = 'no' or 'yes', scale the source covariance matrix R such that trace(leadfield*R*leadfield')/trace(C)=1
%   'connected_components' = number of connected components of the source mesh (1 or 2)
%   'prewhiten'            = 'no' or 'yes', prewhiten the leadfield matrix with the noise covariance matrix C
%
% These options influence the forward computation of the leadfield
%   'reducerank'      = 'no' or number  (default = 3 for EEG, 2 for MEG)
%   'backproject'     = 'yes' or 'no', in the case of a rank reduction this parameter determines whether the result will be backprojected onto the original subspace (default = 'yes')
%   'normalize'       = 'no', 'yes' or 'column' (default = 'no')
%   'normalizeparam'  = parameter for depth normalization (default = 0.5)
%   'weight'          = number or Nx1 vector, weight for each dipole position to compensate for the size of the corresponding patch (default = 1)
%
% This implements
% - Petrov Y (2012) Harmony: EEG/MEG Linear Inverse Source Reconstruction in the
%   Anatomical Basis of Spherical Harmonics. PLoS ONE 7(10): e44439.
%   doi:10.1371/journal.pone.0044439
%
% See also FT_SOURCEANALYSIS, FT_PREPARE_HEADMODEL, FT_PREPARE_SOURCEMODEL

% TODO implement the following options
% - keepleadfield

% Copyright (C) 2015, Luca Ambrogio
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

if mod(nargin-4,2)
  % the first 4 arguments are fixed, the other arguments should come in pairs
  ft_error('invalid number of optional arguments');
end

% get the optional input arguments, or use defaults
noisecov       = ft_getopt(varargin, 'noisecov');
lambda         = ft_getopt(varargin, 'lambda');  % can be empty, it will then be estimated based on SNR
noiselambda    = ft_getopt(varargin, 'noiselambda', []);
snr            = ft_getopt(varargin, 'snr');     % is used to estimate lambda if lambda is not specified
keepfilter     = istrue(ft_getopt(varargin, 'keepfilter', false));
dowhiten       = istrue(ft_getopt(varargin, 'prewhiten',  false));
doscale        = istrue(ft_getopt(varargin, 'scalesourcecov', false));

% construct the low-level options for the leadfield computation as key-value pairs, these are passed to FT_COMPUTE_LEADFIELD
leadfieldopt = {};
leadfieldopt = ft_setopt(leadfieldopt, 'reducerank',     ft_getopt(varargin, 'reducerank'));
leadfieldopt = ft_setopt(leadfieldopt, 'backproject',    ft_getopt(varargin, 'backproject'));
leadfieldopt = ft_setopt(leadfieldopt, 'normalize',      ft_getopt(varargin, 'normalize'));
leadfieldopt = ft_setopt(leadfieldopt, 'normalizeparam', ft_getopt(varargin, 'normalizeparam'));
leadfieldopt = ft_setopt(leadfieldopt, 'weight',         ft_getopt(varargin, 'weight'));

% flags to avoid calling isfield repeatedly in the loop over grid positions (saves a lot of time)
hasmom        = isfield(sourcemodel, 'mom');
hasleadfield  = isfield(sourcemodel, 'leadfield');
hasfilter     = isfield(sourcemodel, 'filter');

if isempty(lambda) && isempty(snr) && ~hasfilter
  ft_error('either lambda or snr should be specified');
elseif ~isempty(lambda) && ~isempty(snr)
  ft_error('either lambda or snr should be specified, not both');
end

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
  if hasmom
    for i=size(sourcemodel.pos,1)
      % compute the leadfield for a fixed dipole orientation
      sourcemodel.leadfield{i} = ft_compute_leadfield(sourcemodel.pos(i,:), sens, headmodel, leadfieldopt{:}) * sourcemodel.mom(:,i);
    end
  else
    for i=1:size(sourcemodel.pos,1)
      % compute the leadfield
      sourcemodel.leadfield{i} = ft_compute_leadfield(sourcemodel.pos(i,:), sens, headmodel, leadfieldopt{:});
    end
  end
end

% Harmony parameters
num_dip      = size(sourcemodel.leadfield{1}, 2);
num_pos      = size(sourcemodel.leadfield, 2);
num_comp     = ft_getopt(varargin, 'connected_components', 1);
num_harm     = ft_getopt(varargin, 'number_harmonics',     5);
filt_ord     = ft_getopt(varargin, 'filter_order',         5);
lc           = ft_getopt(varargin, 'filter_bs',            5);

% Create the leadfield matrix

% count the number of channels and leadfield components
Nchan   = size(sourcemodel.leadfield{1},1);
Nsource = 0;
for i=1:size(sourcemodel.pos,1)
  Nsource = Nsource + size(sourcemodel.leadfield{i}, 2);
end

% concatenate the leadfield components of all sources into one large matrix
lf = zeros(Nchan, Nsource);
n = 1;
for i=1:size(sourcemodel.pos,1)
  cbeg = n;
  cend = n + num_dip - 1;
  lf(:,cbeg:cend) = sourcemodel.leadfield{i};
  n = n + size(sourcemodel.leadfield{i}, 2);
end

% Create brain harmonics matrix

% Mesh harmonics
[dum,H,d] = mesh_spectrum(sourcemodel,num_harm,num_comp);

% Elimination of pathological zeroth frequencies and concatenating
% harmonics numbers
l = [];
for n = 1:num_comp
  d{n}(1) = 0;
  H{n}(:,1) = 1/sqrt(numel(H{n}(:,1)));
  d_cat = [];
  for m = 1:num_dip
    d_cat = [d_cat; d{n}];
  end
  l = [l; -d_cat];
end

% Total harmonics matrix
if num_comp == 1
  H_blk = H{1};
elseif num_comp == 2
  H_blk = blkdiag(H{1},H{2});
else
  ft_error('this function only supports meshes with 1 or 2 connected components')
end

if num_dip == 1
  H_tot = H_blk;
elseif num_dip == 2
  H_tot = blkdiag(H_blk,H_blk);
elseif num_dip == 3
  H_tot = blkdiag(H_blk,H_blk,H_blk);
else
  ft_error('this function only supports source model with 1, 2 or 3 dimensional source values')
end

% Create harmonic leadfield
lf_h = zeros(size(lf,1),num_harm,num_comp,num_dip); %leadifield in the harmonic domain

for k = 1:num_comp
  for j = 1:num_dip
    ind = (j-1) + 1:num_dip:num_dip*num_pos;
    lf_h(:,:,k,j) = lf(:,ind)*H{k};
  end
end

lf_h = reshape(lf_h,size(lf,1),num_comp*num_harm*num_dip);

% Create harmonic source covariance matrix
q = 1./sqrt(1 + (l/lc).^(2*filt_ord));
sourcecov = diag(q);

% Compute the Harmony spatial filters

% count the number of channels and leadfield components
Nchan   = size(sourcemodel.leadfield{1},1);
Nsource = 0;
for i=1:size(sourcemodel.pos,1)
  Nsource = Nsource + size(sourcemodel.leadfield{i}, 2);
end

% Take the rieal part of the noise cross-spectral density matrix
if isreal(noisecov)
  ft_info('taking the real part of the noise cross-spectral density matrix\n');
  noisecov = real(noisecov);
end

% compute the inverse of the forward model, this is where prior information
% on source and noise covariance would be useful
if isempty(noisecov)
  % use an unregularised minimum norm solution, i.e. using the Moore-Penrose pseudoinverse
  ft_warning('computing a unregularised minimum norm solution. This typically does not work due to numerical accuracy problems');
  w = pinv(lf_h);
elseif ~isempty(noisecov)
  ft_info('using the noise covariance for regularisation\n');
  % the noise covariance has been given and can be used to regularise the solution
  if isempty(sourcecov)
    sourcecov = speye(Nsource);
  end
  % rename some variables for consistency with the publications
  A = lf_h;
  R = sourcecov;
  C = noisecov;
  
  if dowhiten
    ft_info('prewhitening the leadfields using the noise covariance\n');
    
    % compute the prewhitening matrix
    if ~isempty(noiselambda)
      ft_info('using a regularized noise covariance matrix\n');
      % note: if different channel types are present, one should probably load the diagonal with channel-type specific stuff
      [U,S,V] = svd(C+eye(size(C))*noiselambda);
    else
      [U,S,V] = svd(C);
    end
    
    Tol     = 1e-12;
    diagS   = diag(S);
    sel     = find(diagS>Tol.*diagS(1));
    P       = diag(1./sqrt(diag(S(sel,sel))))*U(:,sel)'; % prewhitening matrix
    A       = P*A; % prewhitened leadfields
    C       = eye(size(P,1)); % prewhitened noise covariance matrix
  end
  
  if doscale
    % estimate sourcecov such that trace(ARA')/trace(C) = 1 (see
    % http://martinos.org/mne/manual/mne.html. In the case of prewhitening
    % C reduces to I (and then lambda^2 ~ 1/SNR); note that in mixed
    % channel type covariance matrices prewhitening should be applied in
    % order for this to make sense (otherwise the diagonal elements of C
    % have different units)
    ft_info('scaling the source covariance\n');
    scale = trace(A*(R*A'))/trace(C);
    R     = R./scale;
  end
  
  if ~isempty(snr)
    % the regularisation parameter can be estimated from the noise covariance,
    % see equation 6 in Lin et al. 2004
    lambda = trace(A * R * A')/(trace(C)*snr^2);
  end
  
  if dowhiten
    % as documented on MNE website, this is replacing the part of the code below, it gives
    % more stable results numerically.
    Rc      = chol(R, 'lower');
    [U,S,V] = svd(A * Rc, 'econ');
    s  = diag(S);
    ss = s ./ (s.^2 + lambda);
    w  = Rc * V * diag(ss) * U';
    
    % unwhiten the filters to bring them back into signal subspace
    w = w*P;
    
  else
    % equation 5 from Lin et al 2004 (this implements Dale et al 2000, and Liu et al. 2002)
    denom = (A*R*A'+(lambda^2)*C);
    if cond(denom)<1e12
      w = R * A' / denom;
    else
      ft_info('taking pseudo-inverse due to large condition number\n');
      w = R * A' * pinv(denom);
    end
  end
  
end % if empty noisecov

% for each of the timebins, estimate the source strength
w_x = H_tot*w; % The Harm matrix project the Harmony level filters back to the source locations

if isreal(dat)
  ft_info('the input consists of time-series: computing the dipole moments\n')
  mom = w_x * dat;
  mom_ind = 1;
elseif size(dat,1)==size(dat,2)&&sum(sum(dat-dat'))<10^-5*sum(diag(dat))
  ft_info('the input consists of a cross-spectral density: computing source-level power\n')
  pow_tot = real(sum((w_x*dat).*w_x,2));
  pow = 0;
  for j = 1:num_dip
    pow = pow + pow_tot((1 + (j-1)*end/num_dip):(j*end/num_dip));
    mom_ind = 0;
  end
else
  ft_info('the input consists of Fourier coefficients: computing source-level Fourier coefficients\n')
  mom = w_x*dat; % The Harm matrix project the Harmonic activity back to the source locations
  mom_ind = 1;
end

% assign the estimated source strength to each dipole
if mom_ind == 1
  n = 1;
  for i=1:size(sourcemodel.pos,1)
    cbeg = n;
    cend = n + size(sourcemodel.leadfield{i}, 2) - 1;
    estimate.mom{i} = mom(cbeg:cend,:);
    n = n + size(sourcemodel.leadfield{i}, 2);
  end
end

% compute power (over the three orientations) at each location and for each time

if mom_ind == 1
  estimate.pow = nan(size(sourcemodel.pos,1), size(dat,2));
  for i=1:size(sourcemodel.pos,1)
    estimate.pow(i,:) = sum(abs(estimate.mom{i}).^2, 1);
  end
else
  estimate.pow = pow;
end

% deal with keepfilter option
if keepfilter && ~hasfilter
  % spatial filters have been computed, store them in the output
  % re-assign spatial filter to conventional 1 cell per dipole location
  n = 1;
  for i=1:size(sourcemodel.pos,1)
    cbeg = n;
    cend = n + size(sourcemodel.leadfield{i}, 2) - 1;
    estimate.filter{i} = w(cbeg:cend,:);
    n    = n + size(sourcemodel.leadfield{i}, 2);
  end
elseif keepfilter
  estimate.filter = sourcemodel.filter;
end

% deal with noise covariance
% if ~isempty(noisecov) && ~hasfilter
%   % compute estimate of the projected noise
%   n = 1;
%   for i=1:size(sourcemodel.pos,1)
%     cbeg = n;
%     cend = n + size(sourcemodel.leadfield{i}, 2) - 1;
%     estimate.noisecov{i} = w(cbeg:cend,:)*noisecov*w(cbeg:cend,:)';
%     n    = n + size(sourcemodel.leadfield{i}, 2);
%   end
% elseif ~isempty(noisecov)
%   % compute estimate of the projected noise
%   for i=1:size(sourcemodel.pos,1)
%     estimate.noisecov{i} = sourcemodel.filter{i}*noisecov*sourcemodel.filter{i}';
%   end
%end

% reassign the estimated values over the inside and outside grid positions
estimate.inside  = originside;
estimate.pos     = origpos;
if isfield(estimate, 'pow')
  estimate.pow( originside,:) = estimate.pow;
  estimate.pow(~originside,:) = nan;
end
if isfield(estimate, 'mom')
  estimate.mom( originside) = estimate.mom;
  estimate.mom(~originside) = {[]};
end
if isfield(estimate, 'noisecov')
  estimate.noisecov( originside) = estimate.noisecov;
  estimate.noisecov(~originside) = {[]};
end
if isfield(estimate, 'filter')
  estimate.filter( originside) = estimate.filter;
  estimate.filter(~originside) = {[]};
end
