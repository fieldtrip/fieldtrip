function [dipout] = harmony(dip, grad, headmodel, dat, varargin)

% HARMONY computes a linear estimate of the current in a distributed source model
% using a mesh harmonic based low-pass filter.
%
% Use as
%   [dipout] = minimumnormestimate(dip, grad, headmodel, dat, ...)
%
% Optional input arguments should come in key-value pairs and can include
%  'noisecov'             = Nchan x Nchan matrix with noise covariance
%  'noiselambda'          = scalar value, regularisation parameter for the noise covariance matrix (default=0)
%  'filter_order'         = scalar, order of the mesh Butterwirth filter
%  'filter_bs'            = scalar, stop-band of the mesh Butterworth filter
%  'number_harmonics'     = Integer, number of mesh harmonics used (can be empty, the default will then be identity)
%  'lambda'               = scalar, regularisation parameter (can be empty, it will then be estimated from snr)
%  'snr'                  = scalar, signal to noise ratio
%  'reducerank'           = reduce the leadfield rank, can be 'no' or a number (e.g. 2)
%  'normalize'            = normalize the leadfield
%  'normalizeparam'       = parameter for depth normalization (default = 0.5)
%  'keepfilter'           = 'no' or 'yes', keep the spatial filter in the output
%  'prewhiten'            = 'no' or 'yes', prewhiten the leadfield matrix with the noise covariance matrix C
%  'scalesourcecov'       = 'no' or 'yes', scale the source covariance matrix R such that trace(leadfield*R*leadfield')/trace(C)=1
%  'connected_components' = number of connected components of the source mesh (1 or 2)
%
% Note that leadfield normalization (depth regularisation) should be done by scaling
% the leadfields outside this function, e.g. in prepare_leadfield. Note also that
% with precomputed leadfields the normalization parameters will not have an effect.
%
% This implements
%   Petrov Y (2012) Harmony: EEG/MEG Linear Inverse Source Reconstruction in the
%   Anatomical Basis of Spherical Harmonics. PLoS ONE 7(10): e44439.
%   doi:10.1371/journal.pone.0044439

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
% $Id: harmony.m 10972 2015-12-08 13:09:48Z lucamb $

%% get the optional inputs for the MNE method according to Dale et al 2000, and Liu et al. 2002
noisecov       = ft_getopt(varargin, 'noisecov');
lambda         = ft_getopt(varargin, 'lambda');  % can be empty, it will then be estimated based on SNR
noiselambda    = ft_getopt(varargin, 'noiselambda', []);
snr            = ft_getopt(varargin, 'snr');     % is used to estimate lambda if lambda is not specified

% these settings pertain to the forward model, the defaults are set in compute_leadfield
reducerank     = ft_getopt(varargin, 'reducerank');
normalize      = ft_getopt(varargin, 'normalize');
normalizeparam = ft_getopt(varargin, 'normalizeparam', 0.5);
keepfilter     = istrue(ft_getopt(varargin, 'keepfilter', false));
dowhiten       = istrue(ft_getopt(varargin, 'prewhiten',  false));
doscale        = istrue(ft_getopt(varargin, 'scalesourcecov', false));
hasleadfield   = isfield(dip, 'leadfield');
hasfilter      = isfield(dip, 'filter');

if isempty(lambda) && isempty(snr) && ~isfield(dip, 'filter')
  ft_error('either lambda or snr should be specified');
elseif ~isempty(lambda) && ~isempty(snr)
  ft_error('either lambda or snr should be specified, not both');
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

% select only the dipole positions inside the brain for scanning
dip.pos    = dip.pos(originside,:);
dip.inside = true(size(dip.pos,1),1);
if isfield(dip, 'mom'),
  dip.mom = dip.mom(:,originside);
end
if isfield(dip, 'leadfield')
  hasleadfield = 1;
  dip.leadfield = dip.leadfield(originside);
end
if isfield(dip, 'filter')
  hasfilter = 1;
  dip.filter = dip.filter(originside);
end

%if ~isempty(snr) && doscale
%  ft_error('scaling of the source covariance in combination with a specified snr parameter is not allowed');
%end

% compute leadfield
if hasfilter
  % it does not matter whether the leadfield is there or not, it will not be used
  fprintf('using pre-computed spatial filter: some of the specified options will not have an effect\n');
elseif hasleadfield
  % using the computed leadfields
  fprintf('using pre-computed leadfields: some of the specified options will not have an effect\n');
else
  fprintf('computing forward model\n');
  if isfield(dip, 'mom')
    for i=size(dip.pos,1)
      % compute the leadfield for a fixed dipole orientation
      dip.leadfield{i} = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam) * dip.mom(:,i);
    end
  else
      
    for i=1:size(dip.pos,1)
      % compute the leadfield
      dip.leadfield{i} = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam);
    end
  end
end

%% harmony Parameters
num_dip      = size(dip.leadfield{1}, 2);
num_pos      = size(dip.leadfield, 2);
num_comp     = ft_getopt(varargin, 'connected_components');
num_harm     = ft_getopt(varargin, 'number_harmonics');
filt_ord     = ft_getopt(varargin, 'filter_order');
lc           = ft_getopt(varargin, 'filter_bs');

%% Create the leadfield matrix

% count the number of channels and leadfield components
  Nchan   = size(dip.leadfield{1},1);
  Nsource = 0;
  for i=1:size(dip.pos,1)
    Nsource = Nsource + size(dip.leadfield{i}, 2);
  end
  
  % concatenate the leadfield components of all sources into one large matrix
  lf = zeros(Nchan, Nsource);
  n = 1;
  for i=1:size(dip.pos,1)
    cbeg = n;
    cend = n + num_dip - 1;
    lf(:,cbeg:cend) = dip.leadfield{i};
    n = n + size(dip.leadfield{i}, 2);
  end

%% Create brain harmonics matrix

% Mesh harmonics
[dum,H,d] = ct_mesh_spectrum(dip,num_harm,num_comp);

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
  ft_error('This function only supports meshes with 1 or 2 connected components')
end

if num_dip == 1
  H_tot = H_blk;
elseif num_dip == 2
  H_tot = blkdiag(H_blk,H_blk);
elseif num_dip == 3
  H_tot = blkdiag(H_blk,H_blk,H_blk);
else
  ft_error('This function only supports source model with 1,2 or 3 dimensional source values')
end



%% Create harmonic leadifield
lf_h = zeros(size(lf,1),num_harm,num_comp,num_dip); %leadifield in the harmonic domain

for k = 1:num_comp
    for j = 1:num_dip
        ind = (j-1) + 1:num_dip:num_dip*num_pos;
        lf_h(:,:,k,j) = lf(:,ind)*H{k};
    end
end

lf_h = reshape(lf_h,size(lf,1),num_comp*num_harm*num_dip);

%% Create harmonic source covariance matrix

q = 1./sqrt(1 + (l/lc).^(2*filt_ord));
sourcecov = diag(q);

%% Compute the Harmony spatial filters

  
% count the number of channels and leadfield components
Nchan   = size(dip.leadfield{1},1);
Nsource = 0;
for i=1:size(dip.pos,1)
Nsource = Nsource + size(dip.leadfield{i}, 2);
end

% Take the rieal part of the noise cross-spectral density matrix
if isreal(noisecov) == 0
  fprintf('Taking the real part of the noise cross-spectral density matrix\n');
  noisecov = real(noisecov);
end

% compute the inverse of the forward model, this is where prior information
% on source and noise covariance would be useful
if isempty(noisecov)
% use an unregularised minimum norm solution, i.e. using the Moore-Penrose pseudoinverse
warning('computing a unregularised minimum norm solution. This typically does not work due to numerical accuracy problems');
w = pinv(lf_h);
elseif ~isempty(noisecov)
fprintf('computing the solution where the noise covariance is used for regularisation\n');
% the noise covariance has been given and can be used to regularise the solution
if isempty(sourcecov)
  sourcecov = speye(Nsource);
end
% rename some variables for consistency with the publications
A = lf_h;
R = sourcecov;
C = noisecov;

if dowhiten,
  fprintf('prewhitening the leadfields using the noise covariance\n');

  % compute the prewhitening matrix
  if ~isempty(noiselambda)
    fprintf('using a regularized noise covariance matrix\n');
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
  fprintf('scaling the source covariance\n');
  scale = trace(A*(R*A'))/trace(C);
  R     = R./scale;
end

if ~isempty(snr)
  % the regularisation parameter can be estimated from the noise covariance,
  % see equation 6 in Lin et al. 2004
  lambda = trace(A * R * A')/(trace(C)*snr^2);
end

if dowhiten,
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
  %% equation 5 from Lin et al 2004 (this implements Dale et al 2000, and Liu et al. 2002)
  denom = (A*R*A'+(lambda^2)*C);
  if cond(denom)<1e12
    w = R * A' / denom;
  else
    fprintf('taking pseudo-inverse due to large condition number\n');
    w = R * A' * pinv(denom);
  end
end

end % if empty noisecov
  
%% for each of the timebins, estimate the source strength
w_x = H_tot*w; % The Harm matrix project the Harmony level filters back to the source locations

if isreal(dat) == 1
    fprintf('The input are sensors time-series: Computing the dipole moments\n')
    mom = w_x * dat; 
    mom_ind = 1;
elseif size(dat,1)==size(dat,2)&&sum(sum(dat-dat'))<10^-5*sum(diag(dat))
    fprintf('The input is a sensor level cross-spectral density: Computing source level power\n')
    pow_tot = real(sum((w_x*dat).*w_x,2));
    pow = 0;
    for j = 1:num_dip
        pow = pow + pow_tot((1 + (j-1)*end/num_dip):(j*end/num_dip));
        mom_ind = 0;
    end
else
    fprintf('The input is are sensor level Fourier-coefficients: Computing source level Fourier coefficients\n')
    mom = w_x*dat; % The Harm matrix project the Harmonic activity back to the source locations
    mom_ind = 1;
end

% assign the estimated source strength to each dipole
if mom_ind == 1
n = 1;
for i=1:size(dip.pos,1)
  cbeg = n;
  cend = n + size(dip.leadfield{i}, 2) - 1;
  dipout.mom{i} = mom(cbeg:cend,:);
  n = n + size(dip.leadfield{i}, 2);
end
end

%% compute power (over the three orientations) at each location and for each time

if mom_ind == 1
  dipout.pow = nan(size(dip.pos,1), size(dat,2));
  for i=1:size(dip.pos,1)
    dipout.pow(i,:) = sum(abs(dipout.mom{i}).^2, 1);
  end
else
    dipout.pow = pow;
end

% deal with keepfilter option
if keepfilter && ~hasfilter
  % spatial filters have been computed, store them in the output
  % re-assign spatial filter to conventional 1 cell per dipole location
  n = 1;
  for i=1:size(dip.pos,1)
    cbeg = n;
    cend = n + size(dip.leadfield{i}, 2) - 1;
    dipout.filter{i} = w(cbeg:cend,:);
    n    = n + size(dip.leadfield{i}, 2);
  end
elseif keepfilter
  dipout.filter = dip.filter;
end

% deal with noise covariance
% if ~isempty(noisecov) && ~hasfilter
%   % compute estimate of the projected noise
%   n = 1;
%   for i=1:size(dip.pos,1)
%     cbeg = n;
%     cend = n + size(dip.leadfield{i}, 2) - 1;
%     dipout.noisecov{i} = w(cbeg:cend,:)*noisecov*w(cbeg:cend,:)';
%     n    = n + size(dip.leadfield{i}, 2);
%   end
% elseif ~isempty(noisecov)
%   % compute estimate of the projected noise
%   for i=1:size(dip.pos,1)
%     dipout.noisecov{i} = dip.filter{i}*noisecov*dip.filter{i}';
%   end
%end

% wrap it all up, prepare the complete output
dipout.inside  = originside;
dipout.pos     = origpos;

if isfield(dipout, 'mom')
  dipout.mom( originside) = dipout.mom;
  dipout.mom(~originside) = {[]};
end
if isfield(dipout, 'pow')
  dipout.pow( originside,:) = dipout.pow;
  dipout.pow(~originside,:) = nan;
end
if isfield(dipout, 'noisecov')
  dipout.noisecov( originside) = dipout.noisecov;
  dipout.noisecov(~originside) = {[]};
end
if isfield(dipout, 'filter')
  dipout.filter( originside) = dipout.filter;
  dipout.filter(~originside) = {[]};
end






















