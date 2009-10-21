function [dipout] = minimumnormestimate(dip, grad, vol, dat, varargin);

% MINIMUMNORMESTIMATE computes a linear estimate of the current 
% in a distributed source model  
% 
% Use as 
%   [dipout] = minimumnormestimate(dip, grad, vol, dat, ...)
%
% Optional input arguments should come in key-value pairs and can include
%  'noisecov'         = Nchan x Nchan matrix with noise covariance
%  'sourcecov'        = Nsource x Nsource matrix with source covariance (can be empty, the default will then be identity)
%  'lambda'           = scalar, regularisation parameter (can be empty, it will then be estimated from snr)
%  'snr'              = scalar, signal to noise ratio
%  'reducerank'       = reduce the leadfield rank, can be 'no' or a number (e.g. 2)
%  'normalize'        = normalize the leadfield
%  'normalizeparam'   = parameter for depth normalization (default = 0.5)
%
% Note that leadfield normalization (depth regularisation) should be
% done by scaling the leadfields outside this function, e.g. in
% prepare_leadfield.
%
% This implements
% * Dale AM, Liu AK, Fischl B, Buckner RL, Belliveau JW, Lewine JD,
%   Halgren E (2000): Dynamic statistical parametric mapping: combining
%   fMRI and MEG to produce high-resolution spatiotemporal maps of
%   cortical activity. Neuron 26:55-67.
% * Arthur K. Liu, Anders M. Dale, and John W. Belliveau  (2002): Monte
%   Carlo Simulation Studies of EEG and MEG Localization Accuracy.
%   Human Brain Mapping 16:47-62.
% * Fa-Hsuan Lin, Thomas Witzel, Matti S. Hamalainen, Anders M. Dale,
%   John W. Belliveau, and Steven M. Stufflebeam (2004): Spectral
%   spatiotemporal imaging of cortical oscillations and interactions
%   in the human brain.  NeuroImage 23:582-595.

% TODO implement the following options
% - keepleadfield
% - keepfilter
% - keepinverse (i.e. equivalent to keepfilter)

% Copyright (C) 2004-2008, Robert Oostenveld
% 
% $Log: minimumnormestimate.m,v $
% Revision 1.7  2008/04/21 14:50:10  roboos
% added full literature references to documentation
%
% Revision 1.6  2008/03/18 13:18:06  roboos
% added optional settings for compute_leadfield
% fixed typo in equation (would have cased Matlab error)
% updated documentation
%
% Revision 1.5  2007/08/07 08:48:21  roboos
% implemented regularized MNE according to the method described by Dale and Liu
%
% Revision 1.4  2006/06/22 12:05:54  roboos
% changed a comment and a whitespace
%
% Revision 1.3  2006/05/10 08:17:30  roboos
% some changes in comments and documentation
%
% Revision 1.2  2005/10/05 06:32:39  roboos
% minor change in user feedback
%
% Revision 1.1  2004/09/21 13:26:49  roboos
% initial implementation
%

% ensure that these are row-vectors
dip.inside = dip.inside(:)';
dip.outside = dip.outside(:)';

% get the optional inputs for the MNE method according to Dale et al 2000, and Liu et al. 2002
noisecov       = keyval('noisecov',       varargin);
sourcecov      = keyval('sourcecov',      varargin);
lambda         = keyval('lambda',         varargin);  % can be empty, it will then be estimated based on SNR
snr            = keyval('snr',            varargin);  % is used to estimate lambda if lambda is not specified
% these settings pertain to the forward model, the defaults are set in compute_leadfield
reducerank     = keyval('reducerank',     varargin);
normalize      = keyval('normalize',      varargin);
normalizeparam = keyval('normalizeparam', varargin);

if ~isfield(dip, 'leadfield')
  fprintf('computing forward model\n');
  if isfield(dip, 'mom')
    for i=dip.inside
      % compute the leadfield for a fixed dipole orientation
      dip.leadfield{i} = compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam) * dip.mom(:,i);
    end
  else
    for i=dip.inside
      % compute the leadfield
      dip.leadfield{i} = compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam);
    end
  end
  for i=dip.outside
    dip.leadfield{i} = nan;
  end
else
  fprintf('using specified forward model\n');
end

Nchan = length(grad.label);

% count the number of leadfield components for each source
Nsource = 0;
for i=dip.inside
  Nsource = Nsource + size(dip.leadfield{i}, 2);
end

% concatenate the leadfield components of all sources into one large matrix
lf = zeros(Nchan, Nsource);
n = 1;
for i=dip.inside
  cbeg = n;
  cend = n + size(dip.leadfield{i}, 2) - 1;
  lf(:,cbeg:cend) = dip.leadfield{i};
  n = n + size(dip.leadfield{i}, 2);
end

fprintf('computing MNE source reconstruction, this may take some time...\n');
% compute the inverse o the forward model, this is where prior information
% on source and noise covariance would be usefull
if isempty(noisecov)
  % use an unregularised minimum norm solution, i.e. using the Moore-Penrose pseudoinverse
  w = pinv(lf);
else 
  % the noise covariance has been given and can be used to regularise the solution
  if isempty(sourcecov)
    sourcecov = eye(Nsource);
  end
  % rename some variables for consistency with the publications
  A = lf;
  R = sourcecov;
  C = noisecov;
  % the regularisation parameter can be estimated from the noise covariance, see equation 6 in Lin et al. 2004
  if isempty(lambda)
    lambda = trace(A * R * A')/(trace(C)*snr^2);
  end
  % equation 5 from Lin et al 2004 (this implements Dale et al 2000, and Liu et al. 2002)
  w = R * A' * inv( A * R * A' + (lambda^2) * C);
end

% for each of the timebins, estimate the source strength
mom = w * dat;

% re-assign the estimated source strength over the inside and outside dipoles
n = 1;
for i=dip.inside
  cbeg = n;
  cend = n + size(dip.leadfield{i}, 2) - 1;
  dipout.mom{i} = mom(cbeg:cend,:);
  n = n + size(dip.leadfield{i}, 2);
end
dipout.mom(dip.outside) = {nan};

% for convenience also compute power at each location
for i=dip.inside
  dipout.pow(i,:) = sum(dipout.mom{i}.^2, 1);
end
dipout.pow(dip.outside,:) = nan;

% add other descriptive information to the output source model
dipout.pos     = dip.pos;
dipout.inside  = dip.inside;
dipout.outside = dip.outside;

