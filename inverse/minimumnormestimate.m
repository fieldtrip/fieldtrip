function [dipout] = minimumnormestimate(dip, grad, vol, dat, varargin)

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
%  'keepfilter'       = 'no' or 'yes', keep the spatial filter in the
%                       output
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

% Copyright (C) 2004-2008, Robert Oostenveld
% 
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% ensure that these are row-vectors
dip.inside = dip.inside(:)';
dip.outside = dip.outside(:)';

% get the optional inputs for the MNE method according to Dale et al 2000, and Liu et al. 2002
noisecov       = ft_getopt(varargin, 'noisecov');
sourcecov      = ft_getopt(varargin, 'sourcecov');
lambda         = ft_getopt(varargin, 'lambda');  % can be empty, it will then be estimated based on SNR
snr            = ft_getopt(varargin, 'snr');  % is used to estimate lambda if lambda is not specified

% these settings pertain to the forward model, the defaults are set in compute_leadfield
reducerank     = ft_getopt(varargin, 'reducerank');
normalize      = ft_getopt(varargin, 'normalize');
normalizeparam = ft_getopt(varargin, 'normalizeparam');
keepfilter     = istrue(ft_getopt(varargin, 'keepfilter', false));
dowhiten       = istrue(ft_getopt(varargin, 'prewhiten',  false));
doscale        = istrue(ft_getopt(varargin, 'scalesourcecov', false));

if isempty(lambda) && isempty(snr) && ~isfield(dip, 'filter')
  error('either lambda or snr should be specified');
elseif ~isempty(lambda) && ~isempty(snr)
  error('either lambda or snr should be specified, not both');
end

if ~isempty(snr) && doscale
  error('scaling of the source covariance in combination with a specified snr parameter is not allowed');
end

if ~isfield(dip, 'filter')
  
  % compute the leadfields if needed
  if ~isfield(dip, 'leadfield')
    fprintf('computing forward model\n');
    if isfield(dip, 'mom')
      for i=dip.inside
        % compute the leadfield for a fixed dipole orientation
        dip.leadfield{i} = ft_compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam) * dip.mom(:,i);
      end
    else
      for i=dip.inside
        % compute the leadfield
        dip.leadfield{i} = ft_compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam);
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
    warning('computing a unregularised minimum norm solution. This typically does not work due to numerical accuracy problems');
    w = pinv(lf);
  elseif ~isempty(noisecov)
    fprintf('computing the solution where the noise covariance is used for regularisation\n'); 
    % the noise covariance has been given and can be used to regularise the solution
    if isempty(sourcecov)
      sourcecov = speye(Nsource);
    end
    % rename some variables for consistency with the publications
    A = lf;
    R = sourcecov;
    C = noisecov;
    
    if dowhiten,
      fprintf('prewhitening the leadfields using the noise covariance\n');
      
      % compute the prewhitening matrix
      [U,S,V] = svd(C);
      Tol     = 1e-12;
      diagS   = diag(S);
      sel     = find(diagS>Tol.*diagS(1));
      
      P = diag(1./sqrt(diag(S(sel,sel))))*U(:,sel)';
      
      A = P*A;
      C = eye(size(P,1));
    end
    
    if doscale
      % estimate sourcecov such that trace(ARA')/trace(C) = 1 (see
      % http://martinos.org/mne/manual/mne.html. In the case of prewhitening
      % C reduces to I (and then lambda^2 ~ 1/SNR); note that in mixed
      % channel type covariance matrices prewhitening should be applied in
      % order for this to make sense (otherwise the diagonal elements of C
      % have different units)
      fprintf('scaling the source covariance\n');
      scale = trace(A*R*A')/trace(C);
      R     = R./scale;
      
    elseif ~isempty(snr)
      % the regularisation parameter can be estimated from the noise covariance,
      % see equation 6 in Lin et al. 2004
      lambda = trace(A * R * A')/(trace(C)*snr^2);
    end
    
    % equation 5 from Lin et al 2004 (this implements Dale et al 2000, and Liu et al. 2002)
    denom = (A*R*A'+(lambda^2)*C);
    if cond(denom)<1e12
      w = R * A' / denom;
    else
      fprintf('taking pseudo-inverse due to large condition number\n');
      w = R * A' * pinv(denom);
    end
  
    % unwhiten the filters to bring them back into signal subspace
    if dowhiten
      w = w*P;
    end
       
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
else
  
  % if the filter has been pre computed
  fprintf('using pre-computed spatial filters\n');
  dipout.mom = cell(size(dip.pos,1),1);
  for i=dip.inside
    dipout.mom{i} = dip.filter{i} * dat;
  end
  dipout.mom(dip.outside) = {nan};
  
end

% for convenience also compute power (over the three orientations) at each location and for each time
dipout.pow = nan( size(dip.pos,1), size(dat,2));
for i=dip.inside
  dipout.pow(i,:) = sum(dipout.mom{i}.^2, 1);
end

% add other descriptive information to the output source model
dipout.pos     = dip.pos;
dipout.inside  = dip.inside;
dipout.outside = dip.outside;

if (keepfilter || ~isempty(noisecov)) && ~isfield(dip, 'filter')
  % re-assign spatial filter to conventional 1 cell per dipole location
  n = 1;
  for i=dip.inside
    cbeg = n;
    cend = n + size(dip.leadfield{i}, 2) - 1;
    
    if keepfilter
      dipout.filter{i} = w(cbeg:cend,:);
    end
    
    if ~isempty(noisecov)
      dipout.noisecov{i} = w(cbeg:cend,:)*noisecov*w(cbeg:cend,:)';
    end
    
    n = n + size(dip.leadfield{i}, 2);
  end
  
  if keepfilter,         dipout.filter(dip.outside)  = {nan}; end
  if ~isempty(noisecov), dipout.noisecov(dip.outside) = {nan}; end
elseif isfield(dip, 'filter')
  dipout.filter = dip.filter;
  if ~isempty(noisecov)
    for i=dip.inside
      dipout.noisecov{i} = dipout.filter{i}*noisecov*dipout.filter{i}';
    end
    dipout.noisecov(dip.outside) = {nan};
  end
end
