function [dipout] = music(dip, grad, vol, dat, varargin);

% MUSIC source localization using MUltiple SIgnal Classification
%
% This is a signal subspace method, which covers the techniques for
% multiple source localization by using the eigen structure of the
% measured data matrix.
%
% Use as
%   [dipout] = music(dip, grad, vol, dat, ...)
%
% Optional input arguments should be specified as key-value pairs and can be
%   'cov'              = data covariance matrix
%   'numcomponent'     = integer number
%   'feedback'         = 'none', 'gui', 'dial', 'textbar', 'text', 'textcr', 'textnl'
%   'reducerank'       = reduce the leadfield rank, can be 'no' or a number (e.g. 2)
%   'normalize'        = normalize the leadfield
%   'normalizeparam'   = parameter for depth normalization (default = 0.5)
%
% The original reference is
%   J.C. Mosher, P.S. Lewis and R.M. Leahy, "Multiple dipole modeling and
%   localization from spatiotemporal MEG data", IEEE Trans. Biomed.
%   Eng., pp 541-557, June, 1992.

% Copyright (C) 2004-2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% get the optional settings, or use the default value
cov            = keyval('cov',            varargin);
numcomponent   = keyval('numcomponent',   varargin); % this is required, see below
feedback       = keyval('feedback',       varargin); if isempty(feedback), feedback = 'text'; end
% these settings pertain to the forward model, the defaults are set in compute_leadfield
reducerank     = keyval('reducerank',     varargin);
normalize      = keyval('normalize',      varargin);
normalizeparam = keyval('normalizeparam', varargin);

if isempty(numcomponent)
  error('you must specify the number of signal components');
end

% ensure that these are row-vectors
dip.inside = dip.inside(:)';
dip.outside = dip.outside(:)';

Nchan = length(grad.label);
Ndip  = length(dip.inside);

if ~isempty(cov)
  % compute signal and noise subspace from covariance matrix
  [u, s, v] = svd(cov);
else
  % compute signal and noise subspace from average data matrix
  [u, s, v] = svd(dat);
end
% select the noise subspace, c.f. equation 25
us = u(:,(numcomponent+1):end);
ps = us * us';

% allocate space to hold the result
jr = zeros(length(dip.inside)+length(dip.outside),1);

ft_progress('init', feedback, 'computing music metric');
for i=1:length(dip.inside)

  ft_progress(i/length(dip.inside), 'computing music metric %d/%d\n', i, length(dip.inside));
  i = dip.inside(i);

  if isfield(dip, 'leadfield')
    % reuse the leadfield that was previously computed
    lf = dip.leadfield{i};
  elseif isfield(dip, 'mom')
    % compute the leadfield for a fixed dipole orientation
    lf = ft_compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam) * dip.mom(:,i);
  else
    % compute the leadfield
    lf = ft_compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam);
  end

  % compute the MUSIC metric, c.f. equation 26
  jr(i) = (norm(ps * lf)./norm(lf)).^2;
  % as described in the Mosher 1992 paper on page 550, "...the general approach is to
  % evaluare Jr(i) over a fine three-dimensional grid, plot its inverse,
  % and look for p sharp spikes..."

end
ft_progress('close');

% locations outside the head get assigned a nan
jr(dip.outside) = nan;

% assign the output data
dipout.jr = jr(:);  % ensure that it is a column vector

% add other descriptive information to the output source model
dipout.pos     = dip.pos;
dipout.inside  = dip.inside;
dipout.outside = dip.outside;

