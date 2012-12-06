function [dipout] = beamformer_pcc(dip, grad, vol, dat, Cf, varargin)

% BEAMFORMER_PCC implements an experimental beamformer based on partial canonical 
% correlations or coherences.

% Copyright (C) 2005-2008, Robert Oostenveld & Jan-Mathijs Schoffelen
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
refchan        = keyval('refchan',       varargin);
refdip         = keyval('refdip',        varargin);
supchan        = keyval('supchan',       varargin);
supdip         = keyval('supdip',        varargin);
% these settings pertain to the forward model, the defaults are set in compute_leadfield
reducerank     = keyval('reducerank',     varargin);
normalize      = keyval('normalize',      varargin);
normalizeparam = keyval('normalizeparam', varargin);
% these optional settings have defaults
feedback       = keyval('feedback',      varargin); if isempty(feedback),      feedback = 'text';            end
keepcsd        = keyval('keepcsd',       varargin); if isempty(keepcsd),       keepcsd = 'no';               end
keepfilter     = keyval('keepfilter',    varargin); if isempty(keepfilter),    keepfilter = 'no';            end
keepleadfield  = keyval('keepleadfield', varargin); if isempty(keepleadfield), keepleadfield = 'no';         end
keepmom        = keyval('keepmom',       varargin); if isempty(keepmom),       keepmom = 'yes';              end
lambda         = keyval('lambda',        varargin); if isempty(lambda  ),      lambda = 0;                   end
projectnoise   = keyval('projectnoise',  varargin); if isempty(projectnoise),  projectnoise = 'yes';         end
realfilter     = keyval('realfilter',    varargin); if isempty(realfilter),    realfilter = 'yes';           end
fixedori       = ft_getopt(varargin,'fixedori','no');
fixedori       = strcmp(fixedori,      'yes');

% convert the yes/no arguments to the corresponding logical values
keepcsd        = strcmp(keepcsd,       'yes');  % see below
keepfilter     = strcmp(keepfilter,    'yes');
keepleadfield  = strcmp(keepleadfield, 'yes');
keepmom        = strcmp(keepmom,       'yes');
projectnoise   = strcmp(projectnoise,  'yes');
realfilter     = strcmp(realfilter,    'yes');

% the postprocessing of the pcc beamformer always requires the csd matrix
keepcsd = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the dipole positions that are inside/outside the brain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(dip, 'inside') & ~isfield(dip, 'outside');
  insideLogical = ft_inside_vol(dip.pos, vol);
  dip.inside = find(insideLogical);
  dip.outside = find(~dip.inside);
elseif isfield(dip, 'inside') & ~isfield(dip, 'outside');
  dip.outside    = setdiff(1:size(dip.pos,1), dip.inside);
elseif ~isfield(dip, 'inside') & isfield(dip, 'outside');
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
dip.pos     = dip.pos(dip.inside, :);
dip.inside  = 1:size(dip.pos,1);
dip.outside = [];

if ~isempty(refdip)
  rf = ft_compute_leadfield(refdip, grad, vol, 'reducerank', reducerank, 'normalize', normalize);
else
  rf = [];
end

if ~isempty(supdip)
  sf = ft_compute_leadfield(supdip, grad, vol, 'reducerank', reducerank, 'normalize', normalize);
else
  sf = [];
end

refchan  = refchan;               % these can be passed as optional inputs
supchan  = supchan;               % these can be passed as optional inputs
megchan  = setdiff(1:size(Cf,1), [refchan supchan]);
Nrefchan = length(refchan);
Nsupchan = length(supchan);
Nmegchan = length(megchan);
Nchan    = size(Cf,1);            % should equal Nmegchan + Nrefchan + Nsupchan
Cmeg     = Cf(megchan,megchan);   %  the filter uses the csd between all MEG channels

isrankdeficient = (rank(Cmeg)<size(Cmeg,1));

% it is difficult to give a quantitative estimate of lambda, therefore also
% support relative (percentage) measure that can be specified as string (e.g. '10%')
if ~isempty(lambda) && ischar(lambda) && lambda(end)=='%'
    ratio = sscanf(lambda, '%f%%');
    ratio = ratio/100;
    lambda = ratio * trace(Cmeg)/size(Cmeg,1);
end

if projectnoise
  % estimate the noise power, which is further assumed to be equal and uncorrelated over channels
  if isrankdeficient
    % estimated noise floor is equal to or higher than lambda
    noise = lambda;
  else
    % estimate the noise level in the covariance matrix by the smallest singular value
    noise = svd(Cmeg);
    noise = noise(end);
    % estimated noise floor is equal to or higher than lambda
    noise = max(noise, lambda);
  end
end

if realfilter
  % construct the filter only on the real part of the CSD matrix, i.e. filter is real
  invCmeg = pinv(real(Cmeg) + lambda*eye(Nmegchan));
else
  % construct the filter on the complex CSD matrix, i.e. filter contains imaginary component as well
  % this results in a phase rotation of the channel data if the filter is applied to the data
  invCmeg = pinv(Cmeg + lambda*eye(Nmegchan));
end

% start the scanning with the proper metric
ft_progress('init', feedback, 'beaming sources\n');

for i=1:size(dip.pos,1)
  if isfield(dip, 'leadfield') && isfield(dip, 'mom') && size(dip.mom, 1)==size(dip.leadfield{i}, 2)
    % reuse the leadfield that was previously computed and project
    lf = dip.leadfield{i} * dip.mom(:,i);
  elseif  isfield(dip, 'leadfield') &&  isfield(dip, 'mom')
    % reuse the leadfield that was previously computed but don't project
    lf = dip.leadfield{i};
  elseif isfield(dip, 'leadfield') && ~isfield(dip, 'mom'),
    % reuse the leadfield that was previously computed
    lf = dip.leadfield{i};
  elseif ~isfield(dip, 'leadfield') && isfield(dip, 'mom')
    % compute the leadfield for a fixed dipole orientation
    lf = ft_compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam) * dip.mom(:,i);
  else
    % compute the leadfield
    lf = ft_compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam);
  end

  % concatenate scandip, refdip and supdip
  lfa = [lf rf sf];

  if fixedori
    if isempty(refdip) && isempty(supdip) && isempty(refchan) && isempty(supchan) && (size(lf,2)==3)
      % compute the leadfield for the optimal dipole orientation
      % subsequently the leadfield for only that dipole orientation will
      % be used for the final filter computation
      if isfield(dip, 'filter') && size(dip.filter{i},1)~=1
        filt = dip.filter{i};
      else
        filt = pinv(lfa' * invCmeg * lfa) * lfa' * invCmeg;
      end
      [u, s, v] = svd(real(filt * Cmeg * ctranspose(filt)));
      maxpowori = u(:,1);
      eta = s(1,1)./s(2,2);
      lfa  = lfa * maxpowori;
      dipout.ori{i} = maxpowori;
      dipout.eta{i} = eta;
    else
      warning_once('Ignoring ''fixedori''. The fixedori option is supported only if there is ONE dipole for location.')
    end
  end
  
  if isfield(dip, 'filter')
    % use the provided filter
    filt = dip.filter{i};
  else
    % construct the spatial filter
    filt = pinv(lfa' * invCmeg * lfa) * lfa' * invCmeg;                % use PINV/SVD to cover rank deficient leadfield
  end

  % concatenate the source filters with the channel filters
  Ndip = size(lfa, 2);
  filtn = zeros(Ndip+Nrefchan+Nsupchan, Nmegchan+Nrefchan+Nsupchan);
  % this part of the filter relates to the sources
  filtn(1:Ndip,megchan) = filt;
  % this part of the filter relates to the channels 
  filtn((Ndip+1):end,setdiff(1:(Nmegchan+Nrefchan+Nsupchan), megchan)) = eye(Nrefchan+Nsupchan);
  filt = filtn;
  clear filtn

  if keepcsd
    dipout.csd{i} = filt * Cf * ctranspose(filt);
  end
  if projectnoise
    dipout.noisecsd{i} = noise * (filt * ctranspose(filt));
  end
  if keepmom && ~isempty(dat)
    dipout.mom{i} = filt * dat;
  end
  if keepfilter
    dipout.filter{i} = filt;
  end
  if keepleadfield
    dipout.leadfield{i} = lf;
  end

  ft_progress(i/size(dip.pos,1), 'beaming source %d from %d\n', i, size(dip.pos,1));
end % for all dipoles

ft_progress('close');

dipout.inside  = dip.originside;
dipout.outside = dip.origoutside;
dipout.pos     = dip.origpos;

% remember how all components in the output csd should be interpreted
scandiplabel = repmat({'scandip'}, 1, size(lf, 2));    % based on last leadfield
refdiplabel  = repmat({'refdip'},  1, size(rf, 2));
supdiplabel  = repmat({'supdip'},  1, size(sf, 2));
refchanlabel = repmat({'refchan'}, 1, Nrefchan);
supchanlabel = repmat({'supchan'}, 1, Nsupchan);
% concatenate all the labels
dipout.csdlabel = [scandiplabel refdiplabel supdiplabel refchanlabel supchanlabel];

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
if isfield(dipout, 'csd')
  dipout.csd(dipout.inside)  = dipout.csd;
  dipout.csd(dipout.outside) = {[]};
end
if isfield(dipout, 'noisecsd')
  dipout.noisecsd(dipout.inside)  = dipout.noisecsd;
  dipout.noisecsd(dipout.outside) = {[]};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to compute the pseudo inverse. This is the same as the 
% standard Matlab function, except that the default tolerance is twice as
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

