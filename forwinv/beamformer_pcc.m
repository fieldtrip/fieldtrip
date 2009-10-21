function [dipout] = beamformer_pcc(dip, grad, vol, dat, Cf, varargin)

% BEAMFORMER_PCC implements an experimental beamformer based on partial canonical 
% correlations or coherences.
%
% Use and/or distribution of this function outside the F.C. Donders
% Centre is strictly prohibited!
%
% Copyright (C) 2005-2008, Robert Oostenveld & Jan-Mathijs Schoffelen

% $Log: beamformer_pcc.m,v $
% Revision 1.18  2009/01/07 13:12:03  jansch
% *** empty log message ***
%
% Revision 1.17  2009/01/06 10:25:50  roboos
% changed & into &&
%
% Revision 1.16  2008/10/10 13:01:18  release
% fixed bug: switched order of lambda processing (i.e. as "10%") and filter computation (thanks to Joachim).
% removed double warning for rank deficient.
%
% Revision 1.15  2008/08/13 13:47:42  roboos
% updated documentation
%
% Revision 1.14  2008/07/02 16:03:05  roboos
% fixed bug in % lambda
%
% Revision 1.13  2008/07/02 07:57:33  roboos
% allow specification of percentage noise in lambda, relative to trace(cov)/nchans
%
% Revision 1.12  2008/03/18 13:01:17  roboos
% added optional argument normalizeparam, is passed onto compute_leadfield
%
% Revision 1.11  2007/12/11 11:17:49  roboos
% fixed bug in handling of prespecified dipole moment
%
% Revision 1.10  2006/10/12 10:15:36  roboos
% close the progress indicator
%
% Revision 1.9  2006/10/12 10:11:24  roboos
% fixed bug in selecting dipoles on the inside positions only
% output cell-arrays are [] for outside points
% removed catch for mom, fixed dipoles should also work
%
% Revision 1.8  2006/10/12 09:05:11  roboos
% restructured code to behave as a stand-alone function
%
% Revision 1.6  2006/05/04 08:03:37  roboos
% fixed a bug in the concatenation of the filters, which applied when the channel order was EMG-MEG instead of the more usual MEG-EMG
%
% Revision 1.5  2006/03/01 10:47:06  roboos
% test rank of meg part only instead of full CSD matrix
%
% Revision 1.4  2006/03/01 08:07:33  roboos
% added support for realfilter option
%
% Revision 1.3  2005/11/21 09:33:49  roboos
% added handling of optional reducerank/normalize/powmethod input arguments
%
% Revision 1.2  2005/11/08 11:04:36  roboos
% implemented support for normalize and reducerank using compute_leadfield
% changed from optarg structure to keyval function
%
% Revision 1.1  2005/06/08 16:36:27  roboos
% new implementation, still in experimental stage
% DO NOT DISTRIBUTE
%

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
  [dip.inside, dip.outside] = find_inside_vol(dip.pos, vol);
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
  rf = compute_leadfield(refdip, grad, vol, 'reducerank', reducerank, 'normalize', normalize);
else
  rf = [];
end

if ~isempty(supdip)
  sf = compute_leadfield(supdip, grad, vol, 'reducerank', reducerank, 'normalize', normalize);
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
if isrankdeficient && ~isfield(dip, 'filter')
  warning('cross-spectral density matrix is rank deficient')
end

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
progress('init', feedback, 'beaming sources\n');

for i=1:size(dip.pos,1)
  if isfield(dip, 'leadfield') && ~isfield(dip, 'mom'),
    % reuse the leadfield that was previously computed
    lf = dip.leadfield{i};
  elseif isfield(dip, 'leadfield') && isfield(dip, 'mom'),
    lf = dip.leadfield{i} * dip.mom(:,i);
  elseif ~isfield(dip, 'leadfield') && isfield(dip, 'mom')
    % compute the leadfield for a fixed dipole orientation
    lf = compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam) * dip.mom(:,i);
  else
    % compute the leadfield
    lf = compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam);
  end

  % concatenate scandip, refdip and supdip
  lfa = [lf rf sf];

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

  progress(i/size(dip.pos,1), 'beaming source %d from %d\n', i, size(dip.pos,1));
end % for all dipoles

progress('close');

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
%   $Revision: 1.18 $  $Date: 2009/01/07 13:12:03 $
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

