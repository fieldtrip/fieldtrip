function [estimate] = ft_inverse_beamformer_dics(leadfield, Cf, varargin)

% FT_INVERSE_BEAMFORMER_DICS estimates the source power or source
% coherence according to the Dynamic Imaging of Coherent Sources
% method.
%
% Use as
%   estimate = ft_inverse_beamformer_dics(leadfield, Cf, ...)
% where
%   leadfield  = leadfield of the source of interest or a cell-array with leadfields for multiple sources
%   Cf         = cross-spectral density matrix of the data
% and
%   estimate   = structure with the estimated source parameters
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
% This implements Joachim Gross et al. 2001

% Copyright (C) 2003-2010, Robert Oostenveld

% these optional settings do not have defaults
Pr             = keyval('Pr',            varargin);
Cr             = keyval('Cr',            varargin);
refdip         = keyval('refdip',        varargin);
powmethod      = keyval('powmethod',     varargin); % the default for this is set below
realfilter     = keyval('realfilter',    varargin); % the default for this is set below
% these optional settings have defaults
feedback       = keyval('feedback',      varargin); if isempty(feedback),      feedback = 'text';            end
keepcsd        = keyval('keepcsd',       varargin); if isempty(keepcsd),       keepcsd = 'no';               end
keepfilter     = keyval('keepfilter',    varargin); if isempty(keepfilter),    keepfilter = 'no';            end
keepleadfield  = keyval('keepleadfield', varargin); if isempty(keepleadfield), keepleadfield = 'no';         end
lambda         = keyval('lambda',        varargin); if isempty(lambda  ),      lambda = 0;                   end
projectnoise   = keyval('projectnoise',  varargin); if isempty(projectnoise),  projectnoise = 'yes';         end
fixedori       = keyval('fixedori',      varargin); if isempty(fixedori),      fixedori = 'no';              end

% convert the yes/no arguments to the corresponding logical values
% FIXME use istrue
keepcsd        = strcmp(keepcsd,       'yes');
keepfilter     = strcmp(keepfilter,    'yes');
keepleadfield  = strcmp(keepleadfield, 'yes');
projectnoise   = strcmp(projectnoise,  'yes');
fixedori       = strcmp(fixedori,      'yes');

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

% dics has the following sub-methods, which depend on the additional input arguments
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
  error('invalid combination of input arguments for dics');
end

if ~iscell(leadfield)
  % the leadfield specifies a single source
  leadfield = {leadfield};
end
ndipoles = length(leadfield);

if ~isempty(Cr)
  % ensure that the cross-spectral density with the reference signal is a column matrix
  Cr = Cr(:);
end

isrankdeficient = (rank(Cf)<size(Cf,1));
if isrankdeficient
  warning('cross-spectral density matrix is rank deficient')
end

% it is difficult to give a quantitative estimate of lambda, therefore also
% support relative (percentage) measure that can be specified as string (e.g. '10%')
if ~isempty(lambda) && ischar(lambda) && lambda(end)=='%'
  ratio = sscanf(lambda, '%f%%');
  ratio = ratio/100;
  lambda = ratio * trace(Cf)/size(Cf,1);
end

if projectnoise
  % estimate the noise power, which is further assumed to be equal and uncorrelated over channels
  if isrankdeficient
    % estimated noise floor is equal to or higher than lambda
    noise = lambda;
  else
    % estimate the noise level in the covariance matrix by the smallest singular value
    noise = svd(Cf);
    noise = noise(end);
    % estimated noise floor is equal to or higher than lambda
    noise = max(noise, lambda);
  end
end

% the inverse only has to be computed once for all dipoles
if strcmp(realfilter, 'yes')
  % the filter is computed using only the leadfield and the inverse covariance or CSD matrix
  % therefore using the real-valued part of the CSD matrix here ensures a real-valued filter
  invCf = pinv(real(Cf) + lambda * eye(size(Cf)));
else
  invCf = pinv(Cf + lambda * eye(size(Cf)));
end

% start the scanning with the proper metric
ft_progress('init', feedback, 'scanning grid');
switch submethod

  case 'dics_power'
    % only compute power of a dipole at the grid positions
    for i=1:ndipoles
      lf = leadfield{i};

      if fixedori
        % compute the leadfield for the optimal dipole orientation
        % subsequently the leadfield for only that dipole orientation will be used for the final filter computation
        filt = pinv(lf' * invCf * lf) * lf' * invCf;
        [u, s, v] = svd(real(filt * Cf * ctranspose(filt)));
        maxpowori = u(:,1);
        eta = s(1,1)./s(2,2);
        lf  = lf * maxpowori;
        estimate.ori{i} = maxpowori;
        estimate.eta{i} = eta;
      end

      % construct the spatial filter
      filt = pinv(lf' * invCf * lf) * lf' * invCf;                % Gross eqn. 3, use PINV/SVD to cover rank deficient leadfield
      csd = filt * Cf * ctranspose(filt);                         % Gross eqn. 4 and 5

      % assign the output values
      if powlambda1
        estimate.pow(i) = lambda1(csd);                             % compute the power at the dipole location, Gross eqn. 8
      elseif powtrace
        estimate.pow(i) = real(trace(csd));                         % compute the power at the dipole location
      end
      if keepcsd
        estimate.csd{i} = csd;
      end
      if projectnoise
        if powlambda1
          estimate.noise(i) = noise * lambda1(filt * ctranspose(filt));
        elseif powtrace
          estimate.noise(i) = noise * real(trace(filt * ctranspose(filt)));
        end
        if keepcsd
          estimate.noisecsd{i} = noise * filt * ctranspose(filt);
        end
      end
      if keepfilter
        estimate.filter{i} = filt;
      end
      if keepleadfield
        estimate.leadfield{i} = lf;
      end
      ft_progress(i/ndipoles, 'scanning grid %d/%d\n', i, ndipoles);
    end

  case 'dics_refchan'
    % compute cortico-muscular coherence, using reference cross spectral density
    for i=1:ndipoles

      % get the leadfield for this source
      lf = leadfield{i};

      if fixedori
        % compute the leadfield for the optimal dipole orientation
        % subsequently the leadfield for only that dipole orientation will be used for the final filter computation
        filt = pinv(lf' * invCf * lf) * lf' * invCf;
        [u, s, v] = svd(real(filt * Cf * ctranspose(filt)));
        maxpowori = u(:,1);
        lf  = lf * maxpowori;
        estimate.ori{i} = maxpowori;
      end

      % construct the spatial filter
      filt = pinv(lf' * invCf * lf) * lf' * invCf;                     % use PINV/SVD to cover rank deficient leadfield

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

      estimate.pow(i) = pow;
      estimate.coh(i) = coh;
      if keepcsd
        estimate.csd{i} = csd;
      end
      if projectnoise
        if powlambda1
          estimate.noise(i) = noise * lambda1(filt * ctranspose(filt));
        elseif powtrace
          estimate.noise(i) = noise * real(trace(filt * ctranspose(filt)));
        end
        if keepcsd
          estimate.noisecsd{i} = noise * filt * ctranspose(filt);
        end
      end
      if keepfilter
        estimate.filter{i} = filt;
      end
      ft_progress(i/ndipoles, 'scanning grid %d/%d\n', i, ndipoles);
    end

  case 'dics_refdip'

    if fixedori
      error('fixed orientations are not supported for beaming cortico-cortical coherence');
    end

    % get the leadfield of the reference source
    lf1 = refdip;

    % construct the spatial filter for the first (reference) dipole location
    filt1 = pinv(lf1' * invCf * lf1) * lf1' * invCf;       % use PINV/SVD to cover rank deficient leadfield
    if powlambda1
      Pref = lambda1(filt1 * Cf * ctranspose(filt1));      % compute the power at the first dipole location, Gross eqn. 8
    elseif powtrace
      Pref = real(trace(filt1 * Cf * ctranspose(filt1)));  % compute the power at the first dipole location
    end

    for i=1:ndipoles

      % get the leadfield for the second source, i.e. the one that is being scanned
      lf2 = leadfield{i};

      % construct the spatial filter for the second source
      filt2 = pinv(lf2' * invCf * lf2) * lf2' * invCf;     % use PINV/SVD to cover rank deficient leadfield
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

      estimate.pow(i) = pow;
      estimate.coh(i) = coh;
      if keepcsd
        estimate.csd{i} = csd;
      end
      if projectnoise
        if powlambda1
          estimate.noise(i) = noise * lambda1(filt2 * ctranspose(filt2));
        elseif powtrace
          estimate.noise(i) = noise * real(trace(filt2 * ctranspose(filt2)));
        end
        if keepcsd
          estimate.noisecsd{i} = noise * filt2 * ctranspose(filt2);
        end
      end
      ft_progress(i/ndipoles, 'scanning grid %d/%d\n', i, ndipoles);
    end

end % switch submethod

ft_progress('close');

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
% standard Matlab function, except that the default tolerance is twice as
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

