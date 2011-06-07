function [dipout] = dipole_fit(dip, sens, vol, dat, varargin)

% DIPOLE_FIT performs an equivalent current dipole fit with a single
% or a small number of dipoles to explain an EEG or MEG scalp topography.
%
% Use as
%   [dipout] = dipole_fit(dip, sens, vol, dat, ...)
%
% Additional input arguments should be specified as key-value pairs and can include
%   'constr'      = Structure with constraints
%   'display'     = Level of display [ off | iter | notify | final ]
%   'optimfun'    = Function to use [fminsearch | fminunc ]
%   'maxiter'     = Maximum number of function evaluations allowed [ positive integer ]
%   'metric'      = Error measure to be minimised [ rv | var | abs ]
%   'checkinside' = Boolean flag to check whether dipole is inside source compartment [ 0 | 1 ]
%   'weight'      = weight matrix for maximum likelihood estimation, e.g. inverse noise covariance
%
% The following optional input arguments relate to the computation of the leadfields
%   'reducerank'      = 'no' or number
%   'normalize'       = 'no', 'yes' or 'column'
%   'normalizeparam'  = parameter for depth normalization (default = 0.5)
%
% The maximum likelihood estimation implements
%   Lutkenhoner B. "Dipole source localization by means of maximum
%   likelihood estimation I. Theory and simulations" Electroencephalogr Clin
%   Neurophysiol. 1998 Apr;106(4):314-21.

% Copyright (C) 2003-2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% It is neccessary to provide backward compatibility support for the old function call
% in case people want to use it in conjunction with EEGLAB and the dipfit1 plugin.
% old style: function [dipout] = dipole_fit(dip, dat, sens, vol, constr), where constr is optional
% new style: function [dipout] = dipole_fit(dip, sens, vol, dat, varargin), where varargin is in key-value pairs
if nargin==4 && ~isstruct(sens) && isstruct(dat)
  % looks like old style, the order of the input arguments has to be changed
  warning('converting from old style input\n');
  olddat   = sens;
  oldsens  = vol;
  oldvol   = dat;
  dat      = olddat;
  sens     = oldsens;
  vol      = oldvol;
elseif nargin==5  && ~isstruct(sens) && isstruct(dat)
  % looks like old style, the order of the input arguments has to be changed
  % furthermore the additional constraint has to be fixed
  warning('converting from old style input\n');
  olddat   = sens;
  oldsens  = vol;
  oldvol   = dat;
  dat      = olddat;
  sens     = oldsens;
  vol      = oldvol;
  varargin = {'constr', varargin{1}};  % convert into a  key-value pair
else
  % looks like new style, i.e. with optional key-value arguments
  % this is dealt with below
end

constr         = keyval('constr',         varargin); % default is not to have constraints
metric         = keyval('metric',         varargin); if isempty(metric),        metric = 'rv';                end
checkinside    = keyval('checkinside',    varargin); if isempty(checkinside),   checkinside = 0;              end
display        = keyval('display',        varargin); if isempty(display),       display = 'iter';             end
optimfun       = keyval('optimfun',       varargin); if isa(optimfun, 'char'),  optimfun = str2fun(optimfun); end
maxiter        = keyval('maxiter',        varargin);
reducerank     = keyval('reducerank',     varargin); % for leadfield computation
normalize      = keyval('normalize' ,     varargin); % for leadfield computation
normalizeparam = keyval('normalizeparam', varargin); % for leadfield computation
weight         = keyval('weight',         varargin); % for maximum likelihood estimation

if isempty(optimfun)
  % determine whether the Matlab Optimization toolbox is available and can be used
  if ft_hastoolbox('optim')
    optimfun = @fminunc;
  else
    optimfun = @fminsearch;
  end
end

if isempty(maxiter)
  % set a default for the maximum number of iterations, depends on the optimization function
  if isequal(optimfun, @fminunc)
    maxiter = 100;
  else
    maxiter = 500;
  end
end

% determine whether it is EEG or MEG
iseeg = ft_senstype(sens, 'eeg');
ismeg = ft_senstype(sens, 'meg');

if ismeg && iseeg
  % this is something that I might implement in the future
  error('simultaneous EEG and MEG not supported');
elseif iseeg
  % ensure that the potential data is average referenced, just like the model potential
  dat = avgref(dat);
end

% ensure correct dipole position and moment specification
dip = fixdipole(dip);

% reformat the position parameters in case of multiple dipoles, this
% should result in the matrix changing from [x1 y1 z1; x2 y2 z2] to
% [x1 y1 z1 x2 y2 z2] for the constraints to work
numdip = size(dip.pos, 1);
param  = dip.pos';
param  = param(:)';

% add the orientation to the nonlinear parameters
if isfield(constr, 'fixedori') && constr.fixedori
  for i=1:numdip
    % add the orientation to the list of parameters
    [th, phi, r] = cart2sph(dip.mom(1,i), dip.mom(2,i), dip.mom(3,i));
    param = [param th phi];
  end
end

% reduce the number of parameters to be fitted according to the constraints
if isfield(constr, 'mirror')
  param = param(constr.reduce);
end

% set the parameters for the optimization function
if isequal(optimfun, @fminunc)
  options = optimset(...
    'TolFun',1e-9,...
    'TypicalX',ones(size(param)),...
    'LargeScale','off',...
    'HessUpdate','bfgs',...
    'MaxIter',maxiter,...
    'MaxFunEvals',2*maxiter*length(param),...
    'Display',display);
elseif isequal(optimfun, @fminsearch)
  options = optimset(...
    'MaxIter',maxiter,...
    'MaxFunEvals',2*maxiter*length(param),...
    'Display',display);
else
  warning('unknown optimization function "%s", using default parameters', func2str(optimfun));
end

% perform the optimization with either the fminsearch or fminunc function
[param, fval, exitflag, output] = optimfun(@dipfit_error, param, options, dat, sens, vol, constr, metric, checkinside, reducerank, normalize, normalizeparam, weight);

if exitflag==0
  error('Maximum number of iterations exceeded before reaching the minimum, please try with another initial guess.')
end

% do linear optimization of dipole moment parameters
[err, mom] = dipfit_error(param, dat, sens, vol, constr, metric, checkinside, reducerank, normalize, normalizeparam, weight);

% expand the number of parameters according to the constraints
if isfield(constr, 'mirror')
  param = constr.mirror .* param(constr.expand);
end

% get the dipole position and orientation
if isfield(constr, 'fixedori') && constr.fixedori
  numdip = numel(param)/5;
  ori = zeros(3,numdip);
  for i=1:numdip
    th  = param(end-(2*i)+1);
    phi = param(end-(2*i)+2);
    [ori(1,i), ori(2,i), ori(3,i)] = sph2cart(th, phi, 1);
  end
  pos = reshape(param(1:(numdip*3)), 3, numdip)';
else
  numdip = numel(param)/3;
  pos = reshape(param, 3, numdip)';
end

% return the optimal dipole parameters
dipout.pos = pos;
if isfield(constr, 'fixedori') && constr.fixedori
  dipout.mom  = ori;  % dipole orientation as vector
  dipout.ampl = mom;  % dipole strength
else
  dipout.mom  = mom;  % dipole moment as vector or matrix, which represents both the orientation and strength as vector
end

% ensure correct dipole position and moment specification
dipout = fixdipole(dipout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIPFIT_ERROR computes the error between measured and model data
% and can be used for non-linear fitting of dipole position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, mom] = dipfit_error(param, dat, sens, vol, constr, metric, checkinside, reducerank, normalize, normalizeparam, weight)

% flush pending graphics events, ensure that fitting is interruptible
drawnow;
if ~isempty(get(0, 'currentfigure')) && strcmp(get(gcf, 'tag'), 'stop')
  % interrupt the fitting
  close;
  error('USER ABORT');
end;

% expand the number of parameters according to the constraints
if isfield(constr, 'mirror')
  param = constr.mirror .* param(constr.expand);
end

% get the dipole positions and optionally also the orientation
if isfield(constr, 'fixedori') && constr.fixedori
  numdip = numel(param)/5;
  ori = zeros(3,numdip);
  for i=1:numdip
    th  = param(end-(2*i)+1);
    phi = param(end-(2*i)+2);
    [ori(1,i), ori(2,i), ori(3,i)] = sph2cart(th, phi, 1);
  end
  pos = reshape(param(1:(numdip*3)), 3, numdip)';
else
  numdip = numel(param)/3;
  pos = reshape(param, 3, numdip)';
end

% check whether the dipole is inside the source compartment
if checkinside
  inside = ft_inside_vol(pos, vol);
  if ~all(inside)
    error('Dipole is outside the source compartment');
  end
end

% construct the leadfield matrix for all dipoles
lf = ft_compute_leadfield(pos, sens, vol, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam);
if isfield(constr, 'fixedori') && constr.fixedori
  lf = lf * ori;
end

% compute the optimal dipole moment and the model error
if ~isempty(weight)
  % maximum likelihood estimation using the weigth matrix
  mom = pinv(lf'*weight*lf)*lf'*weight*dat;  % Lutkenhoner equation 5
  dif = dat - lf*mom;
  % compute the generalized goodness-of-fit measure
  switch metric
    case 'rv' % relative residual variance
      num   = dif' * weight * dif;
      denom = dat' * weight * dat;
      err   = sum(num(:)) ./ sum(denom(:)); % Lutkenhonner equation 7, except for the gof=1-rv
    case 'var' % residual variance
      num   = dif' * weight * dif';
      err   = sum(num(:));
    otherwise
      error('Unsupported error metric for maximum likelihood dipole fitting');
  end
else
  % ordinary least squares, this is the same as MLE with weight=eye(nchans,nchans)
  mom = pinv(lf)*dat;
  dif = dat - lf*mom;
  % compute the ordinary goodness-of-fit measures
  switch metric
    case 'rv' % relative residual variance
      err = sum(dif(:).^2) / sum(dat(:).^2);
    case 'var' % residual variance
      err = sum(dif(:).^2);
    case 'abs' % absolute difference
      err = sum(abs(dif));
    otherwise
      error('Unsupported error metric for dipole fitting');
  end
end

if ~isreal(err)
  % this happens for complex valued data, i.e. when fitting a dipole to spectrally decomposed data
  % the error function should return a positive valued real number, otherwise fminunc fails
  err = abs(err);
end
