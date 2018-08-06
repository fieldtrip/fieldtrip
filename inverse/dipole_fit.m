function [dipout] = dipole_fit(dip, sens, headmodel, dat, varargin)

% DIPOLE_FIT performs an equivalent current dipole fit with a single
% or a small number of dipoles to explain an EEG or MEG scalp topography.
%
% Use as
%   [dipout] = dipole_fit(dip, sens, headmodel, dat, ...)
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
% The constraints on the source model are specified in a structure
%   constr.symmetry   = boolean, dipole positions are symmetrically coupled to each other
%   constr.fixedori   = boolean, keep dipole orientation fixed over whole data window
%   constr.rigidbody  = boolean, keep relative position of multiple dipoles fixed
%   constr.mirror     = vector, used for symmetric dipole models
%   constr.reduce     = vector, used for symmetric dipole models
%   constr.expand     = vector, used for symmetric dipole models
%   constr.sequential = boolean, fit different dipoles to sequential slices of the data
%
% The maximum likelihood estimation implements
%   Lutkenhoner B. "Dipole source localization by means of maximum
%   likelihood estimation I. Theory and simulations" Electroencephalogr Clin
%   Neurophysiol. 1998 Apr;106(4):314-21.

% Copyright (C) 2003-2016, Robert Oostenveld
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

% It is neccessary to provide backward compatibility support for the old function call
% in case people want to use it in conjunction with EEGLAB and the dipfit1 plugin.
% old style: function [dipout] = dipole_fit(dip, dat, sens, headmodel, constr), where constr is optional
% new style: function [dipout] = dipole_fit(dip, sens, headmodel, dat, varargin), where varargin is in key-value pairs
if nargin==4 && ~isstruct(sens) && isstruct(dat)
  % looks like old style, the order of the input arguments has to be changed
  ft_warning('converting from old style input\n');
  olddat    = sens;
  oldsens   = headmodel;
  oldhdm    = dat;
  dat       = olddat;
  sens      = oldsens;
  headmodel = oldhdm;
elseif nargin==5  && ~isstruct(sens) && isstruct(dat)
  % looks like old style, the order of the input arguments has to be changed
  % furthermore the additional constraint has to be fixed
  ft_warning('converting from old style input\n');
  olddat    = sens;
  oldsens   = headmodel;
  oldhdm    = dat;
  dat       = olddat;
  sens      = oldsens;
  headmodel = oldhdm;
  varargin = {'constr', varargin{1}};  % convert into a  key-value pair
else
  % looks like new style, i.e. with optional key-value arguments
  % this is dealt with below
end

constr         = ft_getopt(varargin, 'constr'         ); % default is not to have constraints
metric         = ft_getopt(varargin, 'metric',    'rv');
checkinside    = ft_getopt(varargin, 'checkinside', false);
display        = ft_getopt(varargin, 'display', 'iter');
optimfun       = ft_getopt(varargin, 'optimfun'       ); if isa(optimfun, 'char'),  optimfun = str2func(optimfun); end
maxiter        = ft_getopt(varargin, 'maxiter'        );
reducerank     = ft_getopt(varargin, 'reducerank'     ); % for leadfield computation
normalize      = ft_getopt(varargin, 'normalize'      ); % for leadfield computation
normalizeparam = ft_getopt(varargin, 'normalizeparam' ); % for leadfield computation
weight         = ft_getopt(varargin, 'weight'         ); % for maximum likelihood estimation

if isfield(constr, 'mirror')
  % for backward compatibility
  constr.symmetry = true;
end
constr.symmetry   = ft_getopt(constr, 'symmetry', false);
constr.fixedori   = ft_getopt(constr, 'fixedori', false);
constr.rigidbody  = ft_getopt(constr, 'rigidbody', false);
constr.sequential = ft_getopt(constr, 'sequential', false);

if isempty(optimfun)
  % determine whether the MATLAB Optimization toolbox is available and can be used
  if ft_hastoolbox('optim')
    optimfun = @fminunc;
  else
    optimfun = @fminsearch;
  end
end

if isempty(maxiter)
  % set a default for the maximum number of iterations, depends on the optimization function
  if isequal(optimfun, @fminunc)
    maxiter = 1000;
  else
    maxiter = 3000;
  end
end

% determine whether it is EEG or MEG
iseeg = ft_senstype(sens, 'eeg');
ismeg = ft_senstype(sens, 'meg');

if ismeg && iseeg
  % this is something that I might implement in the future
  ft_error('simultaneous EEG and MEG not supported');
elseif iseeg
  % ensure that the potential data is average referenced, just like the model potential
  dat = avgref(dat);
end

% ensure correct dipole position and moment specification
dip = fixdipole(dip);

% convert the dipole model parameters into the non-linear parameter vector that will be optimized
[param, constr] = dipolemodel2param(dip.pos, dip.mom, constr);

% determine the scale
scale = ft_scalingfactor(sens.unit, 'cm');

% set the parameters for the optimization function
if isequal(optimfun, @fminunc)
  options = optimset(...
    'TolFun',1e-9,...
    'TypicalX',scale*ones(size(param)),...
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
  ft_warning('unknown optimization function "%s", using default parameters', func2str(optimfun));
end

% perform the optimization with either the fminsearch or fminunc function
[param, fval, exitflag, output] = optimfun(@dipfit_error, param, options, dat, sens, headmodel, constr, metric, checkinside, reducerank, normalize, normalizeparam, weight);

if exitflag==0
  ft_error('Maximum number of iterations exceeded before reaching the minimum, please try with another initial guess.')
end

% do the linear optimization of the dipole moment parameters
% the error is not interesting any more, only the dipole moment is relevant
[err, mom] = dipfit_error(param, dat, sens, headmodel, constr, metric, checkinside, reducerank, normalize, normalizeparam, weight);

% convert the non-linear parameter vector into the dipole model parameters
[pos, ori] = param2dipolemodel(param, constr);

% return the optimal dipole parameters
dipout.pos = pos;
% return the optimal dipole moment and (optionally) the orientation
if ~isempty(ori)
  dipout.mom  = ori;  % dipole orientation as vector
  dipout.ampl = mom;  % dipole strength
else
  dipout.mom  = mom;  % dipole moment as vector or matrix, which represents both the orientation and strength as vector
end

% ensure correct dipole position and moment specification
dipout = fixdipole(dipout);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIPOLEMODEL2PARAM takes the initial guess for the diople model and converts it
% to a set of parameters that needs to be optimized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [param, constr] = dipolemodel2param(pos, ori, constr)

% reformat the position parameters in case of multiple dipoles, this
% should result in the matrix changing from [x1 y1 z1; x2 y2 z2] to
% [x1 y1 z1 x2 y2 z2] for the constraints to work
param = reshape(pos', 1, numel(pos));

% add the orientation to the nonlinear parameters
if constr.fixedori
  numdip = size(pos,1);
  for i=1:numdip
    % add the orientation to the list of parameters
    [th, phi, r] = cart2sph(ori(1,i), ori(2,i), ori(3,i));
    param = [param th phi];
  end
end

if constr.symmetry && constr.rigidbody
  ft_error('simultaneous symmetry and rigidbody constraints are not supported')
  
elseif constr.symmetry
  % reduce the number of parameters to be fitted according to the constraints
  % select a subset, the other sources will be re-added by the const.mirror field
  param = param(constr.reduce);
  
elseif constr.rigidbody
  constr.coilpos = param;    % store the head localizer coil positions
  param = [0 0 0 0 0 0];     % start with an initial translation and rotation of zero
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAM2DIPOLEMODEL takes the parameters and constraints and converts them into a
% diople model for which the leadfield and residual error can be computed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos, ori] = param2dipolemodel(param, constr)

if constr.symmetry && constr.rigidbody
  ft_error('simultaneous symmetry and rigidbody constraints are not supported')
  
elseif constr.symmetry
  param  = constr.mirror .* param(constr.expand);
  
elseif constr.rigidbody
  numdip    = numel(constr.coilpos)/3;
  pos       = reshape(constr.coilpos, 3, numdip); % convert from vector into 3xN matrix
  pos(4,:)  = 1;
  transform = rigidbody(param);                   % this is a 4x4 homogenous transformation matrix
  pos       = transform * pos;                    % apply the homogenous transformation matrix
  param     = reshape(pos(1:3,:), 1, 3*numdip);
  clear pos                                       % the actual pos will be constructed from param further down 
end

if constr.fixedori
  numdip = numel(param)/5;
  ori = zeros(3,numdip);
  for i=1:numdip
    th  = param(end-(2*i)+1);
    phi = param(end-(2*i)+2);
    [ori(1,i), ori(2,i), ori(3,i)] = sph2cart(th, phi, 1);
  end
  pos = reshape(param(1:(numdip*3)), 3, numdip)'; % convert into a Ndip*3 matrix
else
  numdip = numel(param)/3;
  pos = reshape(param, 3, numdip)'; % convert into a Ndip*3 matrix
  ori = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIPFIT_ERROR computes the error between measured and model data
% and can be used for non-linear fitting of dipole position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, mom] = dipfit_error(param, dat, sens, headmodel, constr, metric, checkinside, reducerank, normalize, normalizeparam, weight)

% flush pending graphics events, ensure that fitting is interruptible
drawnow;
if ~isempty(get(0, 'currentfigure')) && strcmp(get(gcf, 'tag'), 'stop')
  % interrupt the fitting
  close;
  ft_error('USER ABORT');
end;

% convert the non-linear parameter vector into the dipole model parameters
[pos, ori] = param2dipolemodel(param, constr);

% check whether the dipole is inside the source compartment
if checkinside
  inside = ft_inside_vol(pos, headmodel);
  if ~all(inside)
    ft_error('Dipole is outside the source compartment');
  end
end

% construct the leadfield matrix for all dipoles
lf = ft_compute_leadfield(pos, sens, headmodel, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam);
if ~isempty(ori)
  lf = lf * ori;
end

% compute the optimal dipole moment and the model error
if ~isempty(weight)
  % maximum likelihood estimation using the weigth matrix
  if constr.sequential
    ft_error('not supported');
  else
    mom = pinv(lf'*weight*lf)*lf'*weight*dat;  % Lutkenhoner equation 5
    dif = dat - lf*mom;
  end
  % compute the generalized goodness-of-fit measure
  switch metric
    case 'rv' % relative residual variance
      num   = dif' * weight * dif;
      denom = dat' * weight * dat;
      err   = sum(num(:)) ./ sum(denom(:)); % Lutkenhonner equation 7, except for the gof=1-rv
    case 'var' % residual variance
      num   = dif' * weight * dif;
      err   = sum(num(:));
    otherwise
      ft_error('Unsupported error metric for maximum likelihood dipole fitting');
  end
else
  % ordinary least squares, this is the same as MLE with weight=eye(nchans,nchans)
  if constr.sequential
    % the number of slices is the same as the number of dipoles
    % each slice has a number of frames (time points) in it
    % so the data can be nchan*ndip or nchan*(ndip*nframe)
    numdip   = numel(pos)/3;
    numframe = size(dat,2)/numdip;
    
    % do a sainty check on the number of frames, see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3119
    assert(numframe>0 && numframe==round(numframe), 'the number of frames should be a positive integer');
    
    mom = zeros(3*numdip, numdip*numframe);
    for i=1:numdip
      dipsel   = (1:3)        + 3*(i-1);         % 1:3 for the first dipole, 4:6 for the second dipole, ...
      framesel = (1:numframe) + numframe*(i-1);  % 1:numframe for the first, (numframe+1):(2*numframe) for the second, ...
      mom(dipsel,framesel) = pinv(lf(:,dipsel))*dat(:,framesel);
    end
  else
    mom = pinv(lf)*dat;
  end
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
      ft_error('Unsupported error metric for dipole fitting');
  end
end

if ~isreal(err)
  % this happens for complex valued data, i.e. when fitting a dipole to spectrally decomposed data
  % the error function should return a positive valued real number, otherwise fminunc fails
  err = abs(err);
end
