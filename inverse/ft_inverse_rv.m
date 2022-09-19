function [estimate] = ft_inverse_rv(sourcemodel, sens, headmodel, dat, varargin)

% FT_INVERSE_RV scan with a single dipole and computes the residual variance
% at each dipole location.
%
% Use as
%   [estimate] = ft_inverse_rv(sourcemodel, sens, headmodel, dat, ...)
% where
%   sourcemodel is the input source model, see FT_PREPARE_SOURCEMODEL
%   sens        is the gradiometer or electrode definition, see FT_DATATYPE_SENS
%   headmodel   is the volume conductor definition, see FT_PREPARE_HEADMODEL
%   dat         is the data matrix with the ERP or ERF
% and
%   estimate    contains the estimated source parameters
%
% Additional input arguments should be specified as key-value pairs and can include
%   'feedback'         = can be 'none', 'gui', 'dial', 'textbar', 'text', 'textcr', 'textnl' (default = 'text')
%
% These options influence the forward computation of the leadfield
%   'reducerank'      = 'no' or number  (default = 3 for EEG, 2 for MEG)
%   'backproject'     = 'yes' or 'no', in the case of a rank reduction this parameter determines whether the result will be backprojected onto the original subspace (default = 'yes')
%   'normalize'       = 'no', 'yes' or 'column' (default = 'no')
%   'normalizeparam'  = parameter for depth normalization (default = 0.5)
%   'weight'          = number or Nx1 vector, weight for each dipole position to compensate for the size of the corresponding patch (default = 1)
%
% See also FT_SOURCEANALYSIS, FT_PREPARE_HEADMODEL, FT_PREPARE_SOURCEMODEL

% Copyright (C) 2004-2006, Robert Oostenveld
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
feedback = ft_getopt(varargin, 'feedback', 'text');

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
hasfilter     = false; % not used here
hassubspace   = isfield(sourcemodel, 'subspace');

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

if hasleadfield
  ft_info('using precomputed leadfields\n');
  sourcemodel.leadfield = sourcemodel.leadfield(originside);
else
  ft_info('computing forward model on the fly\n');
end

if hassubspace
  ft_info('using subspace projection\n');
  sourcemodel.subspace = sourcemodel.subspace(originside);
  % remember the original data prior to the voxel dependant subspace projection
  dat_pre_subspace = dat;
end

% Check whether the data is a time series, fourier coefficients or a
% cross-spectral density

if isreal(dat)
  ft_notice('the input consists of real-valued topographies: computing the dipole moments and variance')
  datatype = 'time';
elseif size(dat,1)==size(dat,2) && sum(sum((abs(dat - dat')<10^-10)))==numel(dat)
  ft_notice('the input consists of a cross-spectral density: computing source-level power')
  datatype = 'csd';
else
  ft_notice('the input consists of complex-valued topographies: computing source-level Fourier coefficients and power')
  datatype = 'fourier';
end

% allocate space to hold the result
estimate = [];
estimate.rv  = nan(size(sourcemodel.pos,1),1);
estimate.pow = nan(size(sourcemodel.pos,1),1);
estimate.mom = cell(size(sourcemodel.pos,1),1);

ft_progress('init', feedback, 'scanning grid');
for i=1:size(sourcemodel.pos,1)
  ft_progress(i/size(sourcemodel.pos,1), 'scanning grid %d/%d\n', i, size(sourcemodel.pos,1));
  
  if hasleadfield && hasmom && size(sourcemodel.mom, 1)==size(sourcemodel.leadfield{i}, 2)
    % reuse the leadfield that was previously computed and project
    lf = sourcemodel.leadfield{i} * sourcemodel.mom(:,i);
  elseif  hasleadfield &&  hasmom
    % reuse the leadfield that was previously computed but don't project
    lf = sourcemodel.leadfield{i};
  elseif  hasleadfield && ~hasmom
    % reuse the leadfield that was previously computed
    lf = sourcemodel.leadfield{i};
  elseif ~hasleadfield &&  hasmom
    % compute the leadfield for a fixed dipole orientation
    lf = ft_compute_leadfield(sourcemodel.pos(i,:), sens, headmodel, leadfieldopt{:}) * sourcemodel.mom(:,i);
  else
    % compute the leadfield
    lf = ft_compute_leadfield(sourcemodel.pos(i,:), sens, headmodel, leadfieldopt{:});
  end
  
  if hassubspace
    % do subspace projection of the forward model leadfield matrix
    lf = sourcemodel.subspace{i} * lf;
    
    % the data and the covariance (or cross-spectral density) become voxel dependent due to the projection
    if strcmp(datatype, 'time') || strcmp(datatype, 'csd')
      dat = sourcemodel.subspace{i} * dat_pre_subspace*sourcemodel.subspace{i}'; % Subspace cross-spectral density
    else
      dat = sourcemodel.subspace{i} * dat_pre_subspace; % Subspace time-series or fourier coefficients
    end
  end
  
  % Projection matrix (Pseudoinverse)
  lfi = pinv(lf);
    
  if strcmp(datatype, 'time') || strcmp(datatype, 'fourier')
    % Compute dipole moment and residual variance
    estimate.mom{i} = lfi * dat;
    estimate.rv(i)  = sum(sum(abs(dat - lf*estimate.mom{i}).^2, 1), 2)./sum(sum(abs(dat).^2, 1), 2);
    % Compute power at each location, this is convenient for plotting
    estimate.pow(i) = mean(sum(abs(estimate.mom{i}(:)).^2, 1));  % FIXME is this normalization correct?
  else
    % Compute power, the data represents a covariance or CSD matrix
    estimate.pow(i) = sum(real(sum((lfi*dat).*lfi,2)));
    % Compute residual power (variance)
    Prj = eye(size(lf,1)) - lf*lfi; % Projector to the orthogonal complement of the model space
    estimate.rv(i)  = sum(real(sum((Prj*dat).*Prj,2)));
  end

end % for each dipole position
ft_progress('close');

% reassign the estimated values over the inside and outside grid positions
estimate.inside  = originside;
estimate.pos     = origpos;
if isfield(sourcemodel, 'mom')
  estimate.mom( originside) = estimate.mom;
  estimate.mom(~originside) = {[]};
end
estimate.rv( originside) = estimate.rv;
estimate.rv(~originside) = nan;
estimate.pow( originside) = estimate.pow;
estimate.pow(~originside) = nan;

