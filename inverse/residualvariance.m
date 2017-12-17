function [dipout] = residualvariance(dip, grad, headmodel, dat, varargin)

% RESIDUALVARIANCE scan with a single dipole and computes the RV
% at each grid location.
%
% Use as
%   [dipout] = residualvariance(dip, grad, headmodel, dat, ...)

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

% get the optional settings, or use default value
feedback = ft_getopt(varargin, 'feedback', 'text');

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
if isfield(dip, 'mom')
  dip.mom = dip.mom(:,originside);
end
if isfield(dip, 'subspace')
  dip.subspace = dip.subspace(originside);
end

if isfield(dip, 'leadfield')
  dip.leadfield = dip.leadfield(originside);
end

if isfield(dip, 'subspace')
  % remember the original data prior to the voxel dependant subspace projection
  dat_pre_subspace = dat;
  fprintf('using subspace projection\n');
end

% Check whether the data is a time series, fourier coefficients or a
% cross-spectral density

if isreal(dat)==1
  fprintf('The input is a time series: computing source level time series and variance')
  datatype = 'time';
elseif size(dat,1)==size(dat,2)&&sum(sum((abs(dat - dat')<10^-10)))==numel(dat)
  fprintf('The input is a cross-spectral density: computing source level power')
  datatype = 'csd';
else
  fprintf('The input are fourier coeffiecients: computing source level fourier coefficients and power')
  datatype = 'fourier';
end

rv  = nan(size(dip.pos,1),1);
pow = nan(size(dip.pos,1),1);
mom = cell(size(dip.pos,1),1);

ft_progress('init', feedback, 'computing inverse');
for i=1:size(dip.pos,1)
  
  ft_progress(i/size(dip.pos,1), 'computing inverse %d/%d\n', i, size(dip.pos,1));
  
  if isfield(dip, 'leadfield')
    % reuse the leadfield that was previously computed
    lf = dip.leadfield{i};
  else
    % compute the leadfield
    lf = ft_compute_leadfield(dip.pos(i,:), grad, headmodel);
  end
  
  if isfield(dip, 'subspace')
    % do subspace projection of the forward model
    % Leadfield matrix
    lf = dip.subspace{i} * lf;
    
    % the data and the covariance (or cross-spectral density) become voxel dependent due to the projection
    if strcmp(datatype, 'time') || strcmp(datatype, 'csd') 
      dat = dip.subspace{i} * dat_pre_subspace*dip.subspace{i}'; %Subspace cross-spectral density
    else
      dat = dip.subspace{i} * dat_pre_subspace; %Subspace time-series or fourier coefficients
    end
  end
  
  % Projection matrix (Pseudoinverse) 
  lfi    = pinv(lf); 
    
  if strcmp(datatype, 'time') || strcmp(datatype, 'fourier') 
    % Compute dipole moment and residual variance
    mom{i} = lfi * dat;
    rv(i)  = sum(sum(abs(dat - lf*mom{i}).^2, 1), 2)./sum(sum(abs(dat).^2, 1), 2);
  
    % for plotting convenience also compute power at each location
    % FIXME is this normalization correct?
    pow(i) = mean(sum(abs(mom{i}(:)).^2, 1));
  else
    % Compute power
    pow(i) = sum(real(sum((lfi*dat).*lfi,2)));
    % Compute residual power (variance)
    Prj = eye(size(lf,1)) - lf*lfi; %Projector to the orthogonal complement of the model space
    rv(i)  = sum(real(sum((Prj*dat).*Prj,2)));
  end
end
ft_progress('close');

% wrap it all up, prepare the complete output
dipout.inside  = originside;
dipout.pos     = origpos;

dipout.rv  = nan(size(originside));
dipout.pow = nan(size(originside));

% assign the output data
dipout.rv(originside)  = rv(:);   % ensure that it is a column vector
dipout.pow(originside) = pow(:);  % ensure that it is a column vector

if strcmp(datatype, 'time') || strcmp(datatype, 'fourier')
  dipout.mom = cell(size(originside));
  dipout.mom( originside) = mom;
  dipout.mom(~originside) = {[]};
end
