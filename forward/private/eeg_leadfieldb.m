function [lf] = eeg_leadfieldb(dippos, elc, vol)

% EEG_LEADFIELDB computes the electric leadfield for a dipole in a volume
% using the boundary element method
%
% Use as
%   [lf] = eeg_leadfieldb(dippos, elc, vol)
% with the input arguments
%   dippos     = position dipole, 1x3 or Nx3
%   elc        = electrode positions, Nx3 (optional, can be empty)
%   vol        = volume conductor model
%
% The volume conductor model is a structure and should have the fields
%   vol.bnd    = structure array with vertices and triangles of each boundary
%   vol.cond   = conductivity for each compartment
%   vol.mat    = system matrix, which can include the electrode interpolation
%
% The compartment boundaries are described by a structure array with
%   vol.bnd(i).pos
%   vol.bnd(i).pos

% Copyright (C) 2003, Robert Oostenveld
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

% do some sanity checks
if ~isfield(vol, 'bnd')
  ft_error('there are no compartment boundaries present');
end

if length(vol.bnd)~=length(vol.cond)
  ft_error('the number of compartments in the volume is inconsistent');
end

if ~isfield(vol, 'mat')
  ft_error('there is no system matrix present');
end

% determine the number of compartments
ncmp = length(vol.bnd);

% the number of rows in the leadfield matrix should either correspond to
% the number of electrodes, to the number of vertices of the skin
% compartment or to the total number of vertices
nelc  = size(elc, 1);
nskin = size(vol.bnd(vol.skin_surface).pos,1);
nall  = 0;
for i=1:ncmp
  nall = nall + size(vol.bnd(i).pos,1);
end
if size(vol.mat,1)==nelc
  % the output leadfield corresponds to the number of electrodes
elseif size(vol.mat,1)==nskin
  % the output leadfield corresponds to the number skin vertices
elseif size(vol.mat,1)==nall
  % the output leadfield corresponds to the total number of vertices
elseif strcmp(ft_voltype(vol), 'openmeeg')
  % this is handled differently, although at the moment I don't know why
else
  ft_error('unexpected size of system matrix')
end

% determine the conductivity of the source compartment
cond = vol.cond(vol.source);

% compute the infinite medium potential on all vertices
switch ft_voltype(vol)
  case 'dipoli'
    % the system matrix was computed using Thom Oostendorp's DIPOLI
    % concatenate the vertices of all compartment boundaries in a single Nx3 matrix
    pos = [];
    for i=1:ncmp
      pos = [pos; vol.bnd(i).pos];
    end
    % dipoli incorporates the conductivity into the system matrix
    lf = inf_medium_leadfield(dippos, pos, 1);
    
  case 'asa'
    % the system matrix was computed using ASA from www.ant-neuro.com
    % concatenate the vertices of all compartment boundaries in a single Nx3 matrix
    pos = [];
    for i=1:ncmp
      pos = [pos; vol.bnd(i).pos];
    end
    % assume that isolated potential approach was used
    lf = inf_medium_leadfield(dippos, pos, cond);
    
  case 'bemcp'
    % the system matrix was computed using code from Christopher Phillips
    cond = [vol.cond 0]; % add the conductivity of air for simplicity
    lf = cell(1,ncmp);
    % loop over boundaries and compute the leadfield for each
    for i=1:ncmp
      co = (cond(i)+cond(i+1))/2 ;
      lf{i} = inf_medium_leadfield(dippos, vol.bnd(i).pos, co);
    end
    % concatenate the leadfields
    lf = cat(1, lf{:});
    
  otherwise
    ft_error('unsupported type of volume conductor (%s)\n', ft_voltype(vol));
end % switch ft_voltype

if isfield(vol, 'mat') && ~ft_voltype(vol, 'openmeeg')
  % compute the bounded medium potential on all vertices
  % this may include the bilinear interpolation from vertices towards electrodes
  lf = vol.mat * lf;
end

