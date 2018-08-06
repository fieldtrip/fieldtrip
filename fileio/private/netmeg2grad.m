function [grad] = netmeg2grad(hdr)

% NETMEG2GRAD converts a NetMEG header to a gradiometer structure
% that can be understood by FieldTrip and Robert Oostenveld's low-level
% forward and inverse routines. This function only works for headers
% that have been read using FT_READ_DATA and NETCDF.

% Copyright (C) 2011 Robert Oostenveld
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

% start with an empty structure
grad = [];

% According to FT_DATATYPE_SENS, the structure for MEG gradiometers and/or magnetometers contains
%      sens.label    = Mx1 cell-array with channel labels
%      sens.chanpos  = Mx3 matrix with channel positions
%      sens.chanori  = Mx3 matrix with channel orientations, used for synthetic planar gradient computation
%      sens.tra      = MxN matrix to combine coils into channels
%      sens.coilpos  = Nx3 matrix with coil positions
%      sens.coilori  = Nx3 matrix with coil orientations
%      sens.balance  = structure containing info about the balancing, See FT_APPLY_MONTAGE

numchan = hdr.orig.Dim.numsensors; 
numcoil = sum(hdr.orig.Var.numelementsinsensor);

% FIXME there are 310 channels in the test file, out of which the first 306
% are MEG gradiometers and magnetometers. I am not sure whether those will
% always be the first channels
grad.label   = hdr.label(1:numchan);

% this is probably a reasonable orientation
grad.chanori = squeeze(hdr.orig.Var.sensorelementsorient(:,1,:));
grad.chanpos = hdr.orig.Var.sensorlocation;

pos = zeros(numcoil, 3);
ori = zeros(numcoil, 3);
tra = zeros(numchan, numcoil);

coil = 1;
for i=1:numchan
  for j=1:hdr.orig.Var.numelementsinsensor(i)
    % hdr.orig.Var.sensorelementradius;
    pos(coil,:) = hdr.orig.Var.sensorelementsloc(i,j,:);
    ori(coil,:) = hdr.orig.Var.sensorelementsorient(i,j,:);
    tra(i,coil) = hdr.orig.Var.coilweight(i,j);
    coil = coil+1;
  end
end

% add the information for all coils
grad.tra     = tra;
grad.coilpos = pos;
grad.coilori = ori;

% FIXME I don't know whether this always applies
grad.unit = 'cm';

% FIXME this might be a neuromag, 4d or ctf systems
grad.type = 'unknown'; 
