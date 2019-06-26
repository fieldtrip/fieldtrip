function [grad] = nimh2grad(hdr)

% NIMH2GRAD constructs a gradiometer definition from the res4 header whish
% is read using the NIMH implementation of ctf_read_res4. The grad
% structure is compatible with FieldTrip and Robert Oostenveld's low-level
% forward and inverse routines.
%
% Use as
%   hdr  = ctf_read_res4(dataset);
%   grad = nimh2grad(hdr;
%
% See also CTF2GRAD, FIF2GRAD

% Copyright (C) 2005, Robert Oostenveld
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

% only work on the MEG channels
if isfield(hdr.sensor.index, 'meg')
  sel = hdr.sensor.index.meg;
else
  sel = hdr.sensor.index.meg_sens;
end

% start with an empty structure
grad.coilpos = [];
grad.coilori = [];
grad.tra = [];
grad.label = {};

for i=1:length(sel)
  pos = hdr.sensor.info(sel(i)).location';
  ori = hdr.sensor.info(sel(i)).orientation';
  numcoils(i) = size(pos,1);
  if size(ori,1)==1 && size(pos,1)==1
    % one coil position with one orientation: magnetometer
    ori = ori;
  elseif size(ori,1)==1 && size(pos,1)==2
    % two coil positions with one orientation: first order gradiometer
    % assume that the orientation of the upper coil is opposite to the lower coil
    ori = [ori; -ori];
  else
    ft_error('do not know how to deal with higher order gradiometer hardware')
  end

  % add this channels coil positions and orientations
  grad.coilpos = [grad.coilpos; pos];
  grad.ori = [grad.ori; ori];
  grad.label{i} = hdr.sensor.info(sel(i)).label;

  % determine the contribution of each coil to each channel's output signal
  if size(pos,1)==1
    % one coil, assume that the orientation is correct, i.e. the weight is +1
    grad.tra(i,end+1) = 1;
  elseif size(pos,1)==2
    % two coils, assume that the orientation for each coil is correct, i.e. the weights are +1 and +1
    grad.tra(i,end+1) = 1;
    grad.tra(i,end+1) = 1;
  else
    ft_error('do not know how to deal with higher order gradiometer hardware')
  end
end

% prefer to have the labels in a column vector
grad.label = grad.label(:);

% reorder the coils, such that the bottom coils are at the first N
% locations and the top coils at the last N positions. This makes it
% easier to use a selection of the coils for topographic plotting
if all(numcoils==2)
  bot = 1:2:sum(numcoils);
  top = 2:2:sum(numcoils);
  grad.coilpos = grad.coilpos([bot top], :);
  grad.ori = grad.ori([bot top], :);
  grad.tra = grad.tra(:, [bot top]);
end
