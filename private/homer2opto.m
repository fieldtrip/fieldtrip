function opto = homer2opto(SD)

% HOMER2OPTO converts the Homer SD structure to a FieldTrip optode structure
%
% See https://www.nitrc.org/plugins/mwiki/index.php/homer2:Homer_Input_Files#NIRS_data_file_format
%
% The Homer SD structure contains the source/detector geometry and has the following fields:
%
% nSrcs    - Number of lasers; scalar variable
% nDets    - Number of detectors; scalar variable
% SrcPos   - Array of probe coordinates of the lasers; dimensions <number of lasers> by 3
% DetPos   - Array of probe coordinates of the detectors; dimensions <number of detectors> by 3
% Lambda   - Wavelengths used for data acquisition; dimensions <number of wavelengths> by 1
% MeasList - List of source/detector/wavelength measurement channels. It’s an array with dimensions, <number of channels> by 4.The meaning of the 4 columns are as follows:
%   Column 1 index of the source from the SD.SrcPos list.
%   Column 2 index of the detector from the SD.DetPos list.
%   Column 3 is unused right now and contains all ones.
%   Column 4 index of the wavelength from SD.Lambda.
%
% The FieldTrip optode structure is defined in FT_DATATYPE_SENS
%
% See also OPTO2HOMER, BTI2GRAD, CTF2GRAD, FIF2GRAD, ITAB2GRAD, MNE2GRAD, NETMEG2GRAD, YOKOGAWA2GRAD, FT_DATATYPE_SENS

% Copyright (C) 2004-2020, Robert Oostenveld
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

M = size(SD.MeasList,1);  % number of channels
N = SD.nSrcs + SD.nDets;  % number of optodes
K = numel(SD.Lambda);     % number of wavelengths

% FIXME positions seem to be in degrees or so, not in mm or cm
if isfield(SD, 'SrcPos3D') && isfield(SD, 'DetPos3D')
  % use the 3D positions if available
  opto.optopos = cat(1, SD.SrcPos3D, SD.DetPos3D);
else
  % use the 2D positions
  opto.optopos = cat(1, SD.SrcPos, SD.DetPos);
end

% give each transmitter and receiver a unique name
opto.optolabel = {};
for i=1:SD.nSrcs
  opto.optolabel{end+1} = sprintf('S%d', i);
end
for i=1:SD.nDets
  opto.optolabel{end+1} = sprintf('D%d', i);
end

opto.optotype   = cat(1, repmat({'transmitter'}, [SD.nSrcs, 1]), repmat({'receiver'}, [SD.nDets, 1]));
opto.wavelength = SD.Lambda(:)'; % this is an 1xK row vector

if isfield(SD, 'SpatialUnit')
  % not sure whether this is always present
  opto.unit = SD.SpatialUnit;
end

opto.label = cell(M,1);
% use the transmitter and receiver numbers and the wavelength to form the the channel names
for i=1:M
  tx = SD.MeasList(i,1);            % transmitter
  rx = SD.MeasList(i,2);            % receiver
  wl = SD.Lambda(SD.MeasList(i,4)); % wavelength in nm
  opto.label{i} = sprintf('S%d-D%d [%dnm]', tx, rx, round(wl));
end

% the following specifies for each of the M channels at which wavelength each of the
% N optodes transmits (positive integer from 1 to K), or receives (negative ingeger
% from 1 to K), or does not contribute at all (zeros)

opto.tra = zeros(M, N);
for i=1:M
  tx = SD.MeasList(i, 1);
  rx = SD.MeasList(i, 2) + SD.nSrcs; % the transmitters are first in the list of optodes
  opto.tra(i, tx) = +SD.MeasList(i, 4);
  opto.tra(i, rx) = -SD.MeasList(i, 4);
end
