function opto = snirf2opto(probe, measurementList)

% SNIRF2OPTO converts the SNIRF probe and measurementList structures to a FieldTrip
% optode structure.
%
% See https://github.com/fNIRS/snirf/blob/master/snirf_specification.md
%
% The FieldTrip optode structure is defined in FT_DATATYPE_SENS
%
% See also OPTO2HOMER, BTI2GRAD, CTF2GRAD, FIF2GRAD, ITAB2GRAD, MNE2GRAD, NETMEG2GRAD, YOKOGAWA2GRAD, FT_DATATYPE_SENS

% Copyright (C) 2020, Robert Oostenveld
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

nSrc = length(probe.sourceLabels);
nDet = length(probe.detectorLabels);

opto = [];
opto.optolabel = cat(1, probe.sourceLabels(:), probe.detectorLabels(:)); % ensure that they are columns

if ~isempty(probe.sourcePos3D)
  % use the 3D positions if available
  opto.optopos = cat(1, probe.sourcePos3D, probe.detectorPos3D);
else
  % use the 2D positions
  opto.optopos = cat(1, probe.sourcePos2D, probe.detectorPos2D);
  opto.optopos(:,3) = 0; % convert the positions to 3D
end

opto.optotype = cat(1, repmat({'transmitter'}, nSrc, 1), repmat({'receiver'}, nDet, 1));
opto.wavelength = probe.wavelengths(:)';

M = length(measurementList);    % number of channels
N = nSrc+nDet;                  % number of optodes
K = numel(probe.wavelengths);   % number of wavelengths

opto.label = cell(M,1);
for i=1:M
  tx = measurementList(i).sourceIndex;                        % transmitter
  rx = measurementList(i).detectorIndex;                      % receiver
  wl = probe.wavelengths(measurementList(i).wavelengthIndex); % wavelength in nm
  opto.label{i} = sprintf('%s-%s [%dnm]', probe.sourceLabels{tx}, probe.detectorLabels{rx}, round(wl));
end

% the following specifies for each of the M channels at which wavelength each of the
% N optodes transmits (positive integer from 1 to K), or receives (negative ingeger
% from 1 to K), or does not contribute at all (zeros)

opto.tra = zeros(M, N);
for i=1:M
  tx = measurementList(i).sourceIndex;
  rx = measurementList(i).detectorIndex + nSrc; % the transmitters are first in the list of optodes
  opto.tra(i, tx) = +measurementList(i).wavelengthIndex;
  opto.tra(i, rx) = -measurementList(i).wavelengthIndex;
end

