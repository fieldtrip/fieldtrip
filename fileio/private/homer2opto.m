function sens = homer2opto(SD)

% HOMER2OPTO converts the Homer SD structure to a FieldTrip optode structure
%
% The Homer SD structure contains
%          Lambda: [780 850]
%          SrcPos: [16x3 double]
%          DetPos: [16x3 double]
%        DummyPos: []
%           nSrcs: 16
%           nDets: 16
%         nDummys: 0
%        MeasList: [76x4 double]
%      SpringList: []
%      AnchorList: {}
%          SrcMap: [2x16 double]
%           vrnum: {'1.0'  '2'}
%            xmin: -14.5263
%            xmax: 254.5263
%            ymin: -14.5263
%            ymax: 224.5263
%     MeasListAct: [76x1 double]
%     MeasListVis: [76x1 double]
%     SpatialUnit: 'mm'
%
% where
%     SD.MeasList  column 1 = transmitter index
%     SD.MeasList  column 2 = receiver index
%     SD.MeasList  column 3 = ?
%     SD.MeasList  column 4 = wavelength index
%
% The FieldTrip optode structure is defined in FT_DATATYPE_SENS
%
% See also BTI2GRAD, CTF2GRAD, FIF2GRAD, ITAB2GRAD, MNE2GRAD, NETMEG2GRAD, YOKOGAWA2GRAD

% Copyright (C) 2004-2016, Robert Oostenveld
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

sens.optopos    = cat(1, SD.SrcPos, SD.DetPos);
sens.optotype   = cat(1, repmat({'transmitter'}, [SD.nSrcs, 1]), repmat({'receiver'}, [SD.nDets, 1]));
sens.wavelength = SD.Lambda;

% FIXME positions seem to be in degrees, not in mm

if isfield(SD, 'SpatialUnit')
  % not sure whether this is always present
  sens.unit = SD.SpatialUnit;
end

% laser strength is not known, this probably does not apply to diode based systems anyway
sens.laserstrength = nan(1,K);

sens.transmits = false(N, K);
% the optodes are sorted: first transmitters, then receivers
for transmitter=1:SD.SrcPos
  sel = SD.MeasList(:,1)==transmitter;
  for wavelength=1:numel(SD.Lambda)
    if any(SD.MeasList(sel,4)==wavelength)
      sens.transmits(transmitter,wavelength) = true;
    end
  end % for all wavelengths
end % for all transmitters

sens.tra = false(N, M);
for chan=1:M
  transmitter = SD.MeasList(chan, 1);
  receiver    = SD.MeasList(chan, 2);
  sens.tra(transmitter, chan) = true;
  sens.tra(receiver,    chan) = true;
end
