function SD = opto2homer(opto)

% OPTO2HOMER constructs a Homer-compatible sensor definition (SD) from a FieldTrip
% opto structure.
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
% See also HOMER2OPTO, FT_DATATYPE_SENS

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

[nchan, noptode] = size(opto.tra);
transmitters = find(strcmp(opto.optotype, 'transmitter'));
receivers    = find(strcmp(opto.optotype, 'receiver'));

SD = [];
SD.nSrcs = length(transmitters);
SD.nDets = length(receivers);
SD.SrcPos = opto.optopos(transmitters,:);
SD.DetPos = opto.optopos(receivers,:);
SD.Lambda = opto.wavelength;

SD.MeasList = zeros(nchan, 4);
for i=1:nchan
  optodes = opto.tra(i,:);
  SD.MeasList(i,1) = find(transmitters==find(optodes>0));     % the number of the transmitter
  SD.MeasList(i,2) = find(receivers==find(optodes<0));        % the number of the receiver
  SD.MeasList(i,4) = optodes(optodes>0);                      % the wavelength
end
