function [stimulus, response, segment, timezero] = read_brainvision_vmrk(filename)

% READ_BRAINVISION_VMRK reads the markers and latencies
% it returns the stimulus/response code and latency in ms.
%
% Use as
%   [stim, resp, segment, timezero] = read_brainvision_vmrk(filename)
%
% This function needs to read the header from a separate file and
% assumes that it is located at the same location.
%
% See also READ_BRAINVISION_VHDR, READ_BRAINVISION_EEG

% original M. Schulte 31.07.2003
% modifications R. Oostenveld 14.08.2003
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

stimulus=[];
response=[];
segment=[];
timezero=[];

% read the header belonging to this marker file
hdr=read_brainvision_vhdr([filename(1:(end-4)) 'vhdr']);

fid=fopen_or_error(filename,'rt');

line=1;
while line~=-1,
  line=fgetl(fid);
  % pause
  if ~isempty(line),
    if contains(line,'Mk')
      if contains(line,'Stimulus')
        [token,rem] = strtok(line,',');
        type=sscanf(rem,',S %i');
        [token,rem] = strtok(rem,',');
        time=(sscanf(rem,', %i')-1)/hdr.Fs*1000;
        stimulus=[stimulus; type time(1)];

      elseif contains(line,'Response')
        [token,rem] = strtok(line,',');
        type=sscanf(rem,',R %i');
        [token,rem] = strtok(rem,',');
        time=(sscanf(rem,', %i')-1)/hdr.Fs*1000;
        response=[response; type, time(1)];

      elseif contains(line,'New Segment')
        [token,rem] = strtok(line,',');
        time=(sscanf(rem,',,%i')-1)/hdr.Fs*1000;
        segment=[segment; time(1)];

      elseif contains(line,'Time 0')
        [token,rem] = strtok(line,',');
        time=(sscanf(rem,',,%i')-1)/hdr.Fs*1000;
        timezero=[timezero; time(1)];

      end
    end
  else
    line=1;
  end
end

fclose(fid);
