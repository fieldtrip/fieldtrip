function event = read_brainvision_vmrk(filename)

% READ_BRAINVISION_VMRK reads the markers and latencies
% it returns the stimulus/response code and latency in ms.
%
% Use as
%   event = read_brainvision_vmrk(filename)
%
% See also READ_BRAINVISION_VHDR, READ_BRAINVISION_EEG

% Copyright (C) 2003, Michael Schulte
% Copyright (C) 2003-2021 Robert Oostenveld
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


fid = fopen_or_error(filename, 'rt');

event = [];
line  = [];

readTime = ft_platform_supports('datetime');

while ischar(line) || isempty(line)
  line = fgetl(fid);
  if ~isempty(line) && ~(isnumeric(line) && line==-1)
    if startsWith(line, 'Mk')
      % this line starts with "Mk", so it probably contains a marker
      tok = split(line, '=');
      if length(tok)~=2
        ft_warning('skipping unexpected formatted line in BrainVision marker file');
      else
        % the line looks like "MkXXX=YYY", which is ok
        % the interesting part now is in the YYY, i.e. the second token
        tok = split(tok{2}, ',');
        if isempty(tok{1})
          tok{1}  = [];
        end
        if isempty(tok{2})
          tok{2}  = [];
        end
        event(end+1).type     = tok{1};
        event(end  ).value    = tok{2};
        event(end  ).sample   = str2double(tok{3});
        event(end  ).duration = str2double(tok{4});
        if numel(tok)>5 && readTime
          try
            event(end).timestamp = datetime(tok{6}, 'InputFormat', 'yyyyMMddHHmmssSSSSSS');
          catch
            ft_warning('skipping invalid datetime in BrainVision marker file');
            readTime = false;
          end
        end
      end
    end
  end
end

fclose(fid);
