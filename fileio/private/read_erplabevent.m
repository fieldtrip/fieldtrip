% read_erplabevent() - import ERPLAB dataset events
%
% Usage:
%     >> event = read_erplabevent(filename, ...);
%
% Inputs:
%     filename - [string] file name
%
% Optional inputs:
%     'header' - FILEIO structure header
%
% Outputs:
%     event - FILEIO toolbox event structure
%
% Modified from read_eeglabevent

%123456789012345678901234567890123456789012345678901234567890123456789012
%
% Copyright (C) 2008 Arnaud Delorme, SCCN, INC, UCSD, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

function event = read_erplabevent(filename, varargin)

if nargin < 1
  help read_erplabheader;
  return;
end

hdr = ft_getopt(varargin, 'header');

if isempty(hdr)
  hdr = read_erplabheader(filename);
end

event    = [];                % these will be the output in FieldTrip format
oldevent = hdr.orig.bindescr;    % these are in ERPLAB format

for index = 1:length(oldevent)
    event(end+1).type     = 'trial';
    event(end  ).sample   = (index-1)*hdr.nSamples + 1;
    event(end  ).value    = oldevent{index};
    event(end  ).offset   = -hdr.nSamplesPre;
    event(end  ).duration =  hdr.nSamples;
end;    
    
    

