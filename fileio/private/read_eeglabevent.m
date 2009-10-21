% read_eeglabevent() - import EEGLAB dataset events
%
% Usage:
%   >> header = read_eeglabevent(filename);
%
% Inputs:
%   filename - [string] file name
%
% Optional inputs:
%   'header' - FILEIO structure header
%
% Outputs:
%   event    - FILEIO toolbox event structure
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2008-

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2008 Arnaud Delorme, SCCN, INC, UCSD, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: read_eeglabevent.m,v $
% Revision 1.4  2009/08/29 05:11:31  josdie
% After consultation with Arno, changed so that value (when present) and type fields in EEGlab .set files are treated as being reversed from value and type fields in FieldTrip files.
%
% Revision 1.3  2009/08/09 01:45:24  josdie
% Changed event value to equal EEGlab's event value rather than type.
%
% Revision 1.2  2009/02/02 20:45:34  josdie
% FieldTrip's .value field now set to EEGlab's .type field.  FieldTrip's .type field set to 'trigger'.  FieldTrip's .duration field set to 0 rather than empty.
%
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.1  2008/04/18 14:04:48  roboos
% new implementation by Arno, shoudl be tested
%

function event = read_eeglabevent(filename, varargin)

if nargin < 1
  help read_eeglabheader;
  return;
end;

header    = keyval('header',     varargin);

if isempty(header)
  header = read_eeglabheader(filename);
end;

event = [];
oldevent = header.orig.event;
for index = 1:length(oldevent)
    event(index).value   = num2str( oldevent(index).type );
    if isfield(oldevent,'code')
        event(index).type   = oldevent(index).code;
    elseif isfield(oldevent,'value')
        event(index).type   = oldevent(index).value;
    else
        event(index).type   = 'trigger';
    end;
if header.nTrials > 1
    event(index).sample = oldevent(index).latency-header.nSamplesPre;
    event(index).offset = header.nSamplesPre;
  else
    event(index).sample = oldevent(index).latency;
    event(index).offset = 0;
  end;
  if isfield(oldevent, 'duration')
    event(index).duration = oldevent(index).duration;
  else
    event(index).duration = 0;
  end;
end;
