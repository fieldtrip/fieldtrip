% mff_importkeys - import MFF keys into structure
%
% Usage:
%   events = mff_importkeys(events, eventind, keys);
%
% Inputs:
%  events   - EEGLAB event structure
%  eventind - event index
%  keys     - MFF keys (in the form of a Java array)
%
% Output:
%  events   - EEGLAB event structure

% This file is part of mffmatlabio.
%
% mffmatlabio is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% mffmatlabio is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with mffmatlabio.  If not, see <https://www.gnu.org/licenses/>.

function events = mff_importkeys(events, eventCount, keylist, eeglabExport)

if nargin < 3
    error('mff_importkeys requires 3 arguments');
end
if nargin < 4
    eeglabExport = true;
end

if isempty(keylist)
    return
end
eventkeycount = keylist.size;
events(eventCount).mffkeys = char(keylist);

keyVals = [];
for q = 0:eventkeycount-1
    theKey = keylist.get(q);
    keyVals(q+1).code = char(theKey.getCode);
    keyVals(q+1).data = char(theKey.getData);
    keyVals(q+1).datatype = char(theKey.getDataType);
    keyVals(q+1).description = char(theKey.getDescription);
    cleancode = keyVals(q+1).code;
    cleancode( cleancode < 48 ) = []; % invalid char
    cleancode( cleancode > 57 & cleancode < 64 ) = []; % invalid char
    try
        events(eventCount).( [ 'mffkey_' cleancode ]) = keyVals(q+1).data;
    catch
        if showWarning
            disp('Warning: issue when converting MFF event key ************');
            showWarning = false;
        end
    end
end
if eeglabExport
    events(eventCount).mffkeysbackup = vararg2str(keyVals); % for exporting
end