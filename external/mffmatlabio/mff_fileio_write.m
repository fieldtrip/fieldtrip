% mff_fileio_write - export Fieldtrip structures to MFF file
%
% Usage:
%   mff_fileio_write(mffFile, header, data, event);
%
% Inputs:
%   mffFile - [string] MFF file name
%   header   - fieldtrip data header 
%   data     - fieldtrip raw data
%   event    - fieldtrip event structure (optional)

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

function mff_fileio_write(filename, hdr, data, evt)

if exist('pop_fileio')
    a = dir(which('pop_fileio.m'));
    b = dir(which('pop_fileio2.m'));
    if isempty(a) || isempty(b) || a.bytes ~= b.bytes
        fprintf([ '-----------------------------------------------------------------------------\n' ...
                  'WARNING: pop_fileio.m and pop_fileio2.m differ; likely you are using\n' ...
                  '         an old version of EEGLAB; The importer will try to import data using\n' ...
                  '         the pop_fileio.m function, and revert to pop_fileio2.m if it crashes\n' ...
                  '-----------------------------------------------------------------------------\n' ]);
    end    
    try
        EEG = pop_fileio(hdr, data, evt); % maybe the wrong version?
    catch
        EEG = pop_fileio2(hdr, data, evt); % backup function included here
    end   
else
    EEG = pop_fileio2(hdr, data, evt); % backup function included here
end    

if isfield(hdr, 'orig')
    if isfield(hdr.orig, 'etc')
        EEG.etc = hdr.orig.etc;
    end
end

mff_export(EEG, filename);

% add code for eeg_emptyset
% add code for pop_fileio
% add code for finputcheck
