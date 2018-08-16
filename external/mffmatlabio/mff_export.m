% mff_export - export EEGLAB structure to MFF file. This function calls
%                 all other function to create MFF structure, export
%                 events, channels and channel coordinates.
%
% Usage:
%   mff_export(EEG, mffFile);
%
% Inputs:
%  EEG     - EEGLAB structure
%  mffFile - filename/foldername for the MFF file (MFF file/folder must
%            already exist)

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

function mff_export(EEG, outputFile)

matVer = ver('MATLAB');
if datenum(matVer.Date) < 735595
    error('This version of Matlab is too old. Use version 2014a or later');
end

% add mff extension if not present
[filePath, fileName] = fileparts( outputFile);
outputFile = fullfile(filePath, [ fileName '.mff' ]);

% delete folder if it exist
if exist(outputFile)
    rmdir(outputFile, 's');
end

mff_createmff(outputFile);
mff_exportinfo(EEG, outputFile);
mff_exportsubject(EEG, outputFile);
mff_exportinfon(EEG, outputFile);
mff_exportsignal(EEG, outputFile);
indtle = mff_exportcategories(EEG, outputFile);
EEG.event(indtle) = []; % remove time locking events
mff_exportevents(EEG, outputFile);
mff_exportcoordinates(EEG, outputFile);
mff_exportsensorlayout(EEG, outputFile);
mff_exportpnsset(EEG, outputFile);
mff_exportepochs(EEG, outputFile);
mff_exportsensorlayout(EEG, outputFile);
