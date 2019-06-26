% pop_mffexport - export MFF file from EEGLAB structure.
%
% Usage:
%   pop_mffexport( EEG); % pop up menu to select file
%   pop_mffexport( EEG, mffFile); % export file
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

function com = pop_mffexport(EEG)

com = '';
if nargin < 1
    help pop_mffexport;
    return;
end

matVer = ver('MATLAB');
if datenum(matVer.Date) < 735595
    error('This version of Matlab is too old. Use version 2014a or later');
end

if nargin < 2
    % pop up window
    % -------------
    [fileName, filePath] = uiputfile('*', 'Enter an EGI .mff file/folder');
    if fileName(1) == 0, return; end
    fileName = fullfile(filePath, fileName);
end

mff_export(EEG, fileName);
com = sprintf('pop_mffexport(EEG, ''%s'');', fileName);

