% pop_mffimport - import MFF file to EEGLAB structure. Note that events
%                 in MFF files are very rich and that you need to specify 
%                 which field may contain the information for EEGLAB to
%                 extract events. This is the purpose of the second input
%                 "typefield" to that function.
%
% Usage:
%   EEG = pop_mffimport; % pop up menu to select file
%   EEG = pop_mffimport(mffFile, typefield); % import file
%
% Input:
%  mffFile - filename/foldername for the MFF file (MFF file/folder must
%            already exist)
%  typefield - [string or cell] MFF field(s) to use for the type field. The 
%            default is the code field of the MFF file. If several fields
%            are provided, the fields are concatenated to create the EEGLAB
%            type field.
%
% Output:
%  EEG     - EEGLAB structure

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

function [EEG, com] = pop_mffimport(fileName, typefield)

com = '';
matVer = ver('MATLAB');
if datenum(matVer.Date) < 735595
    error('This version of Matlab is too old. Use version 2014a or later');
end

EEG = [];

if ~exist('convertlocs.m', 'file')
    error('This function requires to install the EEGLAB software (and start it)');
end

if nargin < 1
    
    % pop up window
    % -------------
    if ismac
        [fileName, filePath] = uigetfile('*', 'Select an EGI .mff file/folder');
        if fileName(1) == 0, return; end
        fileName = fullfile(filePath, fileName);
    else
        fileName = uigetdir('*', 'Select an EGI .mff file/folder');
        if fileName(1) == 0, return; end
    end
end

EEG = mff_import(fileName);

if nargin < 2
    if isempty(EEG.event)
        disp('No event in data file');
    else
        % popup window parameters
        % -----------------------
        eventfields = setdiff(fieldnames(EEG.event), { 'type' 'latency' 'duration' 'urevent' 'epoch' 'begintime' });
        poscode     = strmatch('code', lower(eventfields), 'exact');
        if isempty(poscode), poscode=1; end
        uilist   = { { 'style' 'text' 'String' ['Event type field:' 10 '(you may select multiple)'] } ...
                     { 'style' 'listbox' 'string' eventfields 'value' poscode 'min' 0 'max' 2} };
        geom = { 1 1 };
        result = inputgui( 'geometry', geom, 'uilist', uilist, 'geomvert', [1 3], 'helpcom', 'pophelp(''pop_mffimport'')', 'title', 'Choose event type field -- pop_mffimport()');
        
        if isempty(result), return; end
        typefield = eventfields(result{1});
    end
end

% Use different event fields to populate the EEGLAB type field
% ------------------------------------------------------------
if ~isempty(typefield) || ~(ischar(typefield) && strcmpi(typefield, 'code'))
    if ischar(typefield), typefield = { typefield }; end
    
    % get data for each MFF event field
    dataField = cell(length(EEG.event), length(typefield));
    for iField = 1:length(typefield)
        dataField(:,iField) = { EEG.event.(typefield{iField}) }';
    end
    
    % copy the data into EEGLAB event structure
    strField = cell(1,length(typefield));
    strField(:) = {'_'}; % add space between events
    for iEvent = 1:length(EEG.event)
        if ~isequal(lower(EEG.event(iEvent).type), 'boundary') % not a boundary
            if ~all(cellfun(@isempty, dataField(iEvent,:))) % not a trial type
                tmpType = [ dataField(iEvent,:); strField ];
                tmpType = [ tmpType{:} ];
                tmpType(end) = [];
                EEG.event(iEvent).type = tmpType;
            end
        end
    end
end
EEG = eeg_checkset(EEG,'eventconsistency');

com = sprintf('EEG = pop_mffimport(''%s'', %s);', fileName, vararg2str({typefield}));

