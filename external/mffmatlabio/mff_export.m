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
if exist('OCTAVE_VERSION', 'builtin') == 0 && datenum(matVer.Date) < 735595
    error('This version of Matlab is too old. Use version 2014a or later');
end

% add mff extension if not present
[filePath, fileName] = fileparts( outputFile);
outputFile = fullfile(filePath, [ fileName '.mff' ]);

% check missing channels and reconstruct data 
if isfield(EEG.etc, 'layout')
    
    if isfield(EEG.chaninfo, 'removedchans') && ~isempty(EEG.chaninfo.removedchans)

            mff_path;
            badchans1 = javaObject('com.egi.services.mff.api.ChannelStatus');
            badchans1.setExclusion('badChannels');
            badchans1.setBinIndex(0);    
            jList1 = javaObject('java.util.ArrayList');
            
            if isfield(EEG.etc, 'info2')
            badchans2 = javaObject('com.egi.services.mff.api.ChannelStatus');
            badchans2.setExclusion('badChannels');
            badchans2.setBinIndex(0);    
            jList2 = javaObject('java.util.ArrayList');
            
            % add back the channels, zero them out and add special structure
            for iRm = 1:length(EEG.chaninfo.removedchans)
                posInsert = str2double(EEG.chaninfo.removedchans(iRm).labels(2:end));
                if isnan(posInsert) % PNS chan
                    posInsert = size(EEG.data,1)+1;
                end
                
                % insert data
                EEG.data(posInsert+1:end+1,:) = EEG.data(posInsert:end,:);
                EEG.data(posInsert,:) = 0;

                % insert channel
                EEG.chanlocs(posInsert+1:end+1) = EEG.chanlocs(posInsert:end);
                EEG.chanlocs(posInsert) = EEG.chaninfo.removedchans(iRm);
                chanint = zeros(1,1,'int32'); 
                
                % deal with PNS channels
                if strcmpi(EEG.chanlocs(posInsert).type, 'PNS')
                    alltypes = { EEG.chanlocs(1:posInsert-1).type };
                    alltypes = cellfun(@(x)num2str(x), alltypes, 'uniformoutput', false);
                    inds = strmatch('pns', lower(alltypes), 'exact');
                    chanint(1) = length(inds)+1; % index of PNS channel
                    jList2.add(chanint);
                else
                    chanint(1) = posInsert;
                    jList1.add(chanint);
                end
            end
    
            % save channel list to remove
            EEG.nbchan = size(EEG.data,1);
            if isfield(EEG.etc, 'info1') && ~isempty(jList1.isEmpty)
                badchans1.setChannels(jList1);
                EEG.etc.info1.ChannelStatus = badchans1; 
            end
            if ~isempty(jList2.isEmpty)
                badchans2.setChannels(jList2);
                if isfield(EEG.etc, 'info2')
                    EEG.etc.info2.ChannelStatus = badchans2; 
                else
                    EEG.etc.info1.ChannelStatus = badchans2; % only PNS
                end
            end
        end
    end
end
           
% delete folder if it exist
if exist(outputFile)
    rmdir(outputFile, 's');
end

mff_createmff(outputFile);
if ~isfield(EEG.etc, 'recordingtime')
    EEG.etc.recordingtime = now;
    EEG.etc.timezone = '00:00';
end
mff_exportinfo(EEG, outputFile);
mff_exportsubject(EEG, outputFile);
mff_exportinfon(EEG, outputFile,1);
if isfield(EEG.etc, 'info2')
    mff_exportinfon(EEG, outputFile, 2);
end
mff_exportsignal(EEG, outputFile);
indtle = mff_exportcategories(EEG, outputFile);
EEG.event(indtle) = []; % remove time locking events
mff_exportevents(EEG, outputFile);
mff_exportcoordinates(EEG, outputFile);
mff_exportsensorlayout(EEG, outputFile);
mff_exportpnsset(EEG, outputFile);
mff_exportepochs(EEG, outputFile);
mff_exportsensorlayout(EEG, outputFile);
disp('Done');