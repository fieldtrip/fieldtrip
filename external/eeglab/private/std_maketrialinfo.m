% std_maketrialinfo() - create trial information structure using the 
%                       .epoch structure of EEGLAB datasets
%
% Usage: 
%   >> STUDY = std_maketrialinfo(STUDY, ALLEEG);  
%
% Inputs:
%   STUDY      - EEGLAB STUDY set
%   ALLEEG     - vector of the EEG datasets included in the STUDY structure 
%
% Inputs:
%   STUDY      - EEGLAB STUDY set updated. The fields which is created or
%                updated is STUDY.datasetinfo.trialinfo
%
% Authors: Arnaud Delorme, SCCN/INC/UCSD, April 2010

% Copyright (C) Arnaud Delorme arno@ucsd.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function STUDY = std_maketrialinfo(STUDY, ALLEEG);

%% test if .epoch field exist in ALLEEG structure
epochfield = cellfun(@isempty, { ALLEEG.epoch });
if any(epochfield)
    fprintf('Warning: some datasets are continuous and trial information cannot be created\n');
    return;
end

%% check if conversion of event is necessary
ff = {};
flagConvert = true;
for index = 1:length(ALLEEG), ff = union(ff, fieldnames(ALLEEG(index).event)); end
for iField = 1:length(ff)
    
    fieldChar = zeros(1,length(ALLEEG))*NaN;
    for index = 1:length(ALLEEG)
        if isfield(ALLEEG(index).event, ff{iField}) 
            if ischar(ALLEEG(index).event(1).(ff{iField}))
                 fieldChar(index) = 1;
            else fieldChar(index) = 0;
            end
        end
    end
    if ~all(fieldChar(~isnan(fieldChar)) == 1) && ~all(fieldChar(~isnan(fieldChar)) == 0)
        % need conversion to char
        for index = 1:length(ALLEEG)
            if fieldChar(index) == 0
                if flagConvert, disp('Warning: converting some event fields to strings - this may be slow'); flagConvert = false; end
                for iEvent = 1:length(ALLEEG(index).event)
                    ALLEEG(index).event(iEvent).(ff{iField}) = num2str(ALLEEG(index).event(iEvent).(ff{iField}));
                end
            end
        end
    end
end
                
%% Make trial info
for index = 1:length(ALLEEG)
    tmpevent = ALLEEG(index).event;
    eventlat = abs(eeg_point2lat( [ tmpevent.latency ], [ tmpevent.epoch ], ALLEEG(index).srate, [ALLEEG(index).xmin ALLEEG(index).xmax]));
    events   = ALLEEG(index).event;
    ff = fieldnames(events);
    ff = setdiff_bc(ff, { 'latency' 'urevent' 'epoch' });
    trialinfo = [];
    
    % process time locking event fields
    % ---------------------------------
    indtle    = find(eventlat == 0);
    epochs    = [ events(indtle).epoch ];
    extractepoch = true;
    
    % Double checking and changing threshold
    if length(epochs) < ALLEEG(index).trials
        indtle    = find(eventlat < 0.02);
        epochs    = [ events(indtle).epoch ];
    end
    
    if length(epochs) ~= ALLEEG(index).trials
        % special case where there are not the same number of time-locking
        % event as there are data epochs
        if length(unique(epochs)) ~= ALLEEG(index).trials
            extractepoch = false;
            disp('std_maketrialinfo: not the same number of time-locking events as trials, trial info ignored');
        else
            % pick one event per epoch
            [tmp tmpind] = unique_bc(epochs(end:-1:1)); % reversing the array ensures the first event gets picked
            tmpind = length(epochs)+1-tmpind;
            indtle = indtle(tmpind);
            if length(indtle) ~= ALLEEG(index).trials
                extractepoch = false;
                disp('std_maketrialinfo: not the same number of time-locking events as trials, trial info ignored');
            end
        end
    end
    if extractepoch
        commands = {};
        for f = 1:length(ff)
            eval( [ 'eventvals = {events(indtle).' ff{f} '};' ]);
            %if isnumeric(eventvals{1})
                %eventvals = cellfun(@num2str, eventvals, 'uniformoutput', false);
            %    eventvals = [ eventvals{:} ];
            %end
            commands = { commands{:} ff{f} eventvals };
        end
        trialinfo = struct(commands{:});
        STUDY.datasetinfo(index).trialinfo = trialinfo;
    end
    
%    % same as above but 10 times slower
%     for e = 1:length(ALLEEG(index).event)
%         if eventlat(e) < 0.0005 % time locking event only
%             epoch = events(e).epoch;
%             for f = 1:length(ff)
%                 fieldval  = getfield(events, {e}, ff{f});
%                 trialinfo = setfield(trialinfo, {epoch}, ff{f}, fieldval);
%             end
%         end
%     end
end

    
