function event = read_spmeeg_event(filename, varargin)

% read_spmeeg_event() - import evtns from SPM5 and SPM8 meeg datasets
%
% Usage:
%   >> header = read_spmeeg_event(filename);
%
% Inputs:
%   filename - [string] file name
%
% Outputs:
%   event   - FILEIO toolbox event structure
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% Vladimir Litvak

if nargin < 1
    help read_spmeeg_event;
    return;
end;

header = ft_getopt(varargin, 'header');

if isempty(header)
    header = read_spmeeg_header(filename);
end;

D = header.orig;

event = [];

if isfield(D, 'Radc') % SPM5
    for i = 1:D.Nevents
        if isfield(D, 'events') && isfield(D.events, 'code') && length(D.events.code) == D.Nevents
            value = D.events.code(i);
        else
            value = [];
        end

        event = [event struct('type', 'trial', 'sample', (i-1)*header.nSamples + 1,...
            'value', value, 'offset', -header.nSamplesPre, 'duration', header.nSamples)];
    end

    if D.Nevents == 1 && isfield(D, 'events') && isfield(D.events, 'time')
        event = [event struct('type', 'spm5_event', 'sample', num2cell(D.events.time),...
            'value', num2cell(D.events.code), 'offset', -header.nSamplesPre, 'duration', header.nSamples)];
    end

elseif all(isfield(D, {'type', 'Nsamples', 'Fsample', 'timeOnset'})) % SPM8

    for i = 1:numel(D.trials)
        if isfield(D.trials, 'events')
            cevent = D.trials(i).events;
            if ~isempty(cevent)
                for j = 1:numel(cevent)
                    if ~strcmp(cevent(j).type, 'trial')
                        csample = round((cevent(j).time - D.trials(i).onset)*header.Fs + 1);
                        if csample > 0 && csample <= header.nSamples
                            tmp = rmfield(cevent(j), 'time');
                            if strcmp(D.type, 'continuous') && (D.timeOnset ~= 0)
                                tmp.sample  = (cevent(j).time - D.timeOnset)* header.Fs + 1;
                            else
                                tmp.sample = (i-1)*header.nSamples + 1 +(cevent(j).time - D.trials(i).onset)*header.Fs;
                            end
                            tmp.sample = round(tmp.sample);
                            tmp.offset = 0;
                            event = [event tmp];
                        end
                    end
                end
            end
        end
    end
    
else
    error('Cannot recognize an SPM EEG header format');
end

