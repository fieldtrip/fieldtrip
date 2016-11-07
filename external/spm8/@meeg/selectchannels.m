function chanind = selectchannels(this, channels)
% Method for getting channel indices based on labels and/or types
% FORMAT  res = selectchannels(this, label)
% this       - MEEG object
% channels   - string or cell array of labels that may also include 
%              'all', or types ('EEG', 'MEG' etc.)
%
% res        - vector of channel indices matching labels
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: selectchannels.m 3798 2010-03-24 12:00:07Z vladimir $

if ischar(channels)
    channels = {channels};
end

chanind = [];

for i = 1:numel(channels)
    switch upper(channels{i})
        case 'ALL'
            chanind = [chanind 1:nchannels(this)];
        case 'EOG'
            chanind = [chanind eogchannels(this)];
        case 'ECG'
            chanind = [chanind ecgchannels(this)];
        case 'EMG'
            chanind = [chanind emgchannels(this)];
        case {'EEG', 'MEG', 'MEGMAG', 'MEGGRAD', 'MEGPLANAR', 'REF', 'REFMAG', 'REFGRAD', 'LFP'}
            chanind = [chanind meegchannels(this, upper(channels{i}))];
        otherwise
            chanind = [chanind indchannel(this, channels{i})];
    end
end

chanind = unique(chanind);