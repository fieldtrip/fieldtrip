function [allTimestamps allSnippets allIndices] = findEventTimes(NEV, Channel, Units)

% Given a NEV file, an channel number, and a unit number this function
% will return indices corresponding to those spikes, timestamps
% corresponding to those spikes (optional), and snippets corresponding to
% those spikes (optional). 
%
% Use [timestamps snippets indices] = findEventTimes(NEV, Channel, Unit)
% 
% INPUTS
%
%   NEV:          This corresponds to the NEV file the desired data is
%                 being extracted from.
%
%   Channel:      The channel number the data is being extracted for.
%
%   Units:        The unit numbers the data is being extracted for. This
%                 variable can be one unit or many units passed as an array
%                 of integers.
%                 DEFAULY: If units is not specified all timestamps from
%                 all units will be passed to the calling function. For
%                 noise pass 255.
%
% OUTPUT
%
%   indices:      An array of all indices that correspond to neural data
%                 for channel and unit passed.
%
%   timestamps:   An array of all timestamps that correspond to neural data
%                 for channel and unit passed.
%
%   snippets:     A matrix of all indices that correspond to neural data
%                 for channel and unit passed.
%
%   IF OUTPUT IS NOT SPECIFIED only "timestamps" WILL BE PASSED TO THE
%   CALLING FUNCTION.
%
%   Example: 
%   
%   [timestamps snippets indices] = findEventTimes(NEV, 3, [1,3:5]);
%
%   In the example above, the indices, timestamps, and snippets for
%   channel #3 and units 1, 3, 4, and 5 will be passed to the calling
%   function.
%
%   Kian Torab
%   ktorab@blackrockmicro.com
%   Blackrock Microsystems
%   Salt Lake City, UT
%   Version 2.0.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist('Units', 'var')
    Units = 0:5;
end

if ischar(Channel) && exist('Units', 'var')
    strcmpi(Channel, 'DigitalIn')
        allSnippets = find(NEV.Data.SerialDigitalIO.UnparsedData == Units);
        allTimestamps = NEV.Data.SerialDigitalIO.TimeStamp(allSnippets);
        allIndices = [];
else

    % Find indices that correspond to "Channel"
    ChannelIndices = find([NEV.Data.Spikes.Electrode] == Channel);

    for i = 1:length(Units)
        % Find indices that correspond to "Units" within "Channel" indices
        UnitIndices{i} = find([NEV.Data.Spikes.Unit(ChannelIndices)] == Units(i));

        % Updating the indices so they correspond to the original NEV indices
        indices{i} = ChannelIndices(UnitIndices{i});

        % Finding the timestamps corresponding to the indices
        timestamps{i} = NEV.Data.Spikes.TimeStamp(indices{i});
        % Finding the snippets corresponding to the indices
        if isfield(NEV.Data.Spikes, 'Waveform')
            snippets{i} = NEV.Data.Spikes.Waveform(:,indices{i});
        elseif nargout == 3 && i == 1
            display('Snippet data was not retrieved because the NEV does not contain any snippets.');
            snippets{i} = [];
        end
    end

    allIndices = cell2mat(indices);
    allTimestamps = double(cell2mat(timestamps));
    if exist('snippets', 'var')
        allSnippets = cell2mat(snippets);
    end
end