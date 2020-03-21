%	
%   Read MEF3 meta- and signaldata from a session directory
%	
%   [metadata, data] = readMef3(sessPath, password, channels, rangeType, rangeStart, rangeEnd)
%	
%   sessPath        = path (absolute or relative) to the MEF3 session folder
%   password        = password to the MEF3 data; Pass empty string/variable if not encrypted
%   channels        = a cell array with the names of the channels to return the signal data from. The input
%                     order here will determine the order of channels in the output matrix. If left empty, all
%                     channels will be read abd returned in the order in which they occur in the metadata
%   rangeType       = (optional) modality that is used to define the data-range to read [either 'time' or 'samples']
%   rangeStart      = (optional) start-point for the reading of data (either as a timepoint or samplenumber), pass -1
%                     to start at the first sample of the timeseries
%   rangeEnd        = (optional) end-point to stop the of reading data (either as a timepoint or samplenumber), pass -1
%                     to end at the last sample of the timeseries
%
%   Returns:
%       metadata    = A structing that contains all session/channel/segment metadata
%       data        = A matrix of doubles containing the requested channel(s) signal data. The first
%                     dimension (rows) represents the samples, the second dimension (columns) represents
%                     the channels (in the order of occurance in the 'channels' input argument or metadata)
%
%
%
%   Examples:
%
%       [metadata] = readMef3('./mefSessionData/');                                          % read metadata only  
%       [metadata, signaldata] = readMef3('./mefSessionData/', [], {'Ch02', 'Ch07'});        % read metadata and two channels of data
%       [metadata, signaldata] = readMef3('./mefSessionData/', [], [], 'samples', 0, 1000);  % read all channels, samples 0-1000  
% 
%
%   Note:  When the rangeType is set to 'samples', the function simply returns the samples as they are
%          found (consecutively) in the datafile, without any regard for time or data gaps; Meaning 
%          that, if there is a time-gap between samples, then these will not appear in the result returned.
%          In contrast, the 'time' rangeType will return the data with NaN values in place for the missing samples.
%
%
%   Copyright 2020, Max van den Boom (Multimodal Neuroimaging Lab, Mayo Clinic, Rochester MN)
%   

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [metadata, data] = readMef3(sessPath, password, channels, rangeType, rangeStart, rangeEnd)
    metadata = [];    
    data = [];
    
    % make sure the password argument has a value
    if ~exist('password', 'var') || isempty(password),  password = [];   end
    
    % make sure the session directory is valid and exists
    if ~exist('sessPath', 'var') || isempty(sessPath) || ~ischar(sessPath)
        fprintf(2, ['Error: missing or invalid session directory ''', sessPath, '''\n']);
        return;
    end
    if ~exist(sessPath, 'dir')
        fprintf(2, ['Error: session directory ''', sessPath, ''' could not be found\n']);
        return;
    end
    
    % read all the metadata in the session (including channels and segments)
    try
        metadata = read_mef_session_metadata(sessPath, password);
    catch e
        fprintf(2, [e.message, '\nUnable to read MEF3 metadata\n']);
        return;
    end
    
    % check whether any (meta)data was found
    % note: kind of weak check, but the underlying meflib doesn't give
    % us much more to work with in terms of error handling
    if isempty(metadata) || ...
        ~isfield(metadata, 'earliest_start_time') || metadata.earliest_start_time < 0 || ...
        ~isfield(metadata, 'latest_end_time') || metadata.latest_end_time < 0
        
        fprintf(2, ['Error: no MEF3 (meta)data found in directory ''', sessPath, '''\n']);
        return;
        
    end
    
    % warn if there are no channels
    if metadata.number_of_time_series_channels == 0 && metadata.number_of_video_channels == 0
        warning('on'); warning('backtrace', 'off');
        warning('No channels found in session directory');
    end
    
    % check if signal data should be returned
    if nargout > 1
        
        % check the range type input argument
        if ~exist('rangeType', 'var'),  rangeType = 'samples';  end
        if isempty(rangeType) || ~ischar(rangeType)
            fprintf(2, 'Error: invalid rangeType input argument, should be either ''time'' or ''samples''\n');
            return;
        end
        rangeType = lower(rangeType);
        if ~strcmp(rangeType, 'time') && ~strcmp(rangeType, 'samples')
            fprintf(2, 'Error: invalid rangeType input argument, should be either ''time'' or ''samples''\n');
            return
        end
        
        % check the range start input argument
        if ~exist('rangeStart', 'var'), rangeStart = -1;        end
        if isempty(rangeStart) || ~isnumeric(rangeStart) || length(rangeStart) ~= 1 || (~(rangeStart == -1) && rangeStart < 0)
            fprintf(2, 'Error: invalid rangeStart input argument, should be a single value numeric (either -1 or >= 0)\n');
            return;
        end
        
        % check the range end input argument
        if ~exist('rangeEnd', 'var'),   rangeEnd = -1;          end
        if isempty(rangeEnd) || ~isnumeric(rangeEnd) || length(rangeEnd) ~= 1 || (~(rangeEnd == -1) && rangeEnd < 0)
            fprintf(2, 'Error: invalid rangeEnd input argument, should be a single value numeric (either -1 or >= 0)\n');
            return;
        end
        
        % make sure the channels argument has a value
        if ~exist('channels', 'var') || isempty(channels),  channels = {};  end
            
        % if no channels requested, then request all
        if isempty(channels)
            for i = 1:metadata.number_of_time_series_channels
                channels{end + 1} = metadata.time_series_channels(i).name;
            end
        end
        
        % allow the request of a single channel as a string argument
        if ischar(channels), channels = {channels}; end;
        
        % check the channels input argument
        if ~iscell(channels)
            fprintf(2, ['Error: invalid input argument for ''channels'', should be a cell array containing channel names as string (e.g. {''ch1'', ''ch2'', ''ch3''})\n']);
            return;
        end
        for i = 1:length(channels)
            if ~ischar(channels{i})
                fprintf(2, ['Error: invalid input argument for ''channels'', should be a cell array containing channel names as string (e.g. {''ch1'', ''ch2'', ''ch3''})\n']);
                return;
            end 
        end
        
        % make sure all requested channels exist
        % note: if one or more channels are not found will return an error. Channel selection can be sensitive;
        %       This way prevents unexpected consequences that could arise when instead returning less or empty channels
        channelsFound = ismember(lower(channels), lower({metadata.time_series_channels.name}));
        if sum(channelsFound) < length(channels)
            for i = 1:length(channels)
                if channelsFound(i) == 0
                    fprintf(2, ['Error: requested channel ''', channels{i}, ''' was not found\n']);
                end
            end
            return;
        end
        
        % catch errors
        try
            
            % loop through the channels in the order they are requested
            for i = 1:length(channels)
                channelIndex = find(ismember(lower({metadata.time_series_channels.name}), lower(channels{i})));
                channelPath = [metadata.time_series_channels(channelIndex).path, filesep, metadata.time_series_channels(channelIndex).name, '.', metadata.time_series_channels(channelIndex).extension];

                % read the data of each channel
                % note: we cannot beforehand determine the size of matrix because of potential recording
                %       gaps in the data. Therefore, we cannot allocate memory beforehand. We do assume the
                %       signal sample length is equal over channels. So instead we use horzcat to expand the
                %       matrix linearly in memory (in contrast to vertcat, which would internally require
                %       copying each element and is therefore much slower)
                if i == 1
                    data = read_mef_ts_data(channelPath, password, rangeType, rangeStart, rangeEnd)';
                else
                    data = horzcat(data, read_mef_ts_data(channelPath, password, rangeType, rangeStart, rangeEnd)');
                end

            end
            
        catch e
            
            % error message
            fprintf(2, [e.message, '\nUnable to read MEF3 channel data\n']);
            return;
            
        end
        
    end
    
end
