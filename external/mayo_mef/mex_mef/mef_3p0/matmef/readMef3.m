%	
%   Read MEF3 meta- and signaldata from a session directory
%	
%   [metadata, data] = readMef3(sessPath, password, channels, rangeType, rangeStart, rangeEnd)
%   [metadata, data] = readMef3(sessPath, password, channels, rangeType, ranges)
%	
%   sessPath        = path (absolute or relative) to the MEF3 session folder
%   password        = password to the MEF3 data; Pass empty string/variable if not encrypted
%   channels        = a cell array with the names of the channels to return the signal data from. The order of channels
%                     in this input argument will determine the order of channels in the output matrix. If left empty, all
%                     channels will be read and ordered as in the metadata.time_series_channels (ordered according to 
%                     the 'acquisition_channel_number' metadata variable of each channel)
%   rangeType       = (optional) modality that is used to define the data-range to read [either 'time' or 'samples']
%   rangeStart      = (optional) start-point for the reading of data (depending on the rangeType, defined as an epoch/unix
%                     timestamp or samplenumber). Pass -1 to start at the first sample of the timeseries
%   rangeEnd        = (optional) end-point to stop the of reading data (depending on the rangeType, defined as
%                     an epoch/unix timestamp or samplenumber). Pass -1 to end at the last sample of the timeseries
%   ranges          = (optional) a Nx2 matrix of multiple ranges (start- and end-points) for the reading of data.
%
%
%   Returns:
%       metadata    = A structing that contains all session/channel/segment metadata. Will return empty on failure to read
%       data        = A matrix of doubles containing the requested channel(s) signal data. The first dimension (rows) represents
%                     the channels (ordered based on the 'channels' input argument); the second dimension (columns) represents the
%                     the samples/time. If multiple ranges are given then there is a third dimension, representing the requested ranges/epochs
% 
%
%   Examples:
%
%       % single range/epoch
%       [metadata] = readMef3('./mefSessionData/');                                          % read metadata only  
%       [metadata, signaldata] = readMef3('./mefSessionData/', [], {'Ch02', 'Ch07'});        % read metadata and two channels of data
%       [metadata, signaldata] = readMef3('./mefSessionData/', [], [], 'samples', 0, 1000);  % read all channels, samples 0-1000 
%
%       % multiple ranges/epochs
%       ranges = [[0,    1000]; ...
%                 [1000, 2000]; ...
%                 [5000, 6000]];
%       [metadata, signaldata] = readMef3('./mefSessionData/', [], [], 'samples', ranges);
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
function [metadata, data] = readMef3(sessPath, password, channels, rangeType, varargin)
    metadata = [];    
    data = [];
    
    % set defaults
    if ~exist('password', 'var') || isempty(password),  password = [];   end
    if ~exist('channels', 'var') || isempty(channels),  channels = {};  end
    
    % make sure the session directory is valid and exists
    if ~exist('sessPath', 'var') || isempty(sessPath) || ~ischar(sessPath)
        fprintf(2, 'Error: missing or invalid session directory ''%s''\n', sessPath);
        return;
    end
    if ~exist(sessPath, 'dir')
        fprintf(2, 'Error: session directory ''%s'' could not be found\n', sessPath);
        return;
    end

    % read all the metadata in the session (including channels and segments)
    try
        metadata = read_mef_session_metadata(sessPath, password);
    catch e
        metadata = [];
        fprintf(2, '%s\nUnable to read MEF3 metadata\n', e.message);
        return;
    end
    
    % check whether any (meta)data was found
    % note: kind of weak check, but the underlying meflib doesn't give
    % us much more to work with in terms of error handling
    if isempty(metadata) || ...
        ~isfield(metadata, 'earliest_start_time') || metadata.earliest_start_time < 0 || ...
        ~isfield(metadata, 'latest_end_time') || metadata.latest_end_time < 0
        
        fprintf(2, 'Error: no MEF3 (meta)data found in directory ''%s''\n', sessPath);
        metadata = [];
        return;
        
    end
    
    % warn if there are no channels
    if metadata.number_of_time_series_channels == 0 && metadata.number_of_video_channels == 0
        warning('on'); warning('backtrace', 'off');
        warning('No channels found in session directory');
    end
    
    % check if too many input arguments are given
    if nargin > 6
        fprintf(2, 'Error: too many input arguments\n');
        return;
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
        
        % set range defaults
        if ~exist('rangeStart', 'var'), rangeStart = -1;        end
        if ~exist('rangeEnd', 'var'),   rangeEnd = -1;          end
        ranges = [-1, -1];
        
        % check whether a single range or multiple ranges are given
        if nargin == 5
            % 5 input args, matrix with multiple ranges given

            % retrieve the ranges
            ranges = varargin{1};
            
            % check the range start input argument
            if isempty(ranges) || ~isnumeric(ranges) || size(ranges, 2) ~= 2 || size(ranges, 1) < 1 || any(ranges(:) < 0)
                fprintf(2, 'Error: invalid ''ranges'' input argument, should be a Nx2 matrix with numeric values (>= 0)\n');
                return;
            end
            
            % sort the ranges and calculate their lengths
            ranges = sort(ranges, 2, 'asc');
            rangesDiff = ranges(:, 2) - ranges(:, 1);
            
            % check whether any of the ranges results in a length of 0
            if any(rangesDiff <= 0)
                fprintf(2, 'Error: invalid ''ranges'' input argument, one or more ranges result in a length of 0\n');
                return;
            end
            
            % warn if input range differs in length, data will be padded with nans at the end
            if any(rangesDiff ~= rangesDiff(1))
                warning('on'); warning('backtrace', 'off');
                warning('The ranges that were requested differ in length, data will be padded with nans at the end');
            end
            
        elseif nargin == 6
            % 6 input args, range start and end given
            
            % check the range start input argument
            if isempty(rangeStart) || ~isnumeric(rangeStart) || length(rangeStart) ~= 1 || (~(rangeStart == -1) && rangeStart < 0)
                fprintf(2, 'Error: invalid rangeStart input argument, should be a single value numeric (either -1 or >= 0)\n');
                return;
            end

            % check the range end input argument
            if isempty(rangeEnd) || ~isnumeric(rangeEnd) || length(rangeEnd) ~= 1 || (~(rangeEnd == -1) && rangeEnd < 0)
                fprintf(2, 'Error: invalid rangeEnd input argument, should be a single value numeric (either -1 or >= 0)\n');
                return;
            end

            % transfer the start and end of the range
            ranges = [varargin{1}, varargin{2}];
            
            % if there are no -1's, % order and check length
            if ~any(ranges == -1)
                ranges = sort(ranges, 2, 'asc');
                if ranges(2) - ranges(1) <= 0
                    fprintf(2, 'Error: invalid range, the requested range length is 0\n');
                    return;
                end
            end
            
        end
    end

    
    % 
    % sort the channels in the metadat (using channel->metadata->section_2->acquisition_channel_number)
    % 

    % list the acquisition channel numbers
    acqChNum = [];
    for iChannel = 1:metadata.number_of_time_series_channels
        acqChNum(iChannel) = metadata.time_series_channels(iChannel).metadata.section_2.acquisition_channel_number;
    end

    % sort the channels
    [ordAcqChNum, prevIndex] = sort(acqChNum);

    % check if it starts at one
    if min(acqChNum) ~= 1
        warning('on'); warning('backtrace', 'off');
        warning('The acquisition channel count does not start at 1, check the (metadata) output to see if ordered correctly');
    end

    % check if not consecutive
    if ~isempty(setdiff(min(acqChNum):max(acqChNum), acqChNum))
        warning('on'); warning('backtrace', 'off');
        warning('The acquisition channel count is not consecutive, check the (metadata) output to see if ordered correctly');
    end

    % re-order the channels in the metadata
    for iChannel = 1:length(ordAcqChNum)
        tmpStruct(iChannel) = metadata.time_series_channels(prevIndex(iChannel));
    end
    metadata.time_series_channels = tmpStruct;

    
    %
    % load signal data
    %
    
    % check if signal data should be returned
    if nargout > 1
        
        % include all channel if no specific channels were given
        if isempty(channels)
            for iChannel = 1:metadata.number_of_time_series_channels
                channels{end + 1} = metadata.time_series_channels(iChannel).name;
            end
        end
        
        % allow the request of a single channel as a string argument
        if ischar(channels), channels = {channels}; end;
        
        % check the channels input argument
        if ~iscell(channels)
            fprintf(2, 'Error: invalid input argument for ''channels'', should be a cell array containing channel names as string (e.g. {''ch1'', ''ch2'', ''ch3''})\n');
            return;
        end
        for iChannel = 1:length(channels)
            if ~ischar(channels{iChannel})
                fprintf(2, 'Error: invalid input argument for ''channels'', should be a cell array containing channel names as string (e.g. {''ch1'', ''ch2'', ''ch3''})\n');
                return;
            end 
        end
        
        % make sure all requested channels exist
        % note: if one or more channels are not found will return an error. Channel selection can be sensitive, this
        %       approach prevents unexpected consequences that could arise when - instead - returning less or empty channels
        channelsFound = ismember(lower(channels), lower({metadata.time_series_channels.name}));
        if sum(channelsFound) < length(channels)
            for iChannel = 1:length(channels)
                if channelsFound(iChannel) == 0
                    fprintf(2, 'Error: requested channel ''%s'' was not found\n', channels{iChannel});
                end
            end
            return;
        end
        
        % catch errors
        try
            
            % loop through the channels in the order they are requested
            for iChannel = 1:length(channels)
                channelIndex = find(ismember(lower({metadata.time_series_channels.name}), lower(channels{iChannel})));
                channelPath = [metadata.time_series_channels(channelIndex).path, filesep, metadata.time_series_channels(channelIndex).name, '.', metadata.time_series_channels(channelIndex).extension];
                
                %
                % read the data of each channel
                %

                % loop through the ranges
                for iRange = 1:size(ranges, 1)

                    % read the signal
                    signal = read_mef_ts_data(channelPath, password, rangeType, ranges(iRange, 1), ranges(iRange, 2))';

                    % on the first read, initialize the array
                    % note: we cannot beforehand determine the size of matrix because of potential
                    %       recording gaps in the data. Therefore, we allocate memory here.
                    if iChannel == 1 && iRange == 1
                        data = nan(length(channels), length(signal), size(ranges, 1));
                    end

                    % if signal too large for the data matrix, pad data matrix with nans
                    if length(signal) > size(data, 2)
                       data = padarray(data, [0, (length(signal) - size(data, 2))], nan, 'post');
                    end

                    % store the data 
                    data(iChannel, 1:length(signal), iRange) = signal;

                end
                
            end

        catch e
            
            % error message
            fprintf(2, '%s\nError while reading MEF3 channel data\n', e.message);
            return;
            
        end
        
    end
    
end
