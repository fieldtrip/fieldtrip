
function [streams,fileheader] = load_xdf(filename,varargin)
% Import an XDF file.
% [Streams,FileHeader] = load_xdf(Filename, Options...)
%
% This is a MATLAB importer for mult-stream XDF (Extensible Data Format) recordings. All
% information covered by the XDF 1.0 specification is imported, plus any additional meta-data
% associated with streams or with the container file itself.
%
% See https://github.com/sccn/xdf/ for more information on XDF.
%
% The function supports several further features, such as compressed XDF archives, robust
% time synchronization, support for breaks in the data, as well as some other defects.
%
% In:
%   Filename : name of the file to import (*.xdf or *.xdfz)
%
%   Options... : A list of optional name-value arguments for special use cases. The allowed names 
%                are listed in the following:
%
%                Parameters that control various processing features:
%
%                'Verbose' : Whether to print verbose diagnostics. (default: false)
%
%                'HandleClockSynchronization' : Whether to enable clock synchronization based on
%                                               ClockOffset chunks. (default: true)
%
%                'HandleJitterRemoval' : Whether to perform jitter removal for regularly sampled 
%                                        streams. (default: true)
%
%                'OnChunk' : Function that is called for each chunk of data as it
%                            is being retrieved from the file; the function is allowed to modify the
%                            data (for example, sub-sample it). The four input arguments are 1) the
%                            matrix of [#channels x #samples] values (either numeric or 2d cell
%                            array of strings), 2) the vector of unprocessed local time stamps (one
%                            per sample), 3) the info struct for the stream (same as the .info field
%                            in the final output, buth without the .effective_srate sub-field), and
%                            4) the scalar stream number (1-based integers). The three return values
%                            are 1) the (optionally modified) data, 2) the (optionally modified)
%                            time stamps, and 3) the (optionally modified) header (default: []).
%
%                'DisableVendorSpecifics' : Whether to perform certain vendor or system specific
%                                           operations. One example is the "BrainVision RDA" data
%                                           where it is necessary to pass the time stamps of a newly
%                                           introduced marker identifier channel on to the actual markers
%                                           in order to keep them perfectly in sync with the EEG data.
%                                           'DisableVendorSpecifics' takes the base names of certain
%                                           streams derived from the same source as a cell array of
%                                           strings (e.g. 'BrainVision RDA'). It is also possible to
%                                           disable vendor specifics altogether by providing the
%                                           value 'all'.
%
%                Parameters for advanced failure recovery in clock synchronization:
%
%                'HandleClockResets' : Whether the importer should check for potential resets of the
%                                      clock of a stream (e.g. computer restart during recording, or
%                                      hot-swap). Only useful if the recording system supports
%                                      recording under such circumstances. (default: true)
%
%                'ClockResetThresholdStds' : A clock reset must be accompanied by a ClockOffset
%                                            chunk being delayed by at least this many standard
%                                            deviations from the distribution. (default: 5)
%
%                'ClockResetThresholdSeconds' : A clock reset must be accompanied by a ClockOffset
%                                               chunk being delayed by at least this many seconds.
%                                               (default: 5)
%
%                'ClockResetThresholdOffsetStds' : A clock reset must be accompanied by a
%                                                  ClockOffset difference that lies at least this many
%                                                  standard deviations from the distribution. (default: 10)
%
%                'ClockResetThresholdOffsetSeconds' : A clock reset must be accompanied by a
%                                                     ClockOffset difference that is at least this
%                                                     many seconds away from the median. (default: 1)
%
%                'ClockResetMaxJitter' : Maximum tolerable jitter (in seconds of error) for clock
%                                        reset handling. (default: 5)
%
%                'CorrectStreamLags' : Apply lag correction described by timing spec sheets in the file. 
%                                      (default: true)
%
%                Parameters for jitter removal in the presence of data breaks:
%
%                'JitterBreakThresholdSeconds' : An interruption in a regularly-sampled stream of at least this
%                                                many seconds will be considered as a potential break (if also
%                                                the BreakThresholdSamples is crossed) and multiple segments
%                                                will be returned. Default: 1
%
%                'JitterBreakThresholdSamples' : An interruption in a regularly-sampled stream of at least this
%                                                many samples will be considered as a potential break (if also
%                                                the BreakThresholdSeconds is crossed) and multiple segments
%                                                will be returned. Default: 500
%                
%                Parameters for streams that can drop samples:
%
%                'FrameRateAccuracy' : How accurate the nominal frame rate is. 
%                                      Used for can_drop_samples == true frames to determine
%                                      the maximum possible number of dropped frames. Default: 0.05%
%
% Out:
%   Streams : cell array of structs, one for each stream; the structs have the following content:
%             .time_series field: contains the stream's time series [#Channels x #Samples]
%                                 this matrix is of the type declared in .info.channel_format
%             .time_stamps field: contains the time stamps for each sample (synced across streams)
%
%             .info field: contains the meta-data of the stream (all values are strings)
%               .name: name of the stream
%               .type: content-type of the stream ('EEG','Events', ...)
%               .channel_format: value format ('int8','int16','int32','int64','float32','double64','string')
%               .nominal_srate: nominal sampling rate of the stream (as declared by the device);
%                               zero for streams with irregular sampling rate
%               .effective_srate: effective (measured) sampling rate of the stream, if regular
%                                 (otherwise omitted)
%               .desc: struct with any domain-specific meta-data declared for the stream; see
%                      www.xdf.org for the declared specifications
%
%             .segments field: struct array containing segment ranges for regularly sampled
%                              time series with breaks (not present if the stream is irregular)
%               .index_range: 1st and last index of the segment within the .time_series/.time_stamps
%                             arrays
%               .t_begin: time of the 1st sample in the segment, in seconds
%               .t_end: time of the last sample in the segment, in seconds
%               .duration: duration of the segment, in seconds
%               .num_samples: number of samples in the segment
%               .effective_srate: effective (i.e. measured) sampling rate within the segment
%
%   FileHeader : struct with file header contents in the .info field
%
% Examples:
%   % load the streams contained in a given XDF file
%   streams = load_xdf('C:\Recordings\myrecording.xdf')
%
% License:
%     This file is covered by the BSD license.
%
%     Copyright (c) 2012, Christian Kothe
%     Portions Copyright (c) 2010, Wouter Falkena
%     All rights reserved.
%
%     Redistribution and use in source and binary forms, with or without
%     modification, are permitted provided that the following conditions are
%     met:
%
%         * Redistributions of source code must retain the above copyright
%           notice, this list of conditions and the following disclaimer.
%         * Redistributions in binary form must reproduce the above copyright
%           notice, this list of conditions and the following disclaimer in
%           the documentation and/or other materials provided with the distribution
%
%     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%     AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%     IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%     ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%     LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%     CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%     SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%     INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%     CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%     ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%     POSSIBILITY OF SUCH DAMAGE.
%
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-04-22
%
%                                Contains portions of xml2struct Copyright (c) 2010, Wouter Falkena,
%                                ASTI, TUDelft, 21-08-2010
%
%                                Contains frame snapping routine to handle dropped video frames by 
%                                Matthew Grivich.
%
%                                version 1.13
LIBVERSION = '1.13';
% check inputs
opts = cell2struct(varargin(2:2:end),varargin(1:2:end),2);
if ~isfield(opts,'OnChunk')
    opts.OnChunk = []; end
if ~isfield(opts,'Verbose')
    opts.Verbose = false; end
if ~isfield(opts,'HandleClockSynchronization')
    opts.HandleClockSynchronization = true; end
if ~isfield(opts,'HandleClockResets')
    opts.HandleClockResets = true; end
if ~isfield(opts,'HandleJitterRemoval')
    opts.HandleJitterRemoval = true; end
if ~isfield(opts,'JitterBreakThresholdSeconds')
    opts.JitterBreakThresholdSeconds = 1; end
if ~isfield(opts,'JitterBreakThresholdSamples')
    opts.JitterBreakThresholdSamples = 500; end
if ~isfield(opts,'ClockResetThresholdSeconds')
    opts.ClockResetThresholdSeconds = 5; end
if ~isfield(opts,'ClockResetThresholdStds')
    opts.ClockResetThresholdStds = 5; end
if ~isfield(opts,'ClockResetThresholdOffsetSeconds')
    opts.ClockResetThresholdOffsetSeconds = 1; end
if ~isfield(opts,'ClockResetThresholdOffsetStds')
    opts.ClockResetThresholdOffsetStds = 10; end
if ~isfield(opts,'WinsorThreshold')
    opts.WinsorThreshold = 0.0001; end
if ~isfield(opts,'ClockResetMaxJitter')
    opts.ClockResetMaxJitter = 5; end
if ~isfield(opts,'DisableVendorSpecifics')
    opts.DisableVendorSpecifics = {}; end
if ~isfield(opts,'CorrectStreamLags')
    opts.CorrectStreamLags = true; end
if ~isfield(opts,'FrameRateAccuracy')
    opts.FrameRateAccuracy = .05; end
if ~exist(filename,'file')
    error(['The file "' filename '" does not exist.']); end


% uncompress if necessary (note: "bonus" feature, not part of the XDF 1.0 spec)
[p,n,x] = fileparts(filename);
if strcmp(x,'.xdfz')
    % idea for this type of approach by Michael Kleder
    import com.mathworks.mlwidgets.io.InterruptibleStreamCopier
    src = java.io.FileInputStream(filename);
    flt = java.util.zip.InflaterInputStream(src);
    filename = [p filesep n '_temp_uncompressed' x];
    dst = java.io.FileOutputStream(filename);
    copier = InterruptibleStreamCopier.getInterruptibleStreamCopier;
    copier.copyStream(flt,dst);
    dst.close();
    src.close();
end

streams = {};                                   % cell array of returned streams (in the order of appearance in the file)
idmap = sparse(2^31-1,1);                       % remaps stream id's onto indices in streams
temp = struct();                                % struct array of temporary per-stream information
fileheader = struct();                          % the file header
f = fopen(filename,'r','ieee-le.l64');          % file handle
closer = onCleanup(@()close_file(f,filename));  % object that closes the file when the function exits


% there is a fast C mex file for the inner loop, but it's 
% not necessarily available for every platform
have_mex = exist('load_xdf_innerloop','file');
if ~have_mex
    if opts.Verbose
        disp(['NOTE: apparently you are missing a compiled binary version of the inner loop code.',...
            ' Attempting to download...']);
    end
    
    fname = ['load_xdf_innerloop.' mexext];
    mex_url = ['https://github.com/sccn/xdf/releases/download/v',...
        LIBVERSION, '/', fname];
    [this_path, this_name, this_ext] = fileparts(mfilename('fullpath'));
    try
        have_mex = true;
        websave(fullfile(this_path, fname), mex_url);
    catch ME
        if opts.Verbose
            disp(['Unable to download the compiled binary version for your platform.',...
                ' Using the slow MATLAB code instead.']);
        end
        have_mex = false;
        %rethrow(ME);
    end
end


% ======================
% === parse the file ===
% ======================

% read [MagicCode]
if ~strcmp(fread(f,4,'*char')','XDF:')
    error(['This is not a valid XDF file (' filename ').']); end

% for each chunk...

if opts.Verbose; fprintf('Now reading from %s ...', filename); end;
while 1
    % read [NumLengthBytes], [Length]
    len = double(read_varlen_int(f));
    if ~len
        break; end
    % read [Tag]
    switch fread(f,1,'uint16')
        case 3 % read [Samples] chunk
            try
                % read [StreamId]
                id = idmap(fread(f,1,'uint32'));
                if have_mex
                    % read the chunk data at once
                    data = fread(f,len-6,'*uint8');
                    % run the mex kernel
                    [values,timestamps] = load_xdf_innerloop(data, temp(id).chns, temp(id).readfmt, temp(id).sampling_interval, temp(id).last_timestamp);
                    temp(id).last_timestamp = timestamps(end);
                else % fallback MATLAB implementation
                    % read [NumSampleBytes], [NumSamples]
                    num = read_varlen_int(f);
                    % allocate space
                    timestamps = zeros(1,num);
                    if strcmp(temp(id).readfmt,'*string')
                        values = cell(temp(id).chns,num);
                    else
                        values = zeros(temp(id).chns,num);
                    end
                    % for each sample...
                    for s=1:num
                        % read or deduce time stamp
                        if fread(f,1,'*uint8')
                            timestamps(s) = fread(f,1,'double');
                        else
                            timestamps(s) = temp(id).last_timestamp + temp(id).sampling_interval;
                        end
                        % read the values
                        if strcmp(temp(id).readfmt,'*string')
                            for v = 1:size(values,1)
                                values{v,s} = fread(f,double(read_varlen_int(f)),'*char')'; end
                        else
                            values(:,s) = fread(f,size(values,1),temp(id).readfmt);
                        end
                        temp(id).last_timestamp = timestamps(s);
                    end
                end
                % optionally send through the OnChunk function
                if ~isempty(opts.OnChunk)
                    [values,timestamps,streams{id}] = opts.OnChunk(values,timestamps,streams{id},id); end %#ok<*AGROW>
                % append to the time series...
                temp(id).time_series{end+1} = values;
                temp(id).time_stamps{end+1} = timestamps;
            catch e
                % an error occurred (perhaps a chopped-off file): emit a warning
                % and return the file up to this point
                warning(e.identifier,e.message);
                break;
            end
        case 2 % read [StreamHeader] chunk
            % read [StreamId]
            streamid = fread(f,1,'uint32');
            id = length(streams)+1;
            idmap(streamid) = id; %#ok<SPRIX>
            % read [Content]
            header = parse_xml_struct(fread(f,len-6,'*char')');
            streams{id} = header;
            if opts.Verbose
                fprintf(['  found stream ' header.info.name '\n']); end
            % generate a few temporary fields
            temp(id).chns = str2num(header.info.channel_count); %#ok<*ST2NM>
            temp(id).srate = str2num(header.info.nominal_srate);
            temp(id).last_timestamp = 0;
            temp(id).time_series = {};
            temp(id).time_stamps = {};
            temp(id).clock_times = [];
            temp(id).clock_values = [];
            if temp(id).srate > 0
                temp(id).sampling_interval = 1/temp(id).srate;
            else
                temp(id).sampling_interval = 0;
            end
            % fread parsing format for data values
            temp(id).readfmt = ['*' header.info.channel_format];
            if strcmp(temp(id).readfmt,'*double64') && ~have_mex
                temp(id).readfmt = '*double'; end % for fread()
        case 6 % read [StreamFooter] chunk
            % read [StreamId]
            id = idmap(fread(f,1,'uint32'));
            % read [Content]
            footer = parse_xml_struct(fread(f,len-6,'*char')');
            streams{id} = hlp_superimposedata(footer,streams{id});
        case 1 % read [FileHeader] chunk
            fileheader = parse_xml_struct(fread(f,len-2,'*char')');
        case 4 % read [ClockOffset] chunk
            % read [StreamId]
            id = idmap(fread(f,1,'uint32'));
            % read [CollectionTime]
            temp(id).clock_times(end+1) = fread(f,1,'double');
            % read [OffsetValue]
            temp(id).clock_values(end+1) = fread(f,1,'double');
        otherwise
            % skip other chunk types (Boundary, ...)
            fread(f,len-2,'*uint8');
    end
end
    
% concatenate the signal across chunks
for k=1:length(temp)
    try
        temp(k).time_series = [temp(k).time_series{:}];
        temp(k).time_stamps = [temp(k).time_stamps{:}];
    catch e
        disp(['Could not concatenate time series for stream ' streams{k}.info.name '; skipping.']);
        disp(['Reason: ' e.message]);
        temp(k).time_series = [];
        temp(k).time_stamps = [];
    end
end


% ===================================================================
% === perform (fault-tolerant) clock synchronization if requested ===
% ===================================================================

if opts.HandleClockSynchronization
    if opts.Verbose
        disp('  performing clock synchronization...'); end
    for k=1:length(temp)
        if ~isempty(temp(k).time_stamps)
            try
                clock_times = temp(k).clock_times;
                clock_values = temp(k).clock_values;
                if isempty(clock_times)
                    error('No clock offset values present.'); end
            catch
                disp(['No clock offsets were available for stream "' streams{k}.info.name '"']);
                continue;
            end
 
            % detect clock resets (e.g., computer restarts during recording) if requested
            % this is only for cases where "everything goes wrong" during recording
            % note that this is a fancy feature that is not needed for normal XDF compliance
            if opts.HandleClockResets
                % first detect potential breaks in the synchronization data; this is only necessary when the
                % importer should be able to deal with recordings where the computer that served a stream
                % was restarted or hot-swapped during an ongoing recording, or the clock was reset otherwise
                time_diff = diff(clock_times);
                value_diff = abs(diff(clock_values));
                % points where a glitch in the timing of successive clock measurements happened
                time_glitch = (time_diff < 0 | (((time_diff - median(time_diff)) ./ median(abs(time_diff-median(time_diff)))) > opts.ClockResetThresholdStds & ...
                    ((time_diff - median(time_diff)) > opts.ClockResetThresholdSeconds)));
                % points where a glitch in successive clock value estimates happened
                value_glitch = (value_diff - median(value_diff)) ./ median(abs(value_diff-median(value_diff))) > opts.ClockResetThresholdOffsetStds & ...
                    (value_diff - median(value_diff)) > opts.ClockResetThresholdOffsetSeconds;
                % points where both a time glitch and a value glitch co-occur are treated as resets
                resets_at = time_glitch & value_glitch;
                % determine the [begin,end] index ranges between resets
                if any(resets_at)
                    tmp = find(resets_at)';
                    tmp = [tmp tmp+1]';
                    tmp = [1 tmp(:)' length(resets_at)];
                    ranges = num2cell(reshape(tmp,2,[])',2);
                    if opts.Verbose
                        disp(['  found ' num2str(nnz(resets_at)) ' clock resets in stream ' streams{k}.info.name '.']); end
                else
                    ranges = {[1,length(clock_times)]};
                end
            else
                % otherwise we just assume that there are no clock resets
                ranges = {[1,length(clock_times)]};
            end
            
            % calculate clock offset mappings for each data range
            mappings = {};
            for r=1:length(ranges)
                idx = ranges{r};
                if idx(1) ~= idx(2)
                    % to accomodate the Winsorizing threshold (in seconds) we rescale the data (robust_fit sets it to 1 unit)
                    mappings{r} = robust_fit([ones(idx(2)-idx(1)+1,1) clock_times(idx(1):idx(2))']/opts.WinsorThreshold, clock_values(idx(1):idx(2))'/opts.WinsorThreshold);
                else
                    mappings{r} = [clock_values(idx(1)) 0]; % just one measurement
                end
            end
            
            if length(ranges) == 1
                % apply the correction to all time stamps
                temp(k).time_stamps = temp(k).time_stamps + (mappings{1}(1) + mappings{1}(2)*temp(k).time_stamps);
            else
                % if there are data segments measured with different clocks we need to
                % determine, for any time stamp lying between two segments, to which of the segments it belongs
                clock_segments = zeros(size(temp(k).time_stamps));  % the segment index to which each stamp belongs                
                begin_of_segment = 1;                               % first index into time stamps that belongs to the current segment
                end_of_segment = NaN; %#ok<NASGU>                   % last index into time stamps that belongs to the current segment
                for r=1:length(ranges)-1
                    cur_end_time = clock_times(ranges{r}(2));       % time at which the current segment ends
                    next_begin_time = clock_times(ranges{r+1}(1));  % time at which the next segment begins
                    % get the data that is not yet processed
                    remaining_indices = begin_of_segment:length(temp(k).time_stamps);
                    if isempty(remaining_indices)
                        break; end
                    remaining_data = temp(k).time_stamps(remaining_indices);
                    if next_begin_time > cur_end_time
                        % clock jumps forward: the end of the segment is where the data time stamps
                        % lie closer to the next segment than the current in time
                        end_of_segment = remaining_indices(min(find([abs(remaining_data-cur_end_time) > abs(remaining_data-next_begin_time),true],1)-1,length(remaining_indices)));
                    else
                        % clock jumps backward: the end of the segment is where the data time stamps
                        % jump back by more than the max conceivable jitter (as any negative delta is jitter)
                        end_of_segment = remaining_indices(min(find([diff(remaining_data) < -opts.ClockResetMaxJitter,true],1),length(remaining_indices)));
                    end
                    % assign the segment of data points to the current range
                    % go to next segment
                    clock_segments(begin_of_segment:end_of_segment) = r;
                    begin_of_segment = end_of_segment+1;
                end
                % assign all remaining time stamps to the last segment
                clock_segments(begin_of_segment:end) = length(ranges);                
                % apply corrections on a per-segment basis
                for r=1:length(ranges)
                    temp(k).time_stamps(clock_segments==r) = temp(k).time_stamps(clock_segments==r) + (mappings{r}(1) + mappings{r}(2)*temp(k).time_stamps(clock_segments==r)); end
            end
        end
    end
end


% ===========================================
% === perform jitter removal if requested ===
% ===========================================
if opts.HandleJitterRemoval
    % jitter removal is a bonus feature that yields linearly increasing timestamps from data 
    % where samples had been time stamped with some jitter (e.g., due to operating system
    % delays)
    if opts.Verbose
        disp('  performing jitter removal...'); end
    
              
    for k=1:length(temp)
        
        if ~isempty(temp(k).time_stamps) && temp(k).srate
            
   
            if isfield(streams{k}.info.desc, 'synchronization') && ... 
                isfield(streams{k}.info.desc.synchronization, 'can_drop_samples') && ...
                strcmp(streams{k}.info.desc.synchronization.can_drop_samples, 'true')        
                temp(k).time_stamps = droppedFramesCorrection(temp(k).time_stamps,temp(k).srate, opts.FrameRateAccuracy);            
           else
            
                % identify breaks in the data
                diffs = diff(temp(k).time_stamps);
                breaks_at = diffs > max(opts.JitterBreakThresholdSeconds,opts.JitterBreakThresholdSamples*temp(k).sampling_interval);
                if any(breaks_at)
                    % turn the break mask into a cell array of [begin,end] index ranges
                    tmp = find(breaks_at)';
                    tmp = [tmp tmp+1]';
                    tmp = [1 tmp(:)' length(breaks_at)];
                    ranges = num2cell(reshape(tmp,2,[])',2);
                    if opts.Verbose
                        disp(['  found ' num2str(nnz(breaks_at)) ' data breaks in stream ' streams{k}.info.name '.']); end
                else
                    ranges = {[1,length(temp(k).time_stamps)]};
                end

                % process each segment separately
                segments = repmat(struct(),1,length(ranges));
                for r=1:length(ranges)
                    range = ranges{r};
                    segments(r).num_samples = range(2)-range(1)+1;
                    segments(r).index_range = range;
                    if segments(r).num_samples > 0
                        indices = segments(r).index_range(1):segments(r).index_range(2);
                        % regress out the jitter
                        mapping = temp(k).time_stamps(indices) / [ones(1,length(indices)); indices];
                        temp(k).time_stamps(indices) = mapping(1) + mapping(2) * indices;
                    end
                    % calculate some other meta-data about the segments
                    segments(r).t_begin = temp(k).time_stamps(range(1));
                    segments(r).t_end = temp(k).time_stamps(range(2));
                    segments(r).duration = segments(r).t_end - segments(r).t_begin;
                    segments(r).effective_srate = (segments(r).num_samples - 1) / segments(r).duration;
                end

                % calculate the weighted mean sampling rate over all segments
                temp(k).effective_rate = sum(bsxfun(@times,[segments.effective_srate],[segments.num_samples]/sum([segments.num_samples])));            

                % transfer the information into the output structs
                streams{k}.info.effective_srate = temp(k).effective_rate;
                streams{k}.segments = segments;
            end
        end
    end
else
    % calculate effective sampling rate
    for k=1:length(temp)
        temp(k).effective_srate = (length(temp(k).time_stamps) - 1) / (temp(k).time_stamps(end) - temp(k).time_stamps(1)); end
end

% copy the information into the output
for k=1:length(temp)
    offset_mean = 0;
    if opts.CorrectStreamLags && ...
        isfield(streams{k}.info.desc, 'synchronization') && ... 
        isfield(streams{k}.info.desc.synchronization, 'offset_mean')
                offset_mean = str2num(streams{k}.info.desc.synchronization.offset_mean);        
    end
    streams{k}.time_series = temp(k).time_series;
    streams{k}.time_stamps = temp(k).time_stamps - offset_mean;
end


% =========================================
% === peform vendor specific operations ===
% =========================================

if ~any(strcmp('all',opts.DisableVendorSpecifics))

    % BrainVision RDA
    targetName = 'BrainVision RDA';
    if ~any(strcmp(opts.DisableVendorSpecifics,targetName))

        % find a target EEG stream...
        for k=1:length(streams)

            if strcmp(streams{k}.info.name,targetName) % Is a BrainVision RDA stream?
                mkChan = [];
                for iChan = 1:length( streams{ k }.info.desc.channels.channel ) % Find marker index channel (any channel, not necessary last)
                    if strcmp( streams{ k }.info.desc.channels.channel{ iChan }.label, 'MkIdx' ) && strcmp( streams{ k }.info.desc.channels.channel{ iChan }.type, 'Marker' ) && strcmp( streams{ k }.info.desc.channels.channel{ iChan }.unit, 'Counts (decimal)' )
                        mkChan = iChan;
                        break % Only one marker channel expected
                    end
                end
                if ~isempty( mkChan ) % Has a marker index channel?
                    for m = 1:length( streams ) % find a corresponding indexed marker stream...
                        if strcmp( streams{ m }.info.name, [ targetName ' Markers' ] ) && strcmp( streams{ m }.info.hostname, streams{ k }.info.hostname ) && strncmp( streams{ m }.info.source_id, streams{ k }.info.source_id, length( streams{ k }.info.source_id ) )
                            if opts.Verbose
                                disp( [ '  performing ', targetName, ' specific tasks for stream ', num2str( k ), '...' ] );
                            end
                            streams = ProcessBVRDAindexedMarkers( streams, k, m, mkChan );
                        end
                    end
					
					% Remove marker index channel
					streams{ k }.time_series( mkChan, : ) = []; 
					streams{ k }.info.desc.channels.channel( mkChan ) = [];
					
					% Decrement channel count by 1
					streams{ k }.info.channel_count = num2str( str2num( streams{ k }.info.channel_count ) - 1 );
					
                end
				
            end

        end

    end
    
end


end


% ========================
% === helper functions ===
% ========================

% read a variable-length integer
function num = read_varlen_int(f)
try
    switch fread(f,1,'*uint8')
        case 1
            num = fread(f,1,'*uint8');
        case 4
            num = fread(f,1,'*uint32');
        case 8
            num = fread(f,1,'*uint64');
        otherwise
            error('Invalid variable-length integer encountered.');
    end
catch %#ok<*CTCH>
    num = 0;
end
end


% close the file and delete temporary data
function close_file(f,filename)
fclose(f);
if strfind(filename,'_temp_uncompressed.xdf')
    delete(filename); end
end


% parse a simplified (attribute-free) subset of XML into a MATLAB struct
function result = parse_xml_struct(str)
import org.xml.sax.InputSource
import javax.xml.parsers.*
import java.io.*
tmp = InputSource();
tmp.setCharacterStream(StringReader(str));
result = parseChildNodes(xmlread(tmp));

% this is part of xml2struct (slightly simplified)
    function [children,ptext] = parseChildNodes(theNode)
        % Recurse over node children.
        children = struct;
        ptext = [];
        if theNode.hasChildNodes
            childNodes = theNode.getChildNodes;
            numChildNodes = childNodes.getLength;
            for count = 1:numChildNodes
                theChild = childNodes.item(count-1);
                [text,name,childs] = getNodeData(theChild);
                if (~strcmp(name,'#text') && ~strcmp(name,'#comment'))
                    if (isfield(children,name))
                        if (~iscell(children.(name)))
                            children.(name) = {children.(name)}; end
                        index = length(children.(name))+1;
                        children.(name){index} = childs;
                        if(~isempty(text))
                            children.(name){index} = text; end
                    else
                        children.(name) = childs;
                        if(~isempty(text))
                            children.(name) = text; end
                    end
                elseif (strcmp(name,'#text'))
                    if (~isempty(regexprep(text,'[\s]*','')))
                        if (isempty(ptext))
                            ptext = text;
                        else
                            ptext = [ptext text];
                        end
                    end
                end
            end
        end
    end

% this is part of xml2struct (slightly simplified)
    function [text,name,childs] = getNodeData(theNode)
        % Create structure of node info.
        name = char(theNode.getNodeName);
        if ~isvarname(name)
            name = regexprep(name,'[-]','_dash_');
            name = regexprep(name,'[:]','_colon_');
            name = regexprep(name,'[.]','_dot_');
        end
        [childs,text] = parseChildNodes(theNode);
        if (isempty(fieldnames(childs)))
            try
                text = char(theNode.getData);
            catch
            end
        end
    end
end


function x = robust_fit(A,y,rho,iters)
% Perform a robust linear regression using the Huber loss function.
% x = robust_fit(A,y,rho,iters)
%
% Input:
%   A : design matrix
%   y : target variable
%   rho : augmented Lagrangian variable (default: 1)
%   iters : number of iterations to perform (default: 1000)
%
% Output:
%   x : solution for x
%
% Notes:
%   solves the following problem via ADMM for x:
%     minimize 1/2*sum(huber(A*x - y))
%
% Based on the ADMM Matlab codes also found at:
%   http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-03-04

if ~exist('rho','var')
    rho = 1; end
if ~exist('iters','var')
    iters = 1000; end
Aty = A'*y;
L = sparse(chol(A'*A,'lower')); U = L';
z = zeros(size(y)); u = z;
for k = 1:iters
    x = U \ (L \ (Aty + A'*(z - u)));
    d = A*x - y + u;
    z = rho/(1+rho)*d + 1/(1+rho)*max(0,(1-(1+1/rho)./abs(d))).*d;
    u = d - z;
end
end


function res = hlp_superimposedata(varargin)
% Merge multiple partially populated data structures into one fully populated one.
% Result = hlp_superimposedata(Data1, Data2, Data3, ...)
%
% The function is applicable when you have cell arrays or structs/struct arrays with non-overlapping
% patterns of non-empty entries, where all entries should be merged into a single data structure
% which retains their original positions. If entries exist in multiple data structures at the same
% location, entries of later items will be ignored (i.e. earlier data structures take precedence).
%
% In:
%   DataK : a data structure that should be super-imposed with the others to form a single data
%           structure
%
% Out:
%   Result : the resulting data structure
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-08-19

% first, compactify the data by removing the empty items
compact = varargin(~cellfun('isempty',varargin));
% start with the last data structure, then merge the remaining data structures into it (in reverse
% order as this avoids having to grow arrays incrementally in typical cases)
res = compact{end};
for k=length(compact)-1:-1:1
    res = merge(res,compact{k}); end
end

function A = merge(A,B)
% merge data structures A and B
if iscell(A) && iscell(B)
    % make sure that both have the same number of dimensions
    if ndims(A) > ndims(B)
        B = grow_cell(B,size(A));
    elseif ndims(A) < ndims(B)
        A = grow_cell(A,size(B));
    end
    % make sure that both have the same size
    if all(size(B)==size(A))
        % we're fine
    elseif all(size(B)>=size(A))
        % A is a minor of B: grow A
        A = grow_cell(A,size(B));
    elseif all(size(A)>=size(B))
        % B is a minor of A: grow B
        B = grow_cell(B,size(A));
    else
        % A and B have mixed sizes... grow both as necessary
        M = max(size(A),size(B));
        A = grow_cell(A,M);
        B = grow_cell(B,M);
    end
    % find all non-empty elements in B
    idx = find(~cellfun(@(x)isequal(x,[]),B));
    if ~isempty(idx)
        % check if any of these is occupied in A
        clean = cellfun('isempty',A(idx));
        if ~all(clean)
            % merge all conflicting items recursively
            conflicts = idx(~clean);
            for k=conflicts(:)'
                A{k} = merge(A{k},B{k}); end
            % and transfer the rest
            if any(clean)
                A(idx(clean)) = B(idx(clean)); end
        else
            % transfer all to A
            A(idx) = B(idx);
        end
    end
elseif isstruct(A) && isstruct(B)
    % first make sure that both have the same fields
    fnA = fieldnames(A);
    fnB = fieldnames(B);
    if isequal(fnA,fnB)
        % we're fine
    elseif isequal(sort(fnA),sort(fnB))
        % order doesn't match -- impose A's order on B
        B = orderfields(B,fnA);
    elseif isempty(setdiff(fnA,fnB))
        % B has a superset of A's fields: add the remaining fields to A, and order them according to B
        remaining = setdiff(fnB,fnA);
        for fn = remaining'
            A(1).(fn{1}) = []; end
        A = orderfields(A,fnB);
    elseif isempty(setdiff(fnB,fnA))
        % A has a superset of B's fields: add the remaining fields to B, and order them according to A
        remaining = setdiff(fnA,fnB);
        for fn = remaining'
            B(1).(fn{1}) = []; end
        B = orderfields(B,fnA);
    else
        % A and B have incommensurable fields; add B's fields to A's fields, add A's fields to B's
        % and order according to A's fields
        remainingB = setdiff(fnB,fnA);
        for fn = remainingB'
            A(1).(fn{1}) = []; end
        remainingA = setdiff(fnA,fnB);
        for fn = remainingA'
            B(1).(fn{1}) = []; end
        B = orderfields(B,A);
    end
    % that being established, convert them to cell arrays, merge their cell arrays, and convert back to structs
    merged = merge(struct2cell(A),struct2cell(B));
    A = cell2struct(merged,fieldnames(A),1);
elseif isstruct(A) && ~isstruct(B)
    if ~isempty(B)
        error('One of the sub-items is a struct, and the other one is of a non-struct type.');
    else
        % we retain A
    end
elseif isstruct(B) && ~isstruct(A)
    if ~isempty(A)
        error('One of the sub-items is a struct, and the other one is of a non-struct type.');
    else
        % we retain B
        A = B;
    end
elseif iscell(A) && ~iscell(B)
    if ~isempty(B)
        error('One of the sub-items is a cell array, and the other one is of a non-cell type.');
    else
        % we retain A
    end
elseif iscell(B) && ~iscell(A)
    if ~isempty(A)
        error('One of the sub-items is a cell array, and the other one is of a non-cell type.');
    else
        % we retain B
        A = B;
    end
elseif isempty(A) && ~isempty(B)
    % we retain B
    A = B;
elseif isempty(B) && ~isempty(A)
    % we retain A
elseif ~isequal(A,B)
    % we retain A and warn about dropping B
    disp('Two non-empty (and non-identical) sub-elements occupied the same index; one was dropped. This warning will only be displayed once.');
end
end


function C = grow_cell(C,idx)
% grow a cell array to accomodate a particular index
% (assuming that this index is not contained in the cell array yet)
tmp = sprintf('%i,',idx);
eval(['C{' tmp(1:end-1) '} = [];']);
end

%function written by Matthew Grivich to handle case where frame rate is
%consistent but with dropped events, like a video camera.
function frameTimesModeled = droppedFramesCorrection(frameTimes, nominalFrameRate, frameRateAccuracy)

  %  figure
  %  plot(frameTimes(1:end-1), frameTimes(2:end)-frameTimes(1:end-1));
  %  xlabel('Frame Time (s)');
  %  ylabel('Frame Interval (s)');
  
    

    nFrames = length(frameTimes); %nFrames shown, does not included dropped.
    
    %calculates the maximum conceivable number of dropped frames, given the
    %nominal frame rate and the expected frame rate accuracy.
    maxDropped = round((frameTimes(end)-frameTimes(1))*nominalFrameRate*(1+frameRateAccuracy) - nFrames);
    stds = zeros(0, 1); %initialize array of zero length
    for iteration = 1:maxDropped
        interval = (frameTimes(end)-frameTimes(1))/(nFrames-1+iteration-1);
   %     fprintf('frequency: %1.20f\n',1/interval);
        frameNumbers = zeros(1,length(frameTimes));

        for i=1:length(frameNumbers)

            frameNumbers(i) =  round((frameTimes(i)-frameTimes(1))/interval)+1;

        end
        %remove zero frame intervals
        for i=length(frameNumbers)-1:-1:1
            if(frameNumbers(i) == frameNumbers(i+1))
                frameNumbers(i) = frameNumbers(i) -1;
            end
        end

    
        pf = polyfit(frameNumbers, frameTimes,1);
        %text(0,.0004, sprintf('interval: %1.20f\n', interval));
    %    fprintf('interval after fit: %1.20f\n',pf(1));

        frameTimesModeled = polyval(pf,frameNumbers);
  %     plot(stds)
      %  plot(frameTimesModeled, frameTimesModeled - frameTimes);
  %      text(0,.02, sprintf('dropped: %d\n', iteration-1));
  %      xlabel('Frame Time (s)')
  %      ylabel('Frame Time Modeled - Frame Time (s)');
  %     pause(0.001);

        stds(iteration) = std(frameTimesModeled-frameTimes);
        %if (mean(stds) - min(stds))/std(stds) > 10
        %    break;
        %end


    end

 %   figure
 %   plot(stds);
    dropped = find(stds==min(stds)) - 1;

    interval = (frameTimes(end)-frameTimes(1))/(nFrames-1+dropped);
   % fprintf('interval: %1.20f\n',interval);
    frameNumbers = zeros(1,length(frameTimes));
    
    %Find closest frame.
    for i=1:length(frameNumbers)
        frameNumbers(i) =  round((frameTimes(i)-frameTimes(1))/interval)+1;
    end
    %remove zero frame intervals
    for i=length(frameNumbers)-1:-1:1
        if(frameNumbers(i) >= frameNumbers(i+1))
            frameNumbers(i) = frameNumbers(i+1) -1;
        end
    end


    pf = polyfit(frameNumbers, frameTimes,1);
    %text(0,.0004, sprintf('interval: %1.20f\n', interval));
%    fprintf('interval after fit: %1.20f\n',pf(1));

    frameTimesModeled = polyval(pf,frameNumbers);
%figure     
%        plot(frameTimesModeled, smooth(frameTimesModeled - frameTimes,21));
%        xlabel('Frame Time (s)')
%        ylabel('Frame Time Modeled - Frame Time (s)');
 %      pause(0.5);
      
    
%    fprintf('Interval: %1.20f\n',pf(1));
%    fprintf('Dropped Frames: %d\n', dropped);
    

   % figure
   % plot(frameTimes(1:end-1), pf(1)*(frameNumbers(2:end)-frameNumbers(1:end-1)));
   % xlabel('Frame Time (s)');
   % ylabel('Frame Interval (s)');
 

end


function streams = ProcessBVRDAindexedMarkers( streams, dataStream, mkStream, mkChan )

clearMarkers = [];

for iMrk = 1:length( streams{ mkStream }.time_series )

    % Decode marker
    MrkInfo = regexp( streams{ mkStream }.time_series{ iMrk }, 'mk(?<idx>\d+)=(?<str>.*)', 'names' );

    if ~isempty( MrkInfo.idx )

        % Find corresponding sample in marker index channel
        lat = find( streams{ dataStream }.time_series( mkChan, : ) == str2double( MrkInfo.idx ) );
        offset = streams{ mkStream }.time_stamps( iMrk ) - streams{ dataStream }.time_stamps( lat );

        % Is the index unique (overflow)?
        if length( lat ) > 1
            [ minOffset, minOffsetIdx ] = min( abs( offset ) ); %#ok<ASGLU>
            lat = lat( minOffsetIdx );
            offset = offset( minOffsetIdx );
        end

        % Sanity check
        if offset > 10
            warning( 'Time stamp difference between indexed marker %s (%.3f s) and corresponding sample in marker channel (%.3f s) exceeding threshold.\n', MrkInfo.idx, streams{ mkStream }.time_stamps( iMrk ), streams{ dataStream }.time_stamps( lat ) )
        end

        % Copy time stamp and rewrite marker
        if ~isempty( lat )
            streams{ mkStream }.time_stamps( iMrk ) = streams{ dataStream }.time_stamps( lat );
            streams{ mkStream }.time_series{ iMrk } = MrkInfo.str;
        else
            warning( 'No corresponding sample found in marker channel for indexed marker %s. Removing...', MrkInfo.idx )
            clearMarkers = [ clearMarkers iMrk ]; %#ok<AGROW>
        end

    end

end

% Remove markers without corresponding marker channel sample
streams{ mkStream }.time_stamps( clearMarkers ) = [];
streams{ mkStream }.time_series( clearMarkers ) = [];


end
