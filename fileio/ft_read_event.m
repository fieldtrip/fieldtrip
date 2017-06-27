function [event] = ft_read_event(filename, varargin)

% FT_READ_EVENT reads all events from an EEG/MEG dataset and returns
% them in a well defined structure. It is a wrapper around different
% EEG/MEG file importers, directly supported formats are CTF, Neuromag,
% EEP, BrainVision, Neuroscan, Neuralynx and Nervus/Nicolet.
%
% Use as
%   [event] = ft_read_event(filename, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'dataformat'     string
%   'headerformat'   string
%   'eventformat'    string
%   'header'         header structure, see FT_READ_HEADER
%   'detectflank'    string, can be 'bit', 'up', 'down', 'both', 'peak', 'trough' or 'auto' (default is system specific)
%   'chanindx'       list with channel indices in case of different sampling frequencies (only for EDF)
%   'trigshift'      integer, number of samples to shift from flank to detect trigger value (default = 0)
%   'trigindx'       list with channel numbers for the trigger detection, only for Yokogawa (default is automatic)
%   'triglabel'      list of channel labels for the trigger detection (default is all ADC* channels for Artinis oxy3-files)
%   'threshold'      threshold for analog trigger channels (default is system specific)
%   'blocking'       wait for the selected number of events (default = 'no')
%   'timeout'        amount of time in seconds to wait when blocking (default = 5)
%   'tolerance'      tolerance in samples when merging analogue trigger channels, only for Neuromag (default = 1, meaning
%                    that an offset of one sample in both directions is compensated for)
%
% Furthermore, you can specify optional arguments as key-value pairs
% for filtering the events, e.g. to select only events of a specific
% type, of a specific value, or events between a specific begin and
% end sample. This event filtering is especially usefull for real-time
% processing. See FT_FILTER_EVENT for more details.
%
% Some data formats have trigger channels that are sampled continuously with
% the same rate as the electrophysiological data. The default is to detect
% only the up-going TTL flanks. The trigger events will correspond with the
% first sample where the TTL value is up. This behaviour can be changed
% using the 'detectflank' option, which also allows for detecting the
% down-going flank or both. In case of detecting the down-going flank, the
% sample number of the event will correspond with the first sample at which
% the TTF went down, and the value will correspond to the TTL value just
% prior to going down.
%
% This function returns an event structure with the following fields
%   event.type      = string
%   event.sample    = expressed in samples, the first sample of a recording is 1
%   event.value     = number or string
%   event.offset    = expressed in samples
%   event.duration  = expressed in samples
%   event.timestamp = expressed in timestamp units, which vary over systems (optional)
%
% The event type and sample fields are always defined, other fields can be empty,
% depending on the type of event file. Events are sorted by the sample on
% which they occur. After reading the event structure, you can use the
% following tricks to extract information about those events in which you
% are interested.
%
% Determine the different event types
%   unique({event.type})
%
% Get the index of all trial events
%   find(strcmp('trial', {event.type}))
%
% Make a vector with all triggers that occurred on the backpanel
%   [event(find(strcmp('backpanel trigger', {event.type}))).value]
%
% Find the events that occurred in trial 26
%   t=26; samples_trials = [event(find(strcmp('trial', {event.type}))).sample];
%   find([event.sample]>samples_trials(t) & [event.sample]<samples_trials(t+1))
%
% The list of supported file formats can be found in FT_READ_HEADER.
%
% To use an external reading function, use key-value pair: 'eventformat', FUNCTION_NAME.
% (Function needs to be on the path, and take as input: filename)
%
% See also FT_READ_HEADER, FT_READ_DATA, FT_WRITE_EVENT, FT_FILTER_EVENT

% Copyright (C) 2004-2016 Robert Oostenveld (Nervus by Jan Brogger)
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

global event_queue        % for fcdc_global
persistent sock           % for fcdc_tcp
persistent db_blob        % for fcdc_mysql

if isempty(db_blob)
  db_blob = 0;
end

if iscell(filename)
  ft_warning('concatenating events from %d files', numel(filename));
  % use recursion to read events from multiple files
  
  hdr = ft_getopt(varargin, 'header');
  if isempty(hdr) || ~isfield(hdr, 'orig') || ~iscell(hdr.orig)
    for i=1:numel(filename)
      % read the individual file headers
      hdr{i}  = ft_read_header(filename{i}, varargin{:});
    end
  else
    % use the individual file headers that were read previously
    hdr = hdr.orig;
  end
  nsmp = nan(size(filename));
  for i=1:numel(filename)
    nsmp(i) = hdr{i}.nSamples*hdr{i}.nTrials;
  end
  offset = [0 cumsum(nsmp(1:end-1))];
  
  event = cell(size(filename));
  for i=1:numel(filename)
    varargin = ft_setopt(varargin, 'header', hdr{i});
    event{i} = ft_read_event(filename{i}, varargin{:});
    for j=1:numel(event{i})
      % add the offset due to the previous files
      event{i}(j).sample = event{i}(j).sample + offset(i);
    end
  end
  % return the concatenated events
  event = appendevent(event{:});
  return
end

% optionally get the data from the URL and make a temporary local copy
filename = fetch_url(filename);

% get the options
hdr              = ft_getopt(varargin, 'header');
detectflank      = ft_getopt(varargin, 'detectflank', 'up', true);   % note that emptymeaningful=true
trigshift        = ft_getopt(varargin, 'trigshift');                 % default is assigned in subfunction
trigindx         = ft_getopt(varargin, 'trigindx');                  % this allows to override the automatic trigger channel detection (e.g., useful for Yokogawa)
headerformat     = ft_getopt(varargin, 'headerformat');
dataformat       = ft_getopt(varargin, 'dataformat');
threshold        = ft_getopt(varargin, 'threshold');                 % this is used for analog channels
tolerance        = ft_getopt(varargin, 'tolerance', 1);
checkmaxfilter   = ft_getopt(varargin, 'checkmaxfilter');            % will be passed to ft_read_header
eventformat      = ft_getopt(varargin, 'eventformat');
chanindx         = ft_getopt(varargin, 'chanindx');                  % used for EDF files with variable sampling rate

if isempty(eventformat)
  % only do the autodetection if the format was not specified
  eventformat = ft_filetype(filename);
end

if iscell(eventformat)
  % this happens for datasets specified as cell array for concatenation
  eventformat = eventformat{1};
end

% this allows to read only events in a certain range, supported for selected data formats only
flt_type         = ft_getopt(varargin, 'type');
flt_value        = ft_getopt(varargin, 'value');
flt_minsample    = ft_getopt(varargin, 'minsample');
flt_maxsample    = ft_getopt(varargin, 'maxsample');
flt_mintimestamp = ft_getopt(varargin, 'mintimestamp');
flt_maxtimestamp = ft_getopt(varargin, 'maxtimestamp');
flt_minnumber    = ft_getopt(varargin, 'minnumber');
flt_maxnumber    = ft_getopt(varargin, 'maxnumber');

% this allows blocking reads to avoid having to poll many times for online processing
blocking         = ft_getopt(varargin, 'blocking', false); % true or false
timeout          = ft_getopt(varargin, 'timeout', 5); % seconds

% convert from 'yes'/'no' into boolean
blocking = istrue(blocking);

if any(strcmp(eventformat, {'brainvision_eeg', 'brainvision_dat'}))
  [p, f] = fileparts(filename);
  filename = fullfile(p, [f '.vhdr']);
  eventformat = 'brainvision_vhdr';
end

if strcmp(eventformat, 'brainvision_vhdr')
  % read the headerfile belonging to the dataset and try to determine the corresponding markerfile
  eventformat = 'brainvision_vmrk';
  hdr = read_brainvision_vhdr(filename);
  % replace the filename with the filename of the markerfile
  if ~isfield(hdr, 'MarkerFile') || isempty(hdr.MarkerFile)
    filename = [];
  else
    [p, f] = fileparts(filename);
    filename = fullfile(p, hdr.MarkerFile);
  end
end

% start with an empty event structure
event = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the events with the low-level reading function
% please maintain this list in alphabetical order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch eventformat
  
  case {'4d' '4d_pdf', '4d_m4d', '4d_xyz'}
    if isempty(hdr)
      hdr = ft_read_header(filename, 'headerformat', eventformat);
    end
    
    % read the trigger channel and do flank detection
    trgindx = match_str(hdr.label, 'TRIGGER');
    if isfield(hdr, 'orig') && isfield(hdr.orig, 'config_data') && (strcmp(hdr.orig.config_data.site_name, 'Glasgow') || strcmp(hdr.orig.config_data.site_name, 'Marseille')),
      trigger = read_trigger(filename, 'header', hdr, 'dataformat', '4d', 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', trgindx, 'detectflank', detectflank, 'trigshift', trigshift, 'fix4d8192', true);
    else
      trigger = read_trigger(filename, 'header', hdr, 'dataformat', '4d', 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', trgindx, 'detectflank', detectflank, 'trigshift', trigshift, 'fix4d8192', false);
    end
    event   = appendevent(event, trigger);
    
    respindx = match_str(hdr.label, 'RESPONSE');
    if ~isempty(respindx)
      response = read_trigger(filename, 'header', hdr, 'dataformat', '4d', 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', respindx, 'detectflank', detectflank, 'trigshift', trigshift);
      event    = appendevent(event, response);
    end
    
  case 'bci2000_dat'
    % this requires the load_bcidat mex file to be present on the path
    ft_hastoolbox('BCI2000', 1);
    
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    
    if isfield(hdr.orig, 'signal') && isfield(hdr.orig, 'states')
      % assume that the complete data is stored in the header, this speeds up subsequent read operations
      signal        = hdr.orig.signal;
      states        = hdr.orig.states;
      parameters    = hdr.orig.parameters;
      total_samples = hdr.orig.total_samples;
    else
      [signal, states, parameters, total_samples] = load_bcidat(filename);
    end
    
    list = fieldnames(states);
    % loop over all states and detect the flanks, the following code was taken from read_trigger
    for i=1:length(list)
      channel   = list{i};
      trig      = double(getfield(states, channel));
      pad       = trig(1);
      trigshift = 0;
      begsample = 1;
      
      switch detectflank
        case 'up'
          % convert the trigger into an event with a value at a specific sample
          for j=find(diff([pad trig(:)'])>0)
            event(end+1).type   = channel;
            event(end  ).sample = j + begsample - 1;      % assign the sample at which the trigger has gone down
            event(end  ).value  = trig(j+trigshift);      % assign the trigger value just _after_ going up
          end
        case 'down'
          % convert the trigger into an event with a value at a specific sample
          for j=find(diff([pad trig(:)'])<0)
            event(end+1).type   = channel;
            event(end  ).sample = j + begsample - 1;      % assign the sample at which the trigger has gone down
            event(end  ).value  = trig(j-1-trigshift);    % assign the trigger value just _before_ going down
          end
        case 'both'
          % convert the trigger into an event with a value at a specific sample
          for j=find(diff([pad trig(:)'])>0)
            event(end+1).type   = [channel '_up'];        % distinguish between up and down flank
            event(end  ).sample = j + begsample - 1;      % assign the sample at which the trigger has gone down
            event(end  ).value  = trig(j+trigshift);      % assign the trigger value just _after_ going up
          end
          % convert the trigger into an event with a value at a specific sample
          for j=find(diff([pad trig(:)'])<0)
            event(end+1).type   = [channel '_down'];      % distinguish between up and down flank
            event(end  ).sample = j + begsample - 1;      % assign the sample at which the trigger has gone down
            event(end  ).value  = trig(j-1-trigshift);    % assign the trigger value just _before_ going down
          end
        otherwise
          error('incorrect specification of ''detectflank''');
      end
    end
    
  case 'besa_besa'
    % read the header
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    
    if ~isempty(detectflank) % parse the trigger channel (indicated by chanindx) for events
      event = read_trigger(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', chanindx, 'detectflank', detectflank, 'trigshift', trigshift, 'threshold', threshold);
    elseif issubfield(hdr, 'orig.events') && ~isempty(hdr.orig.events.offsets) % FIXME: add support for reading in events from the datafile
    else
      event = [];
    end
    
  case {'besa_avr', 'besa_swf'}
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    event(end+1).type     = 'average';
    event(end  ).sample   = 1;
    event(end  ).duration = hdr.nSamples;
    event(end  ).offset   = -hdr.nSamplesPre;
    event(end  ).value    = [];
    
  case {'biosemi_bdf', 'bham_bdf'}
    % read the header, required to determine the stimulus channels and trial specification
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    
    % specify the range to search for triggers, default is the complete file
    if ~isempty(flt_minsample)
      begsample = flt_minsample;
    else
      begsample = 1;
    end
    if ~isempty(flt_maxsample)
      endsample = flt_maxsample;
    else
      endsample = hdr.nSamples*hdr.nTrials;
    end
    
    if ~strcmp(detectflank, 'up')
      if strcmp(detectflank, 'both')
        warning('only up-going flanks are supported for Biosemi');
        detectflank = 'up';
      else
        error('only up-going flanks are supported for Biosemi');
        % FIXME the next section on trigger detection should be merged with the
        % READ_CTF_TRIGGER (which also does masking with bit-patterns) into the
        % READ_TRIGGER function
      end
    end
    
    % find the STATUS channel and read the values from it
    schan = find(strcmpi(hdr.label,'STATUS'));
    sdata = ft_read_data(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', schan);
    
    if ft_platform_supports('int32_logical_operations')
      % convert to 32-bit integer representation and only preserve the lowest 24 bits
      sdata = bitand(int32(sdata), 2^24-1);
    else
      % find indices of negative numbers
      bit24i = find(sdata < 0);
      % make number positive and preserve bits 0-22
      sdata(bit24i) = bitcmp(abs(sdata(bit24i))-1,24);
      % re-insert the sign bit on its original location, i.e. bit24
      sdata(bit24i) = sdata(bit24i)+(2^(24-1));
      % typecast the data to ensure that the status channel is represented in 32 bits
      sdata = uint32(sdata);
    end
    
    byte1 = 2^8  - 1;
    byte2 = 2^16 - 1 - byte1;
    byte3 = 2^24 - 1 - byte1 - byte2;
    
    % get the respective status and trigger bits
    trigger = bitand(sdata, bitor(byte1, byte2)); % this is contained in the lower two bytes
    epoch   = int8(bitget(sdata, 16+1));
    cmrange = int8(bitget(sdata, 20+1));
    battery = int8(bitget(sdata, 22+1));
    
    % determine when the respective status bits go up or down
    flank_trigger = diff([0 trigger]);
    flank_epoch   = diff([0 epoch ]);
    flank_cmrange = diff([0 cmrange]);
    flank_battery = diff([0 battery]);
    
    for i=find(flank_trigger>0)
      event(end+1).type   = 'STATUS';
      event(end  ).sample = i + begsample - 1;
      event(end  ).value  = double(trigger(i));
    end
    
    for i=find(flank_epoch==1)
      event(end+1).type   = 'Epoch';
      event(end  ).sample = i;
    end
    
    for i=find(flank_cmrange==1)
      event(end+1).type   = 'CM_in_range';
      event(end  ).sample = i;
    end
    
    for i=find(flank_cmrange==-1)
      event(end+1).type   = 'CM_out_of_range';
      event(end  ).sample = i;
    end
    
    for i=find(flank_battery==1)
      event(end+1).type   = 'Battery_low';
      event(end  ).sample = i;
    end
    
    for i=find(flank_battery==-1)
      event(end+1).type   = 'Battery_ok';
      event(end  ).sample = i;
    end
    
  case {'biosig', 'gdf'}
    % FIXME it would be nice to figure out how sopen/sread return events
    % for all possible fileformats that can be processed with biosig
    %
    % This section of code is opaque with respect to the gdf file being a
    % single file or the first out of a sequence with postfix _1, _2, ...
    % because it uses private/read_trigger which again uses ft_read_data
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    % the following applies to Biosemi data that is stored in the gdf format
    statusindx = find(strcmp(hdr.label, 'STATUS'));
    if length(statusindx)==1
      % represent the rising flanks in the STATUS channel as events
      event = read_trigger(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', statusindx, 'detectflank', 'up', 'trigshift', trigshift, 'fixbiosemi', true);
    else
      warning('BIOSIG does not have a consistent event representation, skipping events')
      event = [];
    end
    
  case 'AnyWave'
    event = read_ah5_markers(hdr, filename);
  case 'brainvision_vmrk'
    fid=fopen(filename,'rt');
    if fid==-1,
      error('cannot open BrainVision marker file')
    end
    line = [];
    while ischar(line) || isempty(line)
      line = fgetl(fid);
      if ~isempty(line) && ~(isnumeric(line) && line==-1)
        if strncmpi(line, 'Mk', 2)
          % this line contains a marker
          tok = tokenize(line, '=', 0);    % do not squeeze repetitions of the separator
          if length(tok)~=2
            warning('skipping unexpected formatted line in BrainVision marker file');
          else
            % the line looks like "MkXXX=YYY", which is ok
            % the interesting part now is in the YYY, i.e. the second token
            tok = tokenize(tok{2}, ',', 0);    % do not squeeze repetitions of the separator
            if isempty(tok{1})
              tok{1}  = [];
            end
            if isempty(tok{2})
              tok{2}  = [];
            end
            event(end+1).type     = tok{1};
            event(end  ).value    = tok{2};
            event(end  ).sample   = str2num(tok{3});
            event(end  ).duration = str2num(tok{4});
          end
        end
      end
    end
    fclose(fid);
    
  case 'bucn_nirs'
    event = read_bucn_nirsevent(filename);
    
  case 'ced_son'
    % check that the required low-level toolbox is available
    ft_hastoolbox('neuroshare', 1);
    orig = read_ced_son(filename,'readevents','yes');
    event = struct('type',     {orig.events.type},...
      'sample',   {orig.events.sample},...
      'value',    {orig.events.value},...
      'offset',   {orig.events.offset},...
      'duration', {orig.events.duration});
    
  case 'ced_spike6mat'
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    
    chanindx = [];
    for i = 1:numel(hdr.orig)
      if ~any(isfield(hdr.orig{i}, {'units', 'scale'}))
        chanindx = [chanindx i];
      end
    end
    if ~isempty(chanindx)
      trigger = read_trigger(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', chanindx, 'detectflank', detectflank);
      event   = appendevent(event, trigger);
    end
    
  case {'ctf_ds', 'ctf_meg4', 'ctf_res4', 'ctf_old'}
    % obtain the dataset name
    if ft_filetype(filename, 'ctf_meg4') ||  ft_filetype(filename, 'ctf_res4')
      filename = fileparts(filename);
    end
    [path, name, ext] = fileparts(filename);
    headerfile = fullfile(path, [name ext], [name '.res4']);
    datafile   = fullfile(path, [name ext], [name '.meg4']);
    classfile  = fullfile(path, [name ext], 'ClassFile.cls');
    markerfile = fullfile(path, [name ext], 'MarkerFile.mrk');
    
    % in case ctf_old was specified as eventformat, the other reading functions should also know about that
    if strcmp(eventformat, 'ctf_old')
      dataformat   = 'ctf_old';
      headerformat = 'ctf_old';
    end
    
    % read the header, required to determine the stimulus channels and trial specification
    if isempty(hdr)
      hdr = ft_read_header(headerfile, 'headerformat', headerformat);
    end
    
    try
      % read the trigger codes from the STIM channel, usefull for (pseudo) continuous data
      % this splits the trigger channel into the lowers and highest 16 bits,
      % corresponding with the front and back panel of the electronics cabinet at the Donders Centre
      [backpanel, frontpanel] = read_ctf_trigger(filename);
      for i=find(backpanel(:)')
        event(end+1).type   = 'backpanel trigger';
        event(end  ).sample = i;
        event(end  ).value  = backpanel(i);
      end
      for i=find(frontpanel(:)')
        event(end+1).type   = 'frontpanel trigger';
        event(end  ).sample = i;
        event(end  ).value  = frontpanel(i);
      end
    end
    
    % determine the trigger channels from the header
    if isfield(hdr, 'orig') && isfield(hdr.orig, 'sensType')
      origSensType = hdr.orig.sensType;
    elseif isfield(hdr, 'orig') && isfield(hdr.orig, 'res4')
      origSensType =  [hdr.orig.res4.senres.sensorTypeIndex];
    else
      origSensType = [];
    end
    % meg channels are 5, refmag 0, refgrad 1, adcs 18, trigger 11 or 20, eeg 9
    trigindx = find(origSensType==11 | origSensType==20);
    if ~isempty(trigindx)
      % read the trigger channel and do flank detection
      trigger = read_trigger(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', trigindx, 'detectflank', detectflank, 'trigshift', trigshift, 'fixctf', true);
      event   = appendevent(event, trigger);
    end
    
    % read the classification file and make an event for each classified trial
    [condNumbers,condLabels] = read_ctf_cls(classfile);
    if ~isempty(condNumbers)
      Ncond = length(condLabels);
      for i=1:Ncond
        for j=1:length(condNumbers{i})
          event(end+1).type     = 'classification';
          event(end  ).value    = condLabels{i};
          event(end  ).sample   = (condNumbers{i}(j)-1)*hdr.nSamples + 1;
          event(end  ).offset   = -hdr.nSamplesPre;
          event(end  ).duration =  hdr.nSamples;
        end
      end
    end
    
    if exist(markerfile,'file')
      % read the marker file and make an event for each marker
      % this depends on the readmarkerfile function that I got from Tom Holroyd
      % I have not tested this myself extensively, since at the FCDC we
      % don't use the marker files
      mrk = readmarkerfile(filename);
      for i=1:mrk.number_markers
        for j=1:mrk.number_samples(i)
          % determine the location of the marker, expressed in samples
          trialnum = mrk.trial_times{i}(j,1);
          synctime = mrk.trial_times{i}(j,2);
          begsample = (trialnum-1)*hdr.nSamples + 1;    % of the trial, relative to the start of the datafile
          endsample = (trialnum  )*hdr.nSamples;        % of the trial, relative to the start of the datafile
          offset    = round(synctime*hdr.Fs);           % this is the offset (in samples) relative to time t=0 for this trial
          offset    = offset + hdr.nSamplesPre;         % and time t=0 corrsponds with the nSamplesPre'th sample
          % store this marker as an event
          event(end+1).type    = mrk.marker_names{i};
          event(end ).value    = [];
          event(end ).sample   = begsample + offset;
          event(end ).duration = 0;
          event(end ).offset   = offset;
        end
      end
    end
    
  case 'ctf_shm'
    % contact Robert Oostenveld if you are interested in real-time acquisition on the CTF system
    % read the events from shared memory
    event = read_shm_event(filename, varargin{:});
    
  case 'dataq_wdq'
    if isempty(hdr)
      hdr     = ft_read_header(filename, 'headerformat', 'dataq_wdq');
    end
    trigger  = read_wdq_data(filename, hdr.orig, 'lowbits');
    [ix, iy] = find(trigger>1); %it seems as if the value of 1 is meaningless
    for i=1:numel(ix)
      event(i).type   = num2str(ix(i));
      event(i).value  = trigger(ix(i),iy(i));
      event(i).sample = iy(i);
    end
    
  case 'edf'
    % read the header
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    
    if ~isempty(detectflank) % parse the trigger channel (indicated by chanindx) for events
      event = read_trigger(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', chanindx, 'detectflank', detectflank, 'trigshift', trigshift, 'threshold', threshold);
    elseif issubfield(hdr, 'orig.annotation') && ~isempty(hdr.orig.annotation) % EDF itself does not contain events, but EDF+ does define an annotation channel
      % read the data of the annotation channel as 16 bit
      evt = read_edf(filename, hdr);
      % undo the faulty calibration
      evt = (evt - hdr.orig.Off(hdr.orig.annotation)) ./ hdr.orig.Cal(hdr.orig.annotation);
      % convert the 16 bit format into the separate bytes
      evt = typecast(int16(evt), 'uint8');
      % construct the Time-stamped Annotations Lists (TAL), see http://www.edfplus.info/specs/edfplus.html#tal
      tal  = tokenize(char(evt), char(0), true);
      
      event = [];
      for i=1:length(tal)
        % the unprintable characters 20 and 21 are used as separators between time, duration and the annotation
        % duration can be skipped in which case its preceding 21 must also be skipped
        
        tok = tokenize(tal{i}, char(20)); % split time and annotation
        if any(tok{1}==21)
          % the time and duration are specified
          dum = tokenize(tok{1}, char(21)); % split time and duration
          time     = str2double(dum{1});
          duration = str2double(dum{2});
        else
          % only the time is specified
          time     = str2double(tok{1});
          duration = [];
        end
        % there can be multiple annotations per time, the last cell is always empty
        for j=2:length(tok)-1
          anot = char(tok{j});
          % represent the annotation as event
          event(end+1).type     = 'annotation';
          event(end ).value     = anot;
          event(end ).sample    = round(time*hdr.Fs + 1); % expressed in samples, first sample in the file is 1
          event(end ).duration  = round(duration*hdr.Fs); % expressed in samples
          event(end ).timestamp = time; % in seconds, relative to the start of the recording
          event(end ).offset    = 0;
        end
      end
    else
      event = [];
    end
    
  case 'eeglab_set'
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    event = read_eeglabevent(filename, 'header', hdr);
    
  case 'eeglab_erp'
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    event = read_erplabevent(filename, 'header', hdr);
    
  case 'eep_avr'
    % check that the required low-level toolbox is available
    ft_hastoolbox('eeprobe', 1);
    % the headerfile and datafile are the same
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    event(end+1).type     = 'average';
    event(end  ).sample   = 1;
    event(end  ).duration = hdr.nSamples;
    event(end  ).offset   = -hdr.nSamplesPre;
    event(end  ).value    = [];
    
  case {'eep_cnt' 'eep_trg'}
    % check that the required low-level toolbox is available
    ft_hastoolbox('eeprobe', 1);
    
    % this requires the header from the cnt file and the triggers from the trg file
    if strcmp(eventformat, 'eep_cnt')
      trgfile = [filename(1:(end-3)), 'trg'];
      cntfile = filename;
    elseif strcmp(eventformat, 'eep_trg')
      cntfile = [filename(1:(end-3)), 'cnt'];
      trgfile = filename;
    end
    
    if exist(trgfile, 'file')
      trg = read_eep_trg(trgfile);
    else
      warning('The corresponding "%s" file was not found, cannot read in trigger information. No events can be read in.', trgfile);
      trg = []; % make it empty, needed below
    end
    
    if isempty(hdr)
      if exist(cntfile, 'file')
        hdr = ft_read_header(cntfile);
      else
        warning('The corresponding "%s" file was not found, cannot read in header information. No events can be read in.', cntfile);
        hdr = []; % remains empty, needed below
      end
    end
    
    if ~isempty(trg) && ~isempty(hdr)
      % translate the EEProbe trigger codes to events
      for i=1:length(trg)
        event(i).type     = 'trigger';
        event(i).sample   = trg(i).offset;
        event(i).value    = trg(i).type;
        event(i).offset   = 0;
        event(i).duration = 0;
      end
    end
    
  case 'egi_egis'
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    fhdr   = hdr.orig.fhdr;
    chdr   = hdr.orig.chdr;
    ename  = hdr.orig.ename;
    cnames = hdr.orig.cnames;
    fcom   = hdr.orig.fcom;
    ftext  = hdr.orig.ftext;
    eventCount=0;
    for cel=1:fhdr(18)
      for trial=1:chdr(cel,2)
        eventCount=eventCount+1;
        event(eventCount).type     = 'trial';
        event(eventCount).sample   = (eventCount-1)*hdr.nSamples + 1;
        event(eventCount).offset   = -hdr.nSamplesPre;
        event(eventCount).duration =  hdr.nSamples;
        event(eventCount).value    =  cnames{cel};
      end
    end
    
  case 'egi_egia'
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    fhdr   = hdr.orig.fhdr;
    chdr   = hdr.orig.chdr;
    ename  = hdr.orig.ename;
    cnames = hdr.orig.cnames;
    fcom   = hdr.orig.fcom;
    ftext  = hdr.orig.ftext;
    eventCount=0;
    for cel=1:fhdr(18)
      for subject=1:chdr(cel,2)
        eventCount=eventCount+1;
        event(eventCount).type     = 'trial';
        event(eventCount).sample   = (eventCount-1)*hdr.nSamples + 1;
        event(eventCount).offset   = -hdr.nSamplesPre;
        event(eventCount).duration =  hdr.nSamples;
        event(eventCount).value    =  ['Sub' sprintf('%03d',subject) cnames{cel}];
      end
    end
    
  case 'egi_sbin'
    if ~exist('segHdr','var')
      [EventCodes, segHdr, eventData] = read_sbin_events(filename);
    end
    if ~exist('header_array','var')
      [header_array, CateNames, CatLengths, preBaseline] = read_sbin_header(filename);
    end
    if isempty(hdr)
      hdr = ft_read_header(filename,'headerformat','egi_sbin');
    end
    version     = header_array(1);
    unsegmented = ~mod(version, 2);
    
    eventCount=0;
    if unsegmented
      [evType,sampNum] = find(eventData);
      for k = 1:length(evType)
        event(k).sample   = sampNum(k);
        event(k).offset   = [];
        event(k).duration = 0;
        event(k).type     = 'trigger';
        event(k).value    = char(EventCodes(evType(k),:));
      end
    else
      for theEvent=1:size(eventData,1)
        for segment=1:hdr.nTrials
          if any(eventData(theEvent,((segment-1)*hdr.nSamples +1):segment*hdr.nSamples))
            eventCount=eventCount+1;
            event(eventCount).sample   = min(find(eventData(theEvent,((segment-1)*hdr.nSamples +1):segment*hdr.nSamples))) +(segment-1)*hdr.nSamples;
            event(eventCount).offset   = -hdr.nSamplesPre;
            event(eventCount).duration =  length(find(eventData(theEvent,((segment-1)*hdr.nSamples +1):segment*hdr.nSamples )>0))-1;
            event(eventCount).type     = 'trigger';
            event(eventCount).value    =  char(EventCodes(theEvent,:));
          end
        end
      end
    end
    
    if ~unsegmented
      for segment=1:hdr.nTrials  % cell information
        eventCount=eventCount+1;
        event(eventCount).type     = 'trial';
        event(eventCount).sample   = (segment-1)*hdr.nSamples + 1;
        event(eventCount).offset   = -hdr.nSamplesPre;
        event(eventCount).duration =  hdr.nSamples;
        event(eventCount).value    =  char([CateNames{segHdr(segment,1)}(1:CatLengths(segHdr(segment,1)))]);
      end
    end;
    
  case {'egi_mff_v1' 'egi_mff'} % this is currently the default
    % The following represents the code that was written by Ingrid, Robert
    % and Giovanni to get started with the EGI mff dataset format. It might
    % not support all details of the file formats.
    % An alternative implementation has been provided by EGI, this is
    % released as fieldtrip/external/egi_mff and referred further down in
    % this function as 'egi_mff_v2'.
    
    if isempty(hdr)
      % use the corresponding code to read the header
      hdr = ft_read_header(filename, 'headerformat', eventformat);
    end
    
    if ~usejava('jvm')
      error('the xml2struct requires MATLAB to be running with the Java virtual machine (JVM)');
      % an alternative implementation which does not require the JVM but runs much slower is
      % available from http://www.mathworks.com/matlabcentral/fileexchange/6268-xml4mat-v2-0
    end
    
    % get event info from xml files
    warning('off', 'MATLAB:REGEXP:deprecated') % due to some small code xml2struct
    xmlfiles = dir( fullfile(filename, '*.xml'));
    disp('reading xml files to obtain event info... This might take a while if many events/triggers are present')
    if isempty(xmlfiles)
      xml=struct([]);
    else
      xml=[];
    end;
    for i = 1:numel(xmlfiles)
      if strcmpi(xmlfiles(i).name(1:6), 'Events')
        fieldname       = strrep(xmlfiles(i).name(1:end-4), ' ', '_');
        filename_xml    = fullfile(filename, xmlfiles(i).name);
        xml.(fieldname) = xml2struct(filename_xml);
      end
    end
    warning('on', 'MATLAB:REGEXP:deprecated')
    
    % construct info needed for FieldTrip Event
    eventNames = fieldnames(xml);
    begTime = hdr.orig.xml.info.recordTime;
    begTime(11) = ' '; begTime(end-5:end) = [];
    begSDV = datenum(begTime);
    % find out if there are epochs in this dataset
    if isfield(hdr.orig.xml,'epochs') && length(hdr.orig.xml.epochs) > 1
      Msamp2offset = nan(2,size(hdr.orig.epochdef,1),1+max(hdr.orig.epochdef(:,2)-hdr.orig.epochdef(:,1)));
      for iEpoch = 1:size(hdr.orig.epochdef,1)
        nSampEpoch = hdr.orig.epochdef(iEpoch,2)-hdr.orig.epochdef(iEpoch,1)+1;
        Msamp2offset(1,iEpoch,1:nSampEpoch) = hdr.orig.epochdef(iEpoch,1):hdr.orig.epochdef(iEpoch,2); %sample number in samples
        Msamp2offset(2,iEpoch,1:nSampEpoch) = hdr.orig.epochdef(iEpoch,3):hdr.orig.epochdef(iEpoch,3)+nSampEpoch-1; %offset in samples
      end
    end
    
    % construct event according to FieldTrip rules
    eventCount = 0;
    for iXml = 1:length(eventNames)
      eval(['eventField=isfield(xml.' eventNames{iXml} ',''event'');'])
      if eventField
        for iEvent = 1:length(xml.(eventNames{iXml}))
          eventTime  = xml.(eventNames{iXml})(iEvent).event.beginTime;
          eventTime(11) = ' '; eventTime(end-5:end) = [];
          if strcmp('-',eventTime(21))
            % event out of range (before recording started): do nothing.
          else
            eventSDV = datenum(eventTime);
            eventOffset = round((eventSDV - begSDV)*24*60*60*hdr.Fs); %in samples, relative to start of recording
            if eventOffset < 0
              % event out of range (before recording started): do nothing
            else
              % calculate eventSample, relative to start of epoch
              if isfield(hdr.orig.xml,'epochs') && length(hdr.orig.xml.epochs) > 1
                SampIndex=[];
                for iEpoch = 1:size(hdr.orig.epochdef,1)
                  [dum,dum2] = intersect(squeeze(Msamp2offset(2,iEpoch,:)), eventOffset);
                  if ~isempty(dum2)
                    EpochNum = iEpoch;
                    SampIndex = dum2;
                  end
                end
                if ~isempty(SampIndex)
                  eventSample = Msamp2offset(1,EpochNum,SampIndex);
                else
                  eventSample=[]; %Drop event if past end of epoch
                end;
              else
                eventSample = eventOffset+1;
              end
              if ~isempty(eventSample)
                eventCount=eventCount+1;
                event(eventCount).type     = eventNames{iXml}(8:end);
                event(eventCount).sample   = eventSample;
                event(eventCount).offset   = 0;
                event(eventCount).duration = str2double(xml.(eventNames{iXml})(iEvent).event.duration)./1000000000*hdr.Fs;
                event(eventCount).value    = xml.(eventNames{iXml})(iEvent).event.code;
                event(eventCount).orig    = xml.(eventNames{iXml})(iEvent).event;
              end;
            end  %if that takes care of non "-" events that are still out of range
          end %if that takes care of "-" events, which are out of range
        end %iEvent
      end;
    end
    
    % add "epoch" events for epoched data, i.e. data with variable length segments
    if (hdr.nTrials==1)
      eventCount=length(event);
      for iEpoch=1:size(hdr.orig.epochdef,1)
        eventCount=eventCount+1;
        event(eventCount).type     = 'epoch';
        event(eventCount).sample   = hdr.orig.epochdef(iEpoch,1);
        event(eventCount).duration = hdr.orig.epochdef(iEpoch,2)-hdr.orig.epochdef(iEpoch,1)+1;
        event(eventCount).offset   = hdr.orig.epochdef(iEpoch,3);
        event(eventCount).value    = [];
      end % for
    end % if
    
    % add "trial" events for segmented data, i.e. data with constant length segments
    if (hdr.nTrials >1) && (size(hdr.orig.epochdef,1)==hdr.nTrials)
      cellNames=cell(hdr.nTrials,1);
      recTime=zeros(hdr.nTrials,1);
      epochTimes=cell(hdr.nTrials,1);
      for iEpoch=1:hdr.nTrials
        epochTimes{iEpoch}=hdr.orig.xml.epochs(iEpoch).epoch.beginTime;
        recTime(iEpoch)=((str2num(hdr.orig.xml.epochs(iEpoch).epoch.beginTime)/1000)/(1000/hdr.Fs))+1;
      end
      for iCat=1:length(hdr.orig.xml.categories)
        theTimes=cell(length(hdr.orig.xml.categories(iCat).cat.segments),1);
        for i=1:length(hdr.orig.xml.categories(iCat).cat.segments)
          theTimes{i}=hdr.orig.xml.categories(iCat).cat.segments(i).seg.beginTime;
        end
        epochIndex=find(ismember(epochTimes,theTimes));
        for i=1:length(epochIndex)
          cellNames{epochIndex(i)}=hdr.orig.xml.categories(iCat).cat.name;
        end
      end
      
      eventCount=length(event);
      for iEpoch=1:hdr.nTrials
        eventCount=eventCount+1;
        event(eventCount).type     = 'trial';
        event(eventCount).sample   = hdr.orig.epochdef(iEpoch,1);
        event(eventCount).offset   = -hdr.nSamplesPre;
        event(eventCount).duration = hdr.nSamples;
        event(eventCount).value    = cellNames{iEpoch};
      end
    end
    
  case 'egi_mff_v2'
    % ensure that the EGI toolbox is on the path
    ft_hastoolbox('egi_mff', 1);
    if isunix && filename(1)~=filesep
      % add the full path to the dataset directory
      filename = fullfile(pwd, filename);
    elseif ispc && filename(2)~=':'
      % add the full path, including drive letter
      filename = fullfile(pwd, filename);
    end
    % pass the header along to speed it up, it will be read on the fly in case it is empty
    event = read_mff_event(filename, hdr);
    % clean up the fields in the event structure
    fn = fieldnames(event);
    fn = setdiff(fn, {'type', 'sample', 'value', 'offset', 'duration', 'timestamp'});
    for i=1:length(fn)
      event = rmfield(event, fn{i});
    end
    
  case 'smi_txt'
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    if isfield(hdr.orig, 'trigger')
      % this is inefficient, since it keeps the complete data in memory
      % but it does speed up subsequent read operations without the user
      % having to care about it
      smi = hdr.orig;
    else
      smi = read_smi_txt(filename);
    end
    timestamp = [smi.trigger(:).timestamp];
    value     = [smi.trigger(:).value];
    % note that in this dataformat the first input trigger can be before
    % the start of the data acquisition
    for i=1:length(timestamp)
      event(end+1).type       = 'Trigger';
      event(end  ).sample     = (timestamp(i)-hdr.FirstTimeStamp)/hdr.TimeStampPerSample + 1;
      event(end  ).timestamp  = timestamp(i);
      event(end  ).value      = value(i);
      event(end  ).duration   = 1;
      event(end  ).offset     = 0;
    end
    
  case 'eyelink_asc'
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    if isfield(hdr.orig, 'input')
      % this is inefficient, since it keeps the complete data in memory
      % but it does speed up subsequent read operations without the user
      % having to care about it
      asc = hdr.orig;
    else
      asc = read_eyelink_asc(filename);
    end
    timestamp = [asc.input(:).timestamp];
    value     = [asc.input(:).value];
    % note that in this dataformat the first input trigger can be before
    % the start of the data acquisition
    for i=1:length(timestamp)
      event(end+1).type       = 'INPUT';
      event(end  ).sample     = (timestamp(i)-hdr.FirstTimeStamp)/hdr.TimeStampPerSample + 1;
      event(end  ).timestamp  = timestamp(i);
      event(end  ).value      = value(i);
      event(end  ).duration   = 1;
      event(end  ).offset     = 0;
    end
    
  case 'fcdc_global'
    event = event_queue;
    
  case 'fcdc_buffer'
    % read from a networked buffer for realtime analysis
    [host, port] = filetype_check_uri(filename);
    
    % SK: the following was intended to speed up, but does not work
    % the buffer server will try to return exact indices, even
    % if the intend here is to filter based on a maximum range.
    % We could change the definition of GET_EVT to comply
    % with filtering, but that might break other existing code.
    
    %if isempty(flt_minnumber) && isempty(flt_maxnumber)
    %   evtsel = [];
    %else
    %   evtsel = [0 2^32-1];
    %   if ~isempty(flt_minnumber)
    %       evtsel(1) = flt_minnumber-1;
    %   end
    %   if ~isempty(flt_maxnumber)
    %       evtsel(2) = flt_maxnumber-1;
    %   end
    %end
    
    if blocking && isempty(flt_minnumber) && isempty(flt_maxnumber)
      warning('disabling blocking because no selection was specified');
      blocking = false;
    end
    
    if blocking
      nsamples = 0; % disable waiting for samples
      if isempty(flt_minnumber)
        nevents = flt_maxnumber;
      elseif isempty(flt_maxnumber)
        nevents = flt_minnumber;
      else
        nevents = max(flt_minnumber, flt_maxnumber);
      end
      available = buffer_wait_dat([nsamples nevents timeout], host, port);
      if available.nevents<nevents
        error('buffer timed out while waiting for %d events', nevents);
      end
    end
    
    try
      event = buffer('get_evt', [], host, port);
    catch
      if strfind(lasterr, 'the buffer returned an error')
        % this happens if the buffer contains no events
        % which in itself is not a problem and should not result in an error
        event = [];
      else
        rethrow(lasterr);
      end
    end
    
  case 'fcdc_buffer_offline'
    if isdir(filename)
      path = filename;
    else
      [path, file, ext] = fileparts(filename);
    end
    if isempty(hdr)
      headerfile = fullfile(path, 'header');
      hdr = read_buffer_offline_header(headerfile);
    end
    eventfile = fullfile(path, 'events');
    event = read_buffer_offline_events(eventfile, hdr);
    
  case 'fcdc_matbin'
    % this is multiplexed data in a *.bin file, accompanied by a MATLAB file containing the header and event
    [path, file, ext] = fileparts(filename);
    filename = fullfile(path, [file '.mat']);
    % read the events from the MATLAB file
    tmp   = load(filename, 'event');
    if isfield(tmp, 'event')
      event = tmp.event;
    else
      event = [];
    end
    
  case 'fcdc_fifo'
    fifo = filetype_check_uri(filename);
    
    if ~exist(fifo,'file')
      warning('the FIFO %s does not exist; attempting to create it', fifo);
      fid = fopen(fifo, 'r');
      system(sprintf('mkfifo -m 0666 %s',fifo));
    end
    
    msg = fread(fid, inf, 'uint8');
    fclose(fid);
    
    try
      event = mxDeserialize(uint8(msg));
    catch
      warning(lasterr);
    end
    
  case 'fcdc_tcp'
    % requires tcp/udp/ip-toolbox
    ft_hastoolbox('TCP_UDP_IP', 1);
    [host, port] = filetype_check_uri(filename);
    if isempty(sock)
      sock=pnet('tcpsocket',port);
    end
    con = pnet(sock, 'tcplisten');
    if con~=-1
      try
        pnet(con,'setreadtimeout',10);
        % read packet
        msg=pnet(con,'readline'); %,1000,'uint8','network');
        if ~isempty(msg)
          event = mxDeserialize(uint8(str2num(msg)));
        end
        %       catch
        %         warning(lasterr);
      end
      pnet(con,'close');
    end
    con = [];
    
  case 'fcdc_udp'
    % requires tcp/udp/ip-toolbox
    ft_hastoolbox('TCP_UDP_IP', 1);
    [host, port] = filetype_check_uri(filename);
    try
      % read from localhost
      udp=pnet('udpsocket',port);
      % Wait/Read udp packet to read buffer
      len=pnet(udp,'readpacket');
      if len>0,
        % if packet larger then 1 byte then read maximum of 1000 doubles in network byte order
        msg=pnet(udp,'read',1000,'uint8');
        if ~isempty(msg)
          event = mxDeserialize(uint8(msg));
        end
      end
    catch
      warning(lasterr);
    end
    % On break or error close connection
    pnet(udp,'close');
    
  case 'fcdc_serial'
    % this code is moved to a separate file
    event = read_serial_event(filename);
    
  case 'fcdc_mysql'
    % check that the required low-level toolbox is available
    ft_hastoolbox('mysql', 1);
    % read from a MySQL server listening somewhere else on the network
    db_open(filename);
    if db_blob
      event = db_select_blob('fieldtrip.event', 'msg');
    else
      event = db_select('fieldtrip.event', {'type', 'value', 'sample', 'offset', 'duration'});
    end
    
  case 'gtec_hdf5'
    % the header mentions trigger channels, but I don't know how they are stored
    warning('event reading for hdf5 has not yet been implemented due to a lack of a good example file');
    
  case 'gtec_mat'
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    if isempty(trigindx)
      % these are probably trigger channels
      trigindx = match_str(hdr.label, {'Display', 'Target'});
    end
    % use a helper function to read the trigger channels and detect the flanks
    % pass all the other users options to the read_trigger function
    event = read_trigger(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', trigindx, 'detectflank', detectflank, 'trigshift', trigshift);
    
  case {'homer_nirs'}
    % Homer files are MATLAB files in disguise
    orig = load(filename, '-mat');
    % each of the columns of orig.s represents an event type
    % negative values have been editted in Homer and should be ignored
    event = [];
    for i=1:size(orig.s,2)
      smp = find(orig.s(:,i)==1);
      for j=1:numel(smp)
        event(end+1).type     = 'trigger';
        event(end  ).value    = i;
        event(end  ).sample   = smp(j);
        event(end  ).duration = [];
        event(end  ).offset   = [];
      end
    end
    
  case {'itab_raw' 'itab_mhd'}
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    for i=1:hdr.orig.nsmpl
      event(end+1).type = 'trigger';
      event(end  ).value    = hdr.orig.smpl(i).type;
      event(end  ).sample   = hdr.orig.smpl(i).start + 1;
      event(end  ).duration = hdr.orig.smpl(i).ntptot;
      event(end  ).offset   = -hdr.orig.smpl(i).ntppre;  % number of samples prior to the trigger
    end
    if isempty(event)
      warning('no events found in the event table, reading the trigger channel(s)');
      trigindx = find(ft_chantype(hdr, 'flag'));
      trigger = read_trigger(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', trigindx, 'detectflank', detectflank, 'trigshift', trigshift);
      event   = appendevent(event, trigger);
    end
    
  case 'matlab'
    % read the events from a normal MATLAB file
    tmp   = load(filename, 'event');
    event = tmp.event;
    
  case {'manscan_mbi', 'manscan_mb2'}
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    if isfield(hdr.orig, 'epochs') && ~isempty(hdr.orig.epochs)
      trlind = [];
      for i = 1:numel(hdr.orig.epochs)
        trlind = [trlind i*ones(1, diff(hdr.orig.epochs(i).samples) + 1)];
      end
    else
      trlind = ones(1, hdr.nSamples);
    end
    if isfield(hdr.orig, 'events')
      for i = 1:numel(hdr.orig.events)
        for j = 1:length(hdr.orig.events(i).samples)
          event(end+1).type   = 'trigger';
          event(end).value    = hdr.orig.events(i).label;
          event(end).sample   = find(cumsum(trlind == hdr.orig.events(i).epochs(j))...
            == hdr.orig.events(i).samples(j), 1, 'first');
        end
      end
    end
    
  case 'mega_neurone'
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    % this is fast but memory inefficient, since the header contains all data and events
    if isfield(hdr.orig, 'event')
      NEURONE = hdr.orig;
    else
      % ensure that this external toolbox is on the path
      ft_hastoolbox('neurone', 1);
      if filename(end)~=filesep
        % it should end with a slash
        filename = [filename filesep];
      end
      NEURONE = readneurone(filename);
    end
    for i=1:numel(NEURONE.event)
      if isnan(str2double(NEURONE.event(i).type))
        % there are a number of event "Types" that can happen on different "SourcePorts"
        event(i).type     = NEURONE.event(i).type;
        event(i).sample   = round(NEURONE.event(i).latency * 0.001 * NEURONE.srate + 1);
        event(i).value    = [];
      else
        % this seems to correspond with an external trigger code, which is best represented numerically
        event(i).type     = 'trigger';
        event(i).sample   = round(NEURONE.event(i).latency * 0.001 * NEURONE.srate + 1);
        event(i).value    = str2double(NEURONE.event(i).type);
      end
    end
    
  case 'micromed_trc'
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    if isfield(hdr, 'orig') && isfield(hdr.orig, 'Trigger_Area') && isfield(hdr.orig, 'Tigger_Area_Length')
      if ~isempty(trigindx)
        trigger = read_trigger(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', trigindx, 'detectflank', detectflank, 'trigshift', trigshift);
        event   = appendevent(event, trigger);
      else
        % routine that reads analog triggers in case no index is specified
        event = read_micromed_event(filename);
      end
      
    else
      error('Not a correct event format')
    end
    
  case {'mpi_ds', 'mpi_dap'}
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    % determine the DAP files that compromise this dataset
    if isdir(filename)
      ls = dir(filename);
      dapfile = {};
      for i=1:length(ls)
        if ~isempty(regexp(ls(i).name, '.dap$', 'once' ))
          dapfile{end+1} = fullfile(filename, ls(i).name);
        end
      end
      dapfile = sort(dapfile);
    elseif iscell(filename)
      dapfile = filename;
    else
      dapfile = {filename};
    end
    % assume that each DAP file is accompanied by a dat file
    % read the trigger values from the separate dat files
    trg = [];
    for i=1:length(dapfile)
      datfile = [dapfile{i}(1:(end-4)) '.dat'];
      trg = cat(1, trg, textread(datfile, '', 'headerlines', 1));
    end
    % construct a event structure, one 'trialcode' event per trial
    for i=1:length(trg)
      event(i).type     = 'trialcode';            % string
      event(i).sample   = (i-1)*hdr.nSamples + 1; % expressed in samples, first sample of file is 1
      event(i).value    = trg(i);                 % number or string
      event(i).offset   = 0;                      % expressed in samples
      event(i).duration = hdr.nSamples;           % expressed in samples
    end
    
  case 'nervus_eeg'
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    % construct a event structure from data in the header
    maxSampleRate = max([hdr.orig.Segments.samplingRate]);
    earliestDateTime = min([hdr.orig.Segments.dateOLE]);
    for i=1:length(hdr.orig.Events)
      event(i).type     = hdr.orig.Events(i).IDStr;   % string
      event(i).value    = hdr.orig.Events(i).label;  % number or string
      event(i).offset   = 0;                         % expressed in samples
      % calculate the sample value of the event, based on the highest
      % sample rate
      event(i).sample   = (hdr.orig.Events(i).dateOLE-earliestDateTime)*3600*24*maxSampleRate;
      if event(i).sample == 0
        event(i).sample = 1;
      elseif event(i).sample > hdr.nSamples
        event(i).sample = hdr.nSamples;
      end
      event(i).duration = hdr.orig.Events(i).duration*maxSampleRate;
    end
    %Add boundary events to indicate segments
    originalEventCount = length(hdr.orig.Events);
    boundaryEventCount = 1;
    for i=2:length(hdr.orig.Segments)
      event(originalEventCount+boundaryEventCount).type = 'boundary';
      event(originalEventCount+boundaryEventCount).value = 'boundary';
      event(originalEventCount+boundaryEventCount).offset = 0;
      event(originalEventCount+boundaryEventCount).duration = 0;
      event(originalEventCount+boundaryEventCount).sample = sum([hdr.orig.Segments(1:(i-1)).sampleCount]);
      boundaryEventCount = boundaryEventCount+1;
    end
    
  case {'neuromag_eve'}
    % previously this was called babysquid_eve, now it is neuromag_eve
    % see also http://bugzilla.fcdonders.nl/show_bug.cgi?id=2170
    [p, f, x] = fileparts(filename);
    evefile = fullfile(p, [f '.eve']);
    fiffile = fullfile(p, [f '.fif']);
    
    if isempty(hdr)
      hdr = ft_read_header(fiffile, 'headerformat', 'neuromag_mne');
    end
    
    [smp, tim, typ, val] = read_neuromag_eve(evefile);
    
    smp = smp + 1; % should be 1-offset
    smp = smp - double(hdr.orig.raw.first_samp); % the recording to disk may start later than the actual acquisition
    
    if ~isempty(smp)
      sample  = num2cell(smp);
      value   = num2cell(val);
      offset  = num2cell(zeros(size(smp)));
      type    = repmat({'unknown'}, size(typ));
      type(typ==0) = {'trigger'};
      if any(typ~=0)
        % see the comments in read_neuromag_eve
        warning('entries in the *.eve file with a type other than 0 are represented as ''unknown''')
      end
      % convert to a structure array
      event = struct('type', type, 'value', value, 'sample', sample, 'offset', offset);
    else
      event = [];
    end
    
  case {'neuromag_fif' 'neuromag_mne' 'neuromag_mex'}
    if strcmp(eventformat, 'neuromag_fif')
      % the default is to use the MNE reader for fif files
      eventformat = 'neuromag_mne';
    end
    if strcmp(eventformat, 'neuromag_mex')
      % check that the required low-level toolbox is available
      ft_hastoolbox('meg-pd', 1);
      if isempty(headerformat), headerformat = eventformat; end
      if isempty(dataformat),   dataformat   = eventformat; end
    elseif strcmp(eventformat, 'neuromag_mne')
      % check that the required low-level toolbox is available
      ft_hastoolbox('mne', 1);
      if isempty(headerformat), headerformat = eventformat; end
      if isempty(dataformat),   dataformat   = eventformat; end
    end
    
    if isempty(hdr)
      hdr = ft_read_header(filename, 'headerformat', headerformat, 'checkmaxfilter', checkmaxfilter);
    end
    
    % note below we've had to include some chunks of code that are only
    % called if the file is an averaged file, or if the file is continuous.
    % These are defined in hdr by read_header for neuromag_mne, but do not
    % exist for neuromag_fif, hence we run the code anyway if the fields do
    % not exist (this is what happened previously anyway).
    
    if strcmp(eventformat, 'neuromag_mex')
      iscontinuous    = 1;
      isaverage       = 0;
      isepoched       = 0;
    elseif strcmp(eventformat, 'neuromag_mne')
      iscontinuous    = hdr.orig.iscontinuous;
      isaverage       = hdr.orig.isaverage;
      isepoched       = hdr.orig.isepoched;
    end
    
    if iscontinuous
      analogindx = find(strcmp(ft_chantype(hdr), 'analog trigger'));
      otherindx  = find(strcmp(ft_chantype(hdr), 'other trigger'));
      binaryindx = find(strcmp(ft_chantype(hdr), 'digital trigger'));
      
      if isempty(binaryindx) && isempty(analogindx) && isempty(otherindx)
        % included in case of problems with older systems and MNE reader:
        % use a predefined set of channel names
        binary     = {'STI 014', 'STI 015', 'STI 016'};
        binaryindx = match_str(hdr.label, binary);
      end
      
      if ~isempty(binaryindx)
        trigger = read_trigger(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', binaryindx, 'detectflank', detectflank, 'trigshift', trigshift, 'fixneuromag', false);
        event   = appendevent(event, trigger);
      end
      if ~isempty(otherindx)
        trigger = read_trigger(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', otherindx, 'detectflank', detectflank, 'trigshift', trigshift, 'fixneuromag', false);
        event   = appendevent(event, trigger);
      end
      if ~isempty(analogindx)
        % add the triggers to the event structure based on trigger channels with the name "STI xxx"
        % there are some issues with noise on these analog trigger
        % channels, on older systems only
        % read the trigger channel and do flank detection
        trigger = read_trigger(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', analogindx, 'detectflank', detectflank, 'trigshift', trigshift, 'fixneuromag', true);
        event   = appendevent(event, trigger);
        
        if ~isempty(trigger)
          % throw out everything that is not a (real) trigger...
          [sel1, sel2] = match_str({trigger.type}, {'STI001', 'STI002', 'STI003', 'STI004', 'STI005', 'STI006', 'STI007', 'STI008'});
          trigger_tmp = trigger(sel1);
          %trigger_bits = size(unique(sel2), 1);
          all_triggers = unique(sel2);
          trigger_bits = max(all_triggers) - min(all_triggers) + 1;
          
          % collect all samples where triggers occured...
          all_samples = unique([trigger_tmp.sample]);
          
          new_triggers = [];
          i = 1;
          while i <= length(all_samples)
            cur_triggers = [];
            j = 0;
            
            % look also in adjacent samples according to tolerance
            while (i+j <= length(all_samples)) && (all_samples(i+j) - all_samples(i) <= tolerance)
              if j >= 1
                fprintf('Fixing trigger at sample %d\n', all_samples(i));
              end %if
              cur_triggers = appendevent(cur_triggers, trigger_tmp([trigger_tmp.sample] == all_samples(i+j)));
              j = j + 1;
            end %while
            
            % construct new trigger field
            if ~isempty(cur_triggers)
              new_triggers(end+1).type = 'Trigger';
              new_triggers(end).sample = all_samples(i);
              [sel1, sel2] = match_str({cur_triggers.type}, {'STI001', 'STI002', 'STI003', 'STI004', 'STI005', 'STI006', 'STI007', 'STI008'});
              new_triggers(end).value = bin2dec((dec2bin(sum(2.^(sel2-1)), trigger_bits)));
            end %if
            
            i = i + j;
          end %while
          event = appendevent(event, new_triggers);
        end %if
        
      end
      
    elseif isaverage
      % the length of each average can be variable
      nsamples = zeros(1, length(hdr.orig.evoked));
      for i=1:length(hdr.orig.evoked)
        nsamples(i)  = size(hdr.orig.evoked(i).epochs, 2);
      end
      begsample = cumsum([1 nsamples]);
      for i=1:length(hdr.orig.evoked)
        event(end+1).type     = 'average';
        event(end  ).sample   = begsample(i);
        event(end  ).value    = hdr.orig.evoked(i).comment;  % this is a descriptive string
        event(end  ).offset   = hdr.orig.evoked(i).first;
        event(end  ).duration = hdr.orig.evoked(i).last - hdr.orig.evoked(i).first + 1;
      end
      
    elseif isepoched
      error('Support for epoched *.fif data is not yet implemented.')
    end
    
    % check whether the *.fif file is accompanied by an *.eve file
    [p, f, x] = fileparts(filename);
    evefile = fullfile(p, [f '.eve']);
    if exist(evefile, 'file')
      eve = ft_read_event(evefile, 'header', hdr);
      if ~isempty(eve)
        fprintf('appending %d events from "%s"\n', length(eve), evefile);
        event = appendevent(event, eve);
      end
    end
    
  case {'neuralynx_ttl' 'neuralynx_bin' 'neuralynx_dma' 'neuralynx_sdma'}
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    
    % specify the range to search for triggers, default is the complete file
    if ~isempty(flt_minsample)
      begsample = flt_minsample;
    else
      begsample = 1;
    end
    if ~isempty(flt_maxsample)
      endsample = flt_maxsample;
    else
      endsample = hdr.nSamples*hdr.nTrials;
    end
    
    if strcmp(eventformat, 'neuralynx_dma')
      % read the Parallel_in channel from the DMA log file
      ttl = read_neuralynx_dma(filename, begsample, endsample, 'ttl');
    elseif strcmp(eventformat, 'neuralynx_sdma')
      % determine the separate files with the trigger and timestamp information
      [p, f, x] = fileparts(filename);
      ttlfile = fullfile(filename, [f '.ttl.bin']);
      tslfile = fullfile(filename, [f '.tsl.bin']);
      tshfile = fullfile(filename, [f '.tsh.bin']);
      if ~exist(ttlfile) && ~exist(tslfile) && ~exist(tshfile)
        % perhaps it is an old splitted dma dataset?
        ttlfile = fullfile(filename, [f '.ttl']);
        tslfile = fullfile(filename, [f '.tsl']);
        tshfile = fullfile(filename, [f '.tsh']);
      end
      if ~exist(ttlfile) && ~exist(tslfile) && ~exist(tshfile)
        % these files must be present in a splitted dma dataset
        error('could not locate the individual ttl, tsl and tsh files');
      end
      % read the trigger values from the separate file
      ttl = read_neuralynx_bin(ttlfile, begsample, endsample);
    elseif strcmp(eventformat, 'neuralynx_ttl')
      % determine the optional files with timestamp information
      tslfile = [filename(1:(end-4)) '.tsl'];
      tshfile = [filename(1:(end-4)) '.tsh'];
      % read the triggers from a separate *.ttl file
      ttl = read_neuralynx_ttl(filename, begsample, endsample);
    elseif strcmp(eventformat, 'neuralynx_bin')
      % determine the optional files with timestamp information
      tslfile = [filename(1:(end-8)) '.tsl.bin'];
      tshfile = [filename(1:(end-8)) '.tsh.bin'];
      % read the triggers from a separate *.ttl.bin file
      ttl = read_neuralynx_bin(filename, begsample, endsample);
    end
    
    ttl = int32(ttl / (2^16));    % parallel port provides int32, but word resolution is int16. Shift the bits and typecast to signed integer.
    d1  = (diff(ttl)~=0);         % determine the flanks, which can be multiple samples long (this looses one sample)
    d2  = (diff(d1)==1);          % determine the onset of the flanks (this looses one sample)
    smp = find(d2)+2;             % find the onset of the flanks, add the two samples again
    val = ttl(smp+5);             % look some samples further for the trigger value, to avoid the flank
    clear d1 d2 ttl
    ind = find(val~=0);           % look for triggers tith a non-zero value, this avoids downgoing flanks going to zero
    smp = smp(ind);               % discard triggers with a value of zero
    val = val(ind);               % discard triggers with a value of zero
    
    if ~isempty(smp)
      % try reading the timestamps
      if strcmp(eventformat, 'neuralynx_dma')
        tsl = read_neuralynx_dma(filename, 1, max(smp), 'tsl');
        tsl = typecast(tsl(smp), 'uint32');
        tsh = read_neuralynx_dma(filename, 1, max(smp), 'tsh');
        tsh = typecast(tsh(smp), 'uint32');
        ts  = timestamp_neuralynx(tsl, tsh);
      elseif exist(tslfile) && exist(tshfile)
        tsl = read_neuralynx_bin(tslfile, 1, max(smp));
        tsl = tsl(smp);
        tsh = read_neuralynx_bin(tshfile, 1, max(smp));
        tsh = tsh(smp);
        ts  = timestamp_neuralynx(tsl, tsh);
      else
        ts = [];
      end
      
      % reformat the values as cell array, since the struct function can work with those
      type      = repmat({'trigger'},size(smp));
      value     = num2cell(val);
      sample    = num2cell(smp + begsample - 1);
      duration  = repmat({[]},size(smp));
      offset    = repmat({[]},size(smp));
      if ~isempty(ts)
        timestamp  = reshape(num2cell(ts),size(smp));
      else
        timestamp  = repmat({[]},size(smp));
      end
      % convert it into a structure array, this can be done in one go
      event = struct('type', type, 'value', value, 'sample', sample, 'timestamp', timestamp, 'offset', offset, 'duration', duration);
      clear type value sample timestamp offset duration
    end
    
    if (strcmp(eventformat, 'neuralynx_bin') || strcmp(eventformat, 'neuralynx_ttl')) && isfield(hdr, 'FirstTimeStamp')
      % the header was obtained from an external dataset which could be at a different sampling rate
      % use the timestamps to redetermine the sample numbers
      fprintf('using sample number of the downsampled file to reposition the TTL events\n');
      % convert the timestamps into samples, keeping in mind the FirstTimeStamp and TimeStampPerSample
      smp = round(double(ts - uint64(hdr.FirstTimeStamp))./hdr.TimeStampPerSample + 1);
      for i=1:length(event)
        % update the sample number
        event(i).sample = smp(i);
      end
    end
    
  case 'neuralynx_ds'
    % read the header of the dataset
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    % the event file is contained in the dataset directory
    if     exist(fullfile(filename, 'Events.Nev'))
      filename = fullfile(filename, 'Events.Nev');
    elseif exist(fullfile(filename, 'Events.nev'))
      filename = fullfile(filename, 'Events.nev');
    elseif exist(fullfile(filename, 'events.Nev'))
      filename = fullfile(filename, 'events.Nev');
    elseif exist(fullfile(filename, 'events.nev'))
      filename = fullfile(filename, 'events.nev');
    end
    % read the events, apply filter is applicable
    nev = read_neuralynx_nev(filename, 'type', flt_type, 'value', flt_value, 'mintimestamp', flt_mintimestamp, 'maxtimestamp', flt_maxtimestamp, 'minnumber', flt_minnumber, 'maxnumber', flt_maxnumber);
    
    % the following code should only be executed if there are events,
    % otherwise there will be an error subtracting an uint64 from an []
    if ~isempty(nev)
      % now get the values as cell array, since the struct function can work with those
      value     = {nev.TTLValue};
      timestamp = {nev.TimeStamp};
      number    = {nev.EventNumber};
      type      = repmat({'trigger'},size(value));
      duration  = repmat({[]},size(value));
      offset    = repmat({[]},size(value));
      sample    = num2cell(round(double(cell2mat(timestamp) - hdr.FirstTimeStamp)/hdr.TimeStampPerSample + 1));
      % convert it into a structure array
      event = struct('type', type, 'value', value, 'sample', sample, 'timestamp', timestamp, 'duration', duration, 'offset', offset, 'number', number);
    end
    
  case {'neuralynx_nev'}
    % instead of the directory containing the combination of nev and ncs files, the nev file was specified
    % do NOT read the header of the dataset
    
    % read the events, apply filter is applicable
    nev = read_neuralynx_nev(filename, 'type', flt_type, 'value', flt_value, 'mintimestamp', flt_mintimestamp, 'maxtimestamp', flt_maxtimestamp, 'minnumber', flt_minnumber, 'maxnumber', flt_maxnumber);
    
    % the following code should only be executed if there are events,
    % otherwise there will be an error subtracting an uint64 from an []
    if ~isempty(nev)
      % now get the values as cell array, since the struct function can work with those
      value     = {nev.TTLValue};
      timestamp = {nev.TimeStamp};
      number    = {nev.EventNumber};
      type      = repmat({'trigger'},size(value));
      duration  = repmat({[]},size(value));
      offset    = repmat({[]},size(value));
      % since the ncs files are not specified, there is no mapping between lfp samples and timestamps
      sample    = repmat({nan}, size(value));
      % convert it into a structure array
      event = struct('type', type, 'value', value, 'sample', sample, 'timestamp', timestamp, 'duration', duration, 'offset', offset, 'number', number);
    end
    
  case 'nmc_archive_k'
    event = read_nmc_archive_k_event(filename);
    
  case 'netmeg'
    warning('FieldTrip:ft_read_event:unsupported_event_format', 'reading of events for the netmeg format is not yet supported');
    event = [];
    
  case 'neuroshare' % NOTE: still under development
    % check that the required neuroshare toolbox is available
    ft_hastoolbox('neuroshare', 1);
    
    tmp = read_neuroshare(filename, 'readevent', 'yes');
    for i=1:length(tmp.event.timestamp)
      event(i).type      = tmp.hdr.eventinfo(i).EventType;
      event(i).value     = tmp.event.data(i);
      event(i).timestamp = tmp.event.timestamp(i);
      event(i).sample    = tmp.event.sample(i);
    end
    
  case 'neuralynx_cds'
    % this is a combined Neuralynx dataset with separate subdirectories for the LFP, MUA and spike channels
    dirlist   = dir(filename);
    %haslfp   = any(filetype_check_extension({dirlist.name}, 'lfp'));
    %hasmua   = any(filetype_check_extension({dirlist.name}, 'mua'));
    %hasspike = any(filetype_check_extension({dirlist.name}, 'spike'));
    %hastsl   = any(filetype_check_extension({dirlist.name}, 'tsl'));   % separate file with original TimeStampLow
    %hastsh   = any(filetype_check_extension({dirlist.name}, 'tsh'));   % separate file with original TimeStampHi
    hasttl    = any(filetype_check_extension({dirlist.name}, 'ttl'));   % separate file with original Parallel_in
    hasnev    = any(filetype_check_extension({dirlist.name}, 'nev'));   % original Events.Nev file
    hasmat    = 0;
    if hasttl
      eventfile = fullfile(filename, dirlist(find(filetype_check_extension({dirlist.name}, 'ttl'))).name);
      % read the header from the combined dataset
      if isempty(hdr)
        hdr = ft_read_header(filename);
      end
      % read the events from the *.ttl file
      event = ft_read_event(eventfile);
      % convert the sample numbers from the dma or ttl file to the downsampled dataset
      % assume that the *.ttl file is sampled at 32556Hz and is aligned with the rest of the data
      for i=1:length(event)
        event(i).sample = round((event(i).sample-1) * hdr.Fs/32556 + 1);
      end
      % elseif hasnev
      % FIXME, do something here
      % elseif hasmat
      % FIXME, do something here
    else
      error('no event file found');
    end
    
    %   The sample number is missingin the code below, since it is not available
    %   without looking in the continuously sampled data files. Therefore
    %   sorting the events (later in this function) based on the sample number
    %   fails and no events can be returned.
    %
    %   case 'neuralynx_nev'
    %     [nev] = read_neuralynx_nev(filename);
    %     % select only the events with a TTL value
    %     ttl = [nev.TTLValue];
    %     sel = find(ttl~=0);
    %     % now get the values as cell array, since teh struct function can work with those
    %     value     = {nev(sel).TTLValue};
    %     timestamp = {nev(sel).TimeStamp};
    %     event = struct('value', value, 'timestamp', timestamp);
    %     for i=1:length(event)
    %       % assign the other fixed elements
    %       event(i).type     = 'trigger';
    %       event(i).offset   = [];
    %       event(i).duration = [];
    %       event(i).sample   = [];
    %     end
    
  case {'neuroprax_eeg', 'neuroprax_mrk'}
    event = [];
    % start reading the markers, which I believe to be more like clinical annotations
    tmp = np_readmarker (filename, 0, inf, 'samples');
    for i = 1:numel(tmp.marker)
      if isempty(tmp.marker{i})
        break;
      end
      event = appendevent(event, struct('type', tmp.markernames(i), 'sample', num2cell(tmp.marker{i}), 'value', {tmp.markertyp(i)}));
    end
    % if present, read the digital triggers which are present as channel in the data
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    trgindx = match_str(hdr.label, 'DTRIG');
    if ~isempty(trgindx)
      trigger = read_trigger(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', trgindx, 'detectflank', detectflank, 'trigshift', trigshift);
      event   = appendevent(event, trigger);
    end
    
  case 'nexstim_nxe'
    event = read_nexstim_event(filename);
    
  case 'nihonkohden_m00'
    % in the data I tested the triggers are marked as DC offsets (deactivation of the DC channel)
    begsample = 1;
    if isempty(hdr)
      hdr = read_nihonkohden_hdr(filename);
    end
    
    if isfield(hdr, 'dat')
      % this is inefficient, since it keeps the complete data in memory
      % but it does speed up subsequent read operations without the user
      % having to care about it
      dat = hdr.dat;
    else
      dat = read_nihonkohden_m00(filename, begsample, hdr.nSamples);
    end
    
    % read the trigger channel and do flank detection
    event_chan = {'DC09','DC10','DC11','DC12'};
    
    trgindx = match_str(hdr.label, event_chan);
    trig = dat(trgindx,:);
    clear dat;
    begsample = 1;
    
    % marking offset trigger latencies
    %     onlat  = (diff([trig(:,1) trig],1,2)>0);
    offlat = (diff([trig trig(:,1)],1,2)<0);
    
    %     onset  = find(sum(double(onlat), 1)>0);
    offset = find(sum(double(offlat),1)>0);
    
    event = [];
    for j=1:size(offset,2);
      value = bin2dec(num2str(flipud(offlat(:,offset(j)))')); % flipup is needed to code bin2dec properly: DC09 = +1, DC10 = +2, DC11 = +4, DC12 = +8
      event(end+1).type   = 'down_flank';                     % distinguish between up and down flank
      event(end  ).sample = offset(j) + begsample - 1;        % assign the sample at which the trigger has gone down
      event(end  ).value  = value;                            % assign the trigger value just _before_ going down
    end
    
  case 'nimh_cortex'
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    cortex = hdr.orig.trial;
    for i=1:length(cortex)
      % add one 'trial' event for every trial and add the trigger events
      event(end+1).type     = 'trial';
      event(end  ).sample   = nan;
      event(end  ).duration = nan;
      event(end  ).offset   = nan;
      event(end  ).value    = i; % use the trial number as value
      for j=1:length(cortex(i).event)
        event(end+1).type     = 'trigger';
        event(end  ).sample   = nan;
        event(end  ).duration = nan;
        event(end  ).offset   = nan;
        event(end  ).value    = cortex(i).event(j);
      end
    end
    
  case 'ns_avg'
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    event(end+1).type     = 'average';
    event(end  ).sample   = 1;
    event(end  ).duration = hdr.nSamples;
    event(end  ).offset   = -hdr.nSamplesPre;
    event(end  ).value    = [];
    
  case {'ns_cnt', 'ns_cnt16', 'ns_cnt32'}
    % read the header, the original header includes the event table
    if isempty(hdr)
      hdr = ft_read_header(filename, 'headerformat', eventformat);
    end
    % translate the event table into known FieldTrip event types
    for i=1:numel(hdr.orig.event)
      event(end+1).type     = 'trigger';
      event(end  ).sample   = hdr.orig.event(i).offset + 1; % +1 was in EEGLAB pop_loadcnt
      event(end  ).value    = hdr.orig.event(i).stimtype;
      event(end  ).offset   = 0;
      event(end  ).duration = 0;
      
      % the code above assumes that all events are stimulus triggers
      % howevere, there are also interesting events possible, such as responses
      if hdr.orig.event(i).stimtype~=0
        event(end+1).type     = 'stimtype';
        event(end  ).sample   = hdr.orig.event(i).offset + 1; % +1 was in EEGLAB pop_loadcnt
        event(end  ).value    = hdr.orig.event(i).stimtype;
        event(end  ).offset   = 0;
        event(end ).duration  = 0;
      elseif hdr.orig.event(i).keypad_accept~=0
        event(end+1).type     = 'keypad_accept';
        event(end  ).sample   = hdr.orig.event(i).offset + 1; % +1 was in EEGLAB pop_loadcnt
        event(end  ).value    = hdr.orig.event(i).keypad_accept;
        event(end  ).offset   = 0;
        event(end  ).duration = 0;
      end
    end
    
  case 'ns_eeg'
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    for i=1:hdr.nTrials
      % the *.eeg file has a fixed trigger value for each trial
      % furthermore each trial has additional fields like accept, correct, response and rt
      
      tmp = read_ns_eeg(filename, i);
      % create an event with the trigger value
      event(end+1).type     = 'trial';
      event(end  ).sample   = (i-1)*hdr.nSamples + 1;
      event(end  ).value    = tmp.sweep.type;  % trigger value
      event(end  ).offset   = -hdr.nSamplesPre;
      event(end  ).duration =  hdr.nSamples;
      
      % create an event with the boolean accept/reject code
      event(end+1).type     = 'accept';
      event(end  ).sample   = (i-1)*hdr.nSamples + 1;
      event(end  ).value    = tmp.sweep.accept;  % boolean value indicating accept/reject
      event(end  ).offset   = -hdr.nSamplesPre;
      event(end  ).duration =  hdr.nSamples;
      
      % create an event with the boolean accept/reject code
      event(end+1).type     = 'correct';
      event(end  ).sample   = (i-1)*hdr.nSamples + 1;
      event(end  ).value    = tmp.sweep.correct;  % boolean value
      event(end  ).offset   = -hdr.nSamplesPre;
      event(end  ).duration =  hdr.nSamples;
      
      % create an event with the boolean accept/reject code
      event(end+1).type     = 'response';
      event(end  ).sample   = (i-1)*hdr.nSamples + 1;
      event(end  ).value    = tmp.sweep.response;  % probably a boolean value
      event(end  ).offset   = -hdr.nSamplesPre;
      event(end  ).duration =  hdr.nSamples;
      
      % create an event with the boolean accept/reject code
      event(end+1).type     = 'rt';
      event(end  ).sample   = (i-1)*hdr.nSamples + 1;
      event(end  ).value    = tmp.sweep.rt;  % time in seconds
      event(end  ).offset   = -hdr.nSamplesPre;
      event(end  ).duration =  hdr.nSamples;
    end
    
  case 'plexon_nex'
    event = read_nex_event(filename);
    
  case 'tmsi_poly5'
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    trgindx = find(strcmp(hdr.chantype, 'trigger'));
    if ~isempty(trgindx)
      % the "Digi" value goes down from 255 to 254
      detectflank = 'downdiff';
      event = read_trigger(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', trgindx, 'detectflank', detectflank, 'trigshift', trigshift);
    end
    
  case {'yokogawa_ave', 'yokogawa_con', 'yokogawa_raw'}
    % check that the required low-level toolbox is available
    if ~ft_hastoolbox('yokogawa', 0);
      ft_hastoolbox('yokogawa_meg_reader', 1);
    end
    % the user should be able to specify the analog threshold
    % the user should be able to specify the trigger channels
    % the user should be able to specify the flank, but the code falls back to 'auto' as default
    if isempty(detectflank)
      detectflank = 'auto';
    end
    event = read_yokogawa_event(filename, 'detectflank', detectflank, 'trigindx', trigindx, 'threshold', threshold);
    
  case 'oxy3'
    ft_hastoolbox('artinis', 1);
    
    event = read_artinis_oxy3(filename, true);
    if isempty(hdr)
      hdr = read_artinis_oxy3(filename);
    end
    
    if isempty(trigindx) % indx gets precedence over labels! numbers before words
      triglabel = ft_getopt(varargin, 'triglabel', 'ADC*');  % this allows subselection of AD channels to be markes as trigger channels (for Artinis oxy3 data)
      trigindx = find(ismember(hdr.label, ft_channelselection(triglabel, hdr.label)));
    end
    
    % read the trigger channel and do flank detection
    trigger = read_trigger(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'threshold', threshold, 'chanindx', trigindx, 'detectflank', detectflank, 'trigshift', trigshift, 'fixartinis', true);
    
    % remove consecutive triggers
    if ~isempty(trigger)
      i = 1;
      last_trigger_sample = trigger(i).sample;
      while i<numel(trigger)
        if strcmp(trigger(i).type, trigger(i+1).type) && trigger(i+1).sample-last_trigger_sample <= tolerance
          [trigger(i).value, idx] = max([trigger(i).value, trigger(i+1).value]);
          fprintf('Merging triggers at sample %d and %d\n', trigger(i).sample, trigger(i+1).sample);
          last_trigger_sample =  trigger(i+1).sample;
          if (idx==2)
            trigger(i).sample = trigger(i+1).sample;
          end
          
          trigger(i+1) = [];
        else
          i=i+1;
          last_trigger_sample = trigger(i).sample;
        end
      end
      
      event = appendevent(event, trigger);
    end
    
  case 'spmeeg_mat'
    if isempty(hdr)
      hdr = ft_read_header(filename);
    end
    event = read_spmeeg_event(filename, 'header', hdr);
    
  otherwise
    % attempt to run eventformat as a function
    % in case using an external read function was desired, this is where it is executed
    % if it fails, the regular unsupported warning message is thrown
    try
      event = feval(eventformat,filename);
    catch
      warning('FieldTrip:ft_read_event:unsupported_event_format','unsupported event format (%s)', eventformat);
      event = [];
    end
end

if ~isempty(hdr) && hdr.nTrials>1 && (isempty(event) || ~any(strcmp({event.type}, 'trial')))
  % the data suggests multiple trials and trial events have not yet been defined
  % make an event for each trial according to the file header
  for i=1:hdr.nTrials
    event(end+1).type     = 'trial';
    event(end  ).sample   = (i-1)*hdr.nSamples + 1;
    event(end  ).offset   = -hdr.nSamplesPre;
    event(end  ).duration =  hdr.nSamples;
    event(end  ).value    =  [];
  end
end

if ~isempty(event)
  % make sure that all required elements are present
  if ~isfield(event, 'type'),     error('type field not defined for each event');     end
  if ~isfield(event, 'sample'),   error('sample field not defined for each event');   end
  if ~isfield(event, 'value'),    for i=1:length(event), event(i).value = [];    end; end
  if ~isfield(event, 'offset'),   for i=1:length(event), event(i).offset = [];   end; end
  if ~isfield(event, 'duration'), for i=1:length(event), event(i).duration = []; end; end
end

% make sure that all numeric values are double
if ~isempty(event)
  for i=1:length(event)
    if isnumeric(event(i).value)
      event(i).value = double(event(i).value);
    end
    event(i).sample    = double(event(i).sample);
    event(i).offset    = double(event(i).offset);
    event(i).duration  = double(event(i).duration);
  end
end

if ~isempty(event)
  % sort the events on the sample on which they occur
  % this has the side effect that events without a sample number are discarded
  sample = [event.sample];
  if ~all(isnan(sample))
    [dum, indx] = sort(sample);
    event = event(indx);
  end
end

% apply the optional filters
event = ft_filter_event(event, varargin{:});

if isempty(event)
  % ensure that it has the correct fields, even if it is empty
  event = struct('type', {}, 'value', {}, 'sample', {}, 'offset', {}, 'duration', {});
end
