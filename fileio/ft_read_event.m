function [event] = read_event(filename, varargin)

% READ_EVENT reads all events from an EEG/MEG dataset and returns
% them in a well defined structure. It is a wrapper around different
% EEG/MEG file importers, directly supported formats are CTF, Neuromag,
% EEP, BrainVision, Neuroscan and Neuralynx.
%
% Use as
%   [event] = read_event(filename, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'eventformat'   string
%   'header'        structure, see READ_HEADER
%   'detectflank'   string, can be 'up', 'down' or 'both' (default = 'up')
%   'trigshift'     integer, number of samples to shift from flank to detect trigger value (default = 0)
%
% Furthermore, you can specify optional arguments as key-value pairs
% for filtering the events, e.g. to select only events of a specific
% type, of a specific value, or events between a specific begin and
% end sample. This event filtering is especially usefull for real-time
% processing. See FILTER_EVENT for more details.
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
% See also READ_HEADER, READ_DATA, WRITE_DATA, WRITE_EVENT, FILTER_EVENT

% Copyright (C) 2004-2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

persistent sock       % for fcdc_tcp

global event_queue        % for fcdc_global
persistent db_blob        % for fcdc_mysql
if isempty(db_blob)
  db_blob = 0;
end

if iscell(filename)
  % use recursion to read from multiple event sources
  event = [];
  for i=1:numel(filename)
    tmp   = read_event(filename{i}, varargin{:});
    event = appendevent(event(:), tmp(:));
  end
  return
end

% get the options
eventformat      = keyval('eventformat',  varargin);
hdr              = keyval('header',       varargin);
detectflank      = keyval('detectflank',  varargin); % up, down or both
trigshift        = keyval('trigshift',    varargin); % default is assigned in subfunction
trigindx         = keyval('trigindx',     varargin); % default is based on chantype helper function
headerformat     = keyval('headerformat', varargin);
dataformat       = keyval('dataformat',   varargin);

% this allows to read only events in a certain range, supported for selected data formats only
flt_type         = keyval('type',         varargin);
flt_value        = keyval('value',        varargin);
flt_minsample    = keyval('minsample',    varargin);
flt_maxsample    = keyval('maxsample',    varargin);
flt_mintimestamp = keyval('mintimestamp', varargin);
flt_maxtimestamp = keyval('maxtimestamp', varargin);
flt_minnumber    = keyval('minnumber', varargin);
flt_maxnumber    = keyval('maxnumber', varargin);


% determine the filetype
if isempty(eventformat)
  eventformat = filetype(filename);
end

% default is to search only for rising or up-going flanks
if isempty(detectflank)
  detectflank = 'up';
end

switch eventformat
  case 'brainvision_vhdr'
    % read the headerfile belonging to the dataset and try to determine the corresponding markerfile
    eventformat = 'brainvision_vmrk';
    hdr = read_brainvision_vhdr(filename);
    % replace the filename with the filename of the markerfile
    if ~isfield(hdr, 'MarkerFile') || isempty(hdr.MarkerFile)
      filename = [];
    else
      [p, f, e] = fileparts(filename);
      filename = fullfile(p, hdr.MarkerFile);
    end
end

% start with an empty event structure
event = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the events with the low-level reading function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch eventformat

  case 'fcdc_global'
    event = event_queue;

  case {'4d' '4d_pdf', '4d_m4d', '4d_xyz'}
    if isempty(hdr)
      hdr = read_header(filename);
    end
    % add the trials to the event structure
    for i=1:hdr.nTrials
      event(end+1).type     = 'trial';
      event(end  ).sample   = (i-1)*hdr.nSamples + 1;
      event(end  ).value    = [];
      event(end  ).offset   = -hdr.nSamplesPre;
      event(end  ).duration = hdr.nSamples;
    end
    % read the trigger channel and do flank detection
    trgindx = match_str(hdr.label, 'TRIGGER');
    if isfield(hdr, 'orig') && isfield(hdr.orig, 'config_data') && strcmp(hdr.orig.config_data.site_name, 'Glasgow'),
      trigger = read_trigger(filename, 'header', hdr, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', trgindx, 'detectflank', detectflank, 'trigshift', trigshift,'fix4dglasgow',1);
    else
      trigger = read_trigger(filename, 'header', hdr, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', trgindx, 'detectflank', detectflank, 'trigshift', trigshift,'fix4dglasgow',0);
    end
    event   = appendevent(event, trigger);

    respindx = match_str(hdr.label, 'RESPONSE');
    if ~isempty(respindx)
      response = read_trigger(filename, 'header', hdr, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', respindx, 'detectflank', detectflank, 'trigshift', trigshift);
      event    = appendevent(event, response);
    end

  case 'bci2000_dat'
    % this requires the load_bcidat mex file to be present on the path
    hastoolbox('BCI2000', 1);
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

  case {'besa_avr', 'besa_swf'}
    if isempty(hdr)
      hdr = read_header(filename);
    end
    event(end+1).type     = 'average';
    event(end  ).sample   = 1;
    event(end  ).duration = hdr.nSamples;
    event(end  ).offset   = -hdr.nSamplesPre;
    event(end  ).value    = [];

  case {'biosemi_bdf', 'bham_bdf'}
    % read the header, required to determine the stimulus channels and trial specification
    if isempty(hdr)
      hdr = read_header(filename);
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
    sdata = read_data(filename, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', schan);

    % find indices of negative numbers
    bit24i = find(sdata < 0);
    % make number positive and preserve bits 0-22
    sdata(bit24i) = bitcmp(abs(sdata(bit24i))-1,24);
    % re-insert the sign bit on its original location, i.e. bit24
    sdata(bit24i) = sdata(bit24i)+(2^(24-1));
    % typecast the data to ensure that the status channel is represented in 32 bits
    sdata = uint32(sdata);

    byte1 = 2^8  - 1;
    byte2 = 2^16 - 1 - byte1;
    byte3 = 2^24 - 1 - byte1 - byte2;

    % get the respective status and trigger bits
    trigger = bitand(sdata, bitor(byte1, byte2)); %  contained in the lower two bytes
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
          tok = tokenize(line, '=', 0);    % do not squeeze repetitions of the seperator
          if length(tok)~=2
            warning('skipping unexpected formatted line in BrainVision marker file');
          else
            % the line looks like "MkXXX=YYY", which is ok
            % the interesting part now is in the YYY, i.e. the second token
            tok = tokenize(tok{2}, ',', 0);    % do not squeeze repetitions of the seperator
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

  case 'ced_son'
    % check that the required low-level toolbox is available
    hastoolbox('neuroshare', 1);
    orig = read_ced_son(filename,'readevents','yes');
    event = struct('type',     {orig.events.type},...
      'sample',   {orig.events.sample},...
      'value',    {orig.events.value},...
      'offset',   {orig.events.offset},...
      'duration', {orig.events.duration});

  case {'ctf_ds', 'ctf_meg4', 'ctf_res4', 'ctf_old'}
    % obtain the dataset name
    if filetype(filename, 'ctf_meg4') ||  filetype(filename, 'ctf_res4')
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
      hdr = read_header(headerfile, 'headerformat', headerformat);
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
    % meg channels are 5, refmag 0, refgrad 1, adcs 18, trigger 11, eeg 9
    trigchanindx = find(origSensType==11);
    if ~isempty(trigchanindx)
      % read the trigger channel and do flank detection
      trigger = read_trigger(filename, 'header', hdr, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', trigchanindx, 'dataformat', dataformat, 'detectflank', detectflank, 'trigshift', trigshift, 'fixctf', 1);
      event   = appendevent(event, trigger);
    end

    % make an event for each trial as defined in the header
    for i=1:hdr.nTrials
      event(end+1).type     = 'trial';
      event(end  ).sample   = (i-1)*hdr.nSamples + 1;
      event(end  ).offset   = -hdr.nSamplesPre;
      event(end  ).duration =  hdr.nSamples;
      event(end  ).value    =  [];
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

  case 'eeglab_set'
    event = read_eeglabevent(filename, 'header', hdr);

  case 'spmeeg_mat'
    event = read_spmeeg_event(filename, 'header', hdr);

  case 'eep_avr'
    % check that the required low-level toolbox is available
    hastoolbox('eeprobe', 1);
    % the headerfile and datafile are the same
    if isempty(hdr)
      hdr = read_header(filename);
    end
    event(end+1).type     = 'average';
    event(end  ).sample   = 1;
    event(end  ).duration = hdr.nSamples;
    event(end  ).offset   = -hdr.nSamplesPre;
    event(end  ).value    = [];

  case 'eep_cnt'
    % check that the required low-level toolbox is available
    hastoolbox('eeprobe', 1);
    % try to read external trigger file in EEP format
    trgfile = [filename(1:(end-3)), 'trg'];
    if exist(trgfile, 'file')
      if isempty(hdr)
        hdr = read_header(filename);
      end
      tmp = read_eep_trg(trgfile);
      % translate the EEProbe trigger codes to events
      for i=1:length(tmp)
        event(i).type     = 'trigger';
        event(i).sample   = round((tmp(i).time/1000) * hdr.Fs) + 1;    % convert from ms to samples
        event(i).value    = tmp(i).code;
        event(i).offset   = 0;
        event(i).duration = 0;
      end
    else
      warning('no triggerfile was found');
    end

  case 'egi_egis'
    if isempty(hdr)
      hdr = read_header(filename);
    end
    fhdr   = hdr.orig.fhdr;
    chdr   = hdr.orig.chdr;
    ename  = hdr.orig.ename;
    cnames = hdr.orig.cnames;
    fcom   = hdr.orig.fcom;
    ftext  = hdr.orig.ftext;
    eventCount=0;
    for cell=1:fhdr(18)
      for trial=1:chdr(cell,2)
        eventCount=eventCount+1;
        event(eventCount).type     = 'trial';
        event(eventCount).sample   = (eventCount-1)*hdr.nSamples + 1;
        event(eventCount).offset   = -hdr.nSamplesPre;
        event(eventCount).duration =  hdr.nSamples;
        event(eventCount).value    =  cnames{cell};
      end
    end

  case 'egi_egia'
    if isempty(hdr)
      hdr = read_header(filename);
    end
    fhdr   = hdr.orig.fhdr;
    chdr   = hdr.orig.chdr;
    ename  = hdr.orig.ename;
    cnames = hdr.orig.cnames;
    fcom   = hdr.orig.fcom;
    ftext  = hdr.orig.ftext;
    eventCount=0;
    for cell=1:fhdr(18)
      for subject=1:chdr(cell,2)
        eventCount=eventCount+1;
        event(eventCount).type     = 'trial';
        event(eventCount).sample   = (eventCount-1)*hdr.nSamples + 1;
        event(eventCount).offset   = -hdr.nSamplesPre;
        event(eventCount).duration =  hdr.nSamples;
        event(eventCount).value    =  ['Sub' sprintf('%03d',subject) cnames{cell}];
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
      hdr = read_header(filename,'headerformat','egi_sbin');
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
                    if event(eventCount).duration == 0
                        event(eventCount).type     = 'trigger';
                    else
                        event(eventCount).type     = 'trial';
                    end;
                    event(eventCount).value    =  char(EventCodes(theEvent,:));
                end
            end
        end
    end
    
    for segment=1:hdr.nTrials  % cell information
      eventCount=eventCount+1;
      event(eventCount).type     = 'trial';
      event(eventCount).sample   = (segment-1)*hdr.nSamples + 1;
      event(eventCount).offset   = -hdr.nSamplesPre;
      event(eventCount).duration =  hdr.nSamples;
      if unsegmented,
          event(eventCount).value    = [];
      else
          event(eventCount).value    =  char([CateNames{segHdr(segment,1)}(1:CatLengths(segHdr(segment,1)))]);
      end
    end

  case 'fcdc_buffer'
    % read from a networked buffer for realtime analysis
    [host, port] = filetype_check_uri(filename);

    % FIXME it should be possible to specify event numbers

    try
      evt = buffer('get_evt', [], host, port);  % indices should be zero-offset
    catch
      if strfind(lasterr, 'the buffer returned an error')
        evt = [];
      else
        rethrow(lasterr);
      end
    end

    type = {
      'char'
      'uint8'
      'uint16'
      'uint32'
      'uint64'
      'int8'
      'int16'
      'int32'
      'int64'
      'single'
      'double'
      };

    wordsize = {
      1 % 'char'
      1 % 'uint8'
      2 % 'uint16'
      4 % 'uint32'
      8 % 'uint64'
      1 % 'int8'
      2 % 'int16'
      4 % 'int32'
      8 % 'int64'
      4 % 'single'
      8 % 'double'
      };

    for i=1:length(evt)
      % convert the field "type" into the Matlab representation
      this_type = type{evt(i).type_type+1};
      this_size = wordsize{evt(i).type_type+1} * evt(i).type_numel;
      sel = 1:this_size;
      if strcmp(this_type, 'char')
        event(i).type = char(evt(i).buf(sel));
      else
        event(i).type = typecast(evt(i).buf(sel), this_type);
      end

      % convert the field "value" into the Matlab representation
      this_type = type{evt(i).value_type+1};
      this_size = wordsize{evt(i).value_type+1} * evt(i).value_numel;
      sel = sel(end) + (1:this_size);
      if strcmp(this_type, 'char')
        event(i).value = char(evt(i).buf(sel));
      else
        event(i).value = typecast(evt(i).buf(sel), this_type);
      end

      % the other fields are simple, because they have a fixed type and only a single elements
      event(i).sample   = evt(i).sample;
      event(i).offset   = evt(i).offset;
      event(i).duration = evt(i).duration;
    end

  case 'fcdc_matbin'
    % this is multiplexed data in a *.bin file, accompanied by a matlab file containing the header and event
    [path, file, ext] = fileparts(filename);
    filename = fullfile(path, [file '.mat']);
    % read the events from the Matlab file
    tmp   = load(filename, 'event');
    event = tmp.event;


  case 'fcdc_fifo'

    
    fifo = filetype_check_uri(filename);
    
    if ~exist(fifo,'file')
      warning('the FIFO %s does not exist; attempting to create it', fifo);
      system(sprintf('mkfifo -m 0666 %s',fifo));
    end
    
    fid = fopen(fifo, 'r');
    msg = fread(fid, inf, 'uint8');
    fclose(fid);
    
    try
      event = mxDeserialize(uint8(msg));
    catch
      warning(lasterr);
    end    

  case 'fcdc_tcp'
    % requires tcp/udp/ip-toolbox
    hastoolbox('TCP_UDP_IP', 1);
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
    hastoolbox('TCP_UDP_IP', 1);
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
    % this code is moved to a seperate file
    event = read_serial_event(filename);

  case 'fcdc_mysql'
    % read from a MySQL server listening somewhere else on the network
    db_open(filename);
    if db_blob
      event = db_select_blob('fieldtrip.event', 'msg');
    else
      event = db_select('fieldtrip.event', {'type', 'value', 'sample', 'offset', 'duration'});
    end

  case 'itab_raw'
    error('suppoport for events in this fileformat is not yet implemented')

  case 'matlab'
    % read the events from a normal Matlab file
    tmp   = load(filename, 'event');
    event = tmp.event;

  case {'mpi_ds', 'mpi_dap'}
    if isempty(hdr)
      hdr = read_header(filename);
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


  case {'neuromag_fif' 'neuromag_mne' 'neuromag_mex'}
    if strcmp(eventformat, 'neuromag_fif')
      % the default is to use the MNE reader for fif files
      eventformat = 'neuromag_mne';
    end
    if strcmp(eventformat, 'neuromag_mex')
      % check that the required low-level toolbox is available
      hastoolbox('meg-pd', 1);
      if isempty(headerformat), headerformat = eventformat; end
      if isempty(dataformat),   dataformat   = eventformat; end
    elseif strcmp(eventformat, 'neuromag_mne')
      % check that the required low-level toolbox is available
      hastoolbox('mne', 1);
      if isempty(headerformat), headerformat = eventformat; end
      if isempty(dataformat),   dataformat   = eventformat; end
    end

    if isempty(hdr)
      hdr = read_header(filename, 'headerformat', headerformat);
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
      analogindx = find(strcmp(chantype(hdr), 'analog trigger'));
      binaryindx = find(strcmp(chantype(hdr), 'digital trigger'));
      
      
      if isempty(binaryindx)&&isempty(analogindx) 
        % included in case of problems with older systems and MNE reader:
        % use a predefined set of channel names
        binary     = {'STI 014', 'STI 015', 'STI 016'};
        binaryindx = match_str(hdr.label, binary);
      end
      
      if ~isempty(binaryindx)
        trigger = read_trigger(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', binaryindx, 'detectflank', detectflank, 'trigshift', trigshift, 'fixneuromag', 0);
        event   = appendevent(event, trigger);
      end
      if ~isempty(analogindx)
        % add the triggers to the event structure based on trigger channels with the name "STI xxx"
        % there are some issues with noise on these analog trigger
        % channels, on older systems only
        % read the trigger channel and do flank detection
        trigger = read_trigger(filename, 'header', hdr, 'dataformat', dataformat, 'begsample', flt_minsample, 'endsample', flt_maxsample, 'chanindx', analogindx, 'detectflank', detectflank, 'trigshift', trigshift, 'fixneuromag', 1);
        event   = appendevent(event, trigger);
      end

      if hdr.nTrials>1
        % make an event for each trial as defined in the header
        for i=1:hdr.nTrials
          event(end+1).type     = 'trial';
          event(end  ).sample   = (i-1)*hdr.nSamples + 1;
          event(end  ).offset   = -hdr.nSamplesPre;
          event(end  ).duration =  hdr.nSamples;
          event(end  ).value    =  [];
        end
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


  case {'neuralynx_ttl' 'neuralynx_bin' 'neuralynx_dma' 'neuralynx_sdma'}
    if isempty(hdr)
      hdr = read_header(filename);
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
      % determine the seperate files with the trigger and timestamp information
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
      % read the trigger values from the seperate file
      ttl = read_neuralynx_bin(ttlfile, begsample, endsample);
    elseif strcmp(eventformat, 'neuralynx_ttl')
      % determine the optional files with timestamp information
      tslfile = [filename(1:(end-4)) '.tsl'];
      tshfile = [filename(1:(end-4)) '.tsh'];
      % read the triggers from a seperate *.ttl file
      ttl = read_neuralynx_ttl(filename, begsample, endsample);
    elseif strcmp(eventformat, 'neuralynx_bin')
      % determine the optional files with timestamp information
      tslfile = [filename(1:(end-8)) '.tsl.bin'];
      tshfile = [filename(1:(end-8)) '.tsh.bin'];
      % read the triggers from a seperate *.ttl.bin file
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
      hdr = read_header(filename);
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

  case 'neuralynx_cds'
    % this is a combined Neuralynx dataset with seperate subdirectories for the LFP, MUA and spike channels
    dirlist   = dir(filename);
    %haslfp   = any(filetype_check_extension({dirlist.name}, 'lfp'));
    %hasmua   = any(filetype_check_extension({dirlist.name}, 'mua'));
    %hasspike = any(filetype_check_extension({dirlist.name}, 'spike'));
    %hastsl   = any(filetype_check_extension({dirlist.name}, 'tsl'));   % seperate file with original TimeStampLow
    %hastsh   = any(filetype_check_extension({dirlist.name}, 'tsh'));   % seperate file with original TimeStampHi
    hasttl    = any(filetype_check_extension({dirlist.name}, 'ttl'));   % seperate file with original Parallel_in
    hasnev    = any(filetype_check_extension({dirlist.name}, 'nev'));   % original Events.Nev file
    hasmat    = 0;
    if hasttl
      eventfile = fullfile(filename, dirlist(find(filetype_check_extension({dirlist.name}, 'ttl'))).name);
      % read the header from the combined dataset
      if isempty(hdr)
        hdr = read_header(filename);
      end
      % read the events from the *.ttl file
      event = read_event(eventfile);
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
    tmp = np_readmarker (filename, 0, inf, 'samples');    
    event = [];
    for i = 1:numel(tmp.marker)
        if isempty(tmp.marker{i})
            break;
        end
        event = [event struct('type', tmp.markernames(i),...
            'sample', num2cell(tmp.marker{i}),...
            'value', {tmp.markertyp(i)})];
    end
        
  case 'nexstim_nxe'
    event = read_nexstim_event(filename);

  case 'nimh_cortex'
    if isempty(hdr)
      hdr = read_header(filename);
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
      hdr = read_header(filename);
    end
    event(end+1).type     = 'average';
    event(end  ).sample   = 1;
    event(end  ).duration = hdr.nSamples;
    event(end  ).offset   = -hdr.nSamplesPre;
    event(end  ).value    = [];

  case {'ns_cnt', 'ns_cnt16', 'ns_cnt32'}
    % read the header, the original header includes the event table
    if isempty(hdr)
      hdr = read_header(filename, 'headerformat', eventformat);
    end
    % translate the event table into known FieldTrip event types
    for i=1:numel(hdr.orig.event)
      event(i).type     = 'trigger';
      event(i).sample   = hdr.orig.event(i).offset + 1; % +1 was in EEGLAB pop_loadcnt
      event(i).value    = hdr.orig.event(i).stimtype;
      event(i).offset   = 0;
      event(i).duration = 0;
    end

  case 'ns_eeg'
    if isempty(hdr)
      hdr = read_header(filename);
    end
    for i=1:hdr.nTrials
      % the *.eeg file has a fixed trigger value for each trial
      % furthermore each trial has the label 'accept' or 'reject'
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
    end

  case 'plexon_nex'
    event = read_nex_event(filename);

  case {'yokogawa_ave', 'yokogawa_con', 'yokogawa_raw'}
    % check that the required low-level toolbox is available
    hastoolbox('yokogawa', 1);
    % allow the user to specify custom trigger channels
    event = read_yokogawa_event(filename, 'trigindx', trigindx);
    
  case 'nmc_archive_k'
    event = read_nmc_archive_k_event(filename);

  case 'neuroshare' % NOTE: still under development
    % check that the required neuroshare toolbox is available
    hastoolbox('neuroshare', 1);
    
    tmp = read_neuroshare(filename, 'readevent', 'yes');
    for i=1:length(tmp.hdr.eventinfo)
      event(i).type      = tmp.hdr.eventinfo(i).EventType;
      event(i).value     = tmp.event.data(i);
      event(i).timestamp = tmp.event.timestamp(i);
      event(i).sample    = tmp.event.sample(i);
    end
    
  otherwise
    error('unsupported event format');
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
  [dum, indx] = sort([event.sample]);
  event = event(indx);
% else
%   warning(sprintf('no events found in %s', filename));
end

% apply the optional filters
event = filter_event(event, varargin{:});


