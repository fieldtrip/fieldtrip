function result = events_tsv(filename, hdr)

% BIDS_TSV reads events from a BIDS _events.tsv file
%
% Use as
%   hdr = bids_tsv(filename)
% to read the header information, or use as
%   evt = bids_tsv(filename, hdr)
% to read the event information.
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT,
% QUALISYS_TSV, MOTION_C3D

opts = detectImportOptions(filename,'filetype','text', 'Delimiter', {'\t'} );

needhdr = nargin==1;
needevt = nargin==2;

if needhdr
  megjson = [filename(1:(end-length('events.tsv'))) 'meg.json'];
  eegjson = [filename(1:(end-length('events.tsv'))) 'eeg.json'];
  ieegjson = [filename(1:(end-length('events.tsv'))) 'ieeg.json'];
  nirsjson = [filename(1:(end-length('events.tsv'))) 'nirs.json'];
  if exist(megjson, 'file')
    jsonfile = megjson;
  elseif exist(eegjson, 'file')
    jsonfile = eegjson;
  elseif exist(ieegjson, 'file')
    jsonfile = ieegjson;
  elseif exist(nirsjson, 'file')
      jsonfile = nirsjson;
  else
    ft_error('cannot find header information correspoding to %s', filename);
  end
  ft_info('reading header information from %s', jsonfile);
  
  orig = read_json(jsonfile);
  
  hdr.Fs           = orig.SamplingFrequency;
  hdr.nSamples     = round(orig.RecordingDuration*orig.SamplingFrequency);
  hdr.nSamplesPre  = 0;   % assume continuous
  hdr.nTrials      = 1;   % assume continuous
  hdr.nChans       = 0;   % number of channels
  hdr.label        = {};  % Nx1 cell-array with the label of each channel
  hdr.chantype     = {};  % Nx1 cell-array with the channel type, see FT_CHANTYPE
  hdr.chanunit     = {};  % Nx1 cell-array with the physical units, see FT_CHANUNIT
  
  fn = fieldnames(orig);
  sel = find(endsWith(fn, 'ChannelCount'));
  for i=1:numel(sel)
    hdr.nChans = hdr.nChans + orig.(fn{sel(i)});
  end
  for i=1:hdr.nChans
    hdr.label{i} = num2str(i);
    hdr.chantype{i} = 'unknown';
    hdr.chanunit{i} = 'unknown';
  end
  
  % keep the original header details
  hdr.orig = orig;
  
  % return the header
  result = hdr;
end % if needhdr

if needevt
  % this gets the numeric stuff as doubles right away, yet might erroneously
  % typecast intended strings (e.g. in event_value) into doubles.
  tsv       = readtable(filename, opts);
  tsv_naive = readtable(filename, 'filetype', 'text', 'Delimiter', {'\t'});
  
  tsv       = table2struct(tsv);
  tsv_naive = table2struct(tsv_naive);
  
  fn = fieldnames(tsv);
  checkflag = false(1,numel(fn));
  for k = 1:numel(fn)
    checkflag(1,k) = ~ischar(tsv(1).(fn{k})) && ischar(tsv_naive(1).(fn{k}));
  end
  
  % heuristic: the element can stay double if all elements of the corresponding
  % element in tsv_naive contain only numbers, or nan
  for k = find(checkflag)
    for m = 1:numel(tsv_naive)
      if ~is_a_number(tsv_naive(m).(fn{k})) && ~is_a_nan(tsv_naive(m).(fn{k}))
        tsv(m).(fn{k}) = tsv_naive(m).(fn{k});
      end
    end
  end
  
  % a FieldTrip event structure should contain
  % type, value, sample, offset and duration as fields
  %
  % a BIDS event table contains
  % event_type, event_value
  
  event = tsv;
  %event = struct([]);
  %event = copyfields(tsv, event, {'type' 'value' 'sample' 'offset' 'duration'});
  if ~isempty(event)
    if ~isfield(event, 'type')
      try
        for k = 1:numel(event)
          event(k).type = tsv(k).event_type;
        end
      end
    end
    if ~isfield(event, 'value')
      try
        for k = 1:numel(event)
          event(k).value = tsv(k).event_value;
        end
      end
    end
    
    if ~isfield(event, 'sample')
      if nargin>2
        % construct the sample number from the onset (in seconds)
        for i=1:numel(event)
          event(i).sample1 = round(event(i).onset*hdr.Fs) + 1;
        end
      else
        % we don't know the sampling rate, so cannot determine the sample number
        event(1).sample = [];
      end
    end
    
    if ~isfield(event, 'offset')
      event(1).offset = [];
    end
    
    if ~isfield(event, 'duration')
      event(1).duration = [];
    end
  end
  
  % return the events
  result = event;
end % if needevt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function json = read_json(filename)
ft_info('reading %s\n', filename);
ft_hastoolbox('jsonlab', 1);
json = loadjson(filename);
json = ft_struct2char(json); % convert strings into char-arrays

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool = is_a_number(input)
pattern = '[0-9\.]';
bool    = length(regexp(input, pattern))==length(input);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool = is_a_nan(input)
pattern = 'nan';
bool    = ~isempty(regexp(lower(input), pattern));
