function event = events_tsv(filename)

% helper function to extract events from a BIDS events.tsv file

opts = detectImportOptions(filename,'filetype','text');

% this gets the numeric stuff as doubles right away, yet might erroneously
% typecast intended strings (e.g. in event_value) into doubles.
tsv       = readtable(filename, opts); 
tsv_naive = readtable(filename, 'filetype', 'text');

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
    try,
      for k = 1:numel(event)
        event(k).type = tsv(k).event_type;
      end
    end
  end
  if ~isfield(event, 'value')
    try,
      for k = 1:numel(event)
        event(k).value = tsv(k).event_value;
      end
    end
  end
  if ~isfield(event, 'offset')
    event(1).offset = [];
  end
  if ~isfield(event, 'sample')
    event(1).sample = [];
  end
  if ~isfield(event, 'duration')
    event(1).duration = [];
  end
end


function bool = is_a_number(input)

pattern = '[0-9\.]';
bool    = length(regexp(input, pattern))==length(input);

function bool = is_a_nan(input)

pattern = 'nan';
bool    = ~isempty(regexp(lower(input), pattern));

