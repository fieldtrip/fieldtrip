function event = events_tsv(filename)

% helper function to extract events from a BIDS events.tsv file

opts = detectImportOptions(filename,'filetype','text');
tsv  = readtable(filename, opts); % this gets the numeric stuff as doubles right away
tsv  = table2struct(tsv);

% a FieldTrip event structure should contain
% type, value, sample, offset and duration as fields
%
% a BIDS event table contains
% event_type, event_value

event = struct([]);
event = copyfields(tsv, event, {'type' 'value' 'sample' 'offset' 'duration'});
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
