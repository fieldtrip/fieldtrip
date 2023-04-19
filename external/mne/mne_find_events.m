function [eventlist] = mne_find_events(fname, stim_channel, consecutive, output)
%
%
%   [eventlist] = mne_find_events(fname, stim_channel, consecutive, output)
%
%   Find event from raw file
%
%   fname        - string; .fiff raw data file name
%   stim_channel - int; the channel that record event
%   consecutive  - bool | 'increasing'
%                  If True, consider instances where the value of the events
%                  channel changes without first returning to zero as multiple
%                  events. If False, report only instances where the value of the
%                  events channel changes from/to zero. If 'increasing', report
%                  adjacent events only when the second event code is greater than
%                  the first.
%   output       - 'onset' | 'offset' | 'step'
%                  Whether to report when events start, when events end, or both.
%
%   eventlist    - size = (n_events, 3)
%                  The first column contains the event time in samples and the third
%                  column contains the event id. If output = 'onset' or 'step', the
%                  second column contains the value of the stim channel immediately
%                  before the event/step. For output = 'offset', the second column
%                  contains the value of the stim channel after the event offset.
%
%   Authors: Fu-Te Wong (zuxfoucault@gmail.com),
%            Chien-Chung Chen / Visual Neuroscience Lab, National Taiwan University
%   Version 1.0 2017/9/17
%   License: BSD (3-clause)

raw = fiff_setup_read_raw(fname);
[data, times] = fiff_read_raw_segment(raw);
pick = stim_channel; % stim_channel

changed = diff(data(pick, :));
changed_idx = find(changed ~= 0);
if length(changed_idx) == 0,
    eventlist = zeros(1,3);
    return
end

pre_step = data(pick, changed_idx);
changed_idx = changed_idx + 1;
post_step = data(pick, changed_idx);
changed_idx = changed_idx + double(raw.first_samp - 1);
eventlist = cat(2, changed_idx', pre_step', post_step');


if strcmpi(consecutive, 'increasing'),
    onsets = eventlist(:, 3) > eventlist(:, 2);
    offsets = (onsets | eventlist(:, 3) == 0) & eventlist(:, 2) > 0;
end

if strcmpi(consecutive, 'True'),
    onsets = eventlist(:, 3) > 0;
    offsets = eventlist(:, 2) > 0;
end

if strcmpi(consecutive, 'False'),
    onsets = eventlist(:, 2) == 0;
    offsets = eventlist(:, 3) == 0;
end

onset_idx = find(onsets);
offset_idx = find(offsets);

if length(onset_idx) == 0 | length(offset_idx) == 0,
    eventlist = zeros(1,3);
    return
end

% Delete orphaned onsets/offsets
if onset_idx(1) > offset_idx(1),
    disp('Removing orphaned offset at the beginning of the file.');
    offset_idx = offset_idx(2:end);
end

if onset_idx(end) > offset_idx(end),
    disp('Removing orphaned offset at the beginning of the file.');
    onset_idx = onset_idx(1:end-1);
end

if strcmpi(output, 'onset'),
    eventlist = eventlist(onset_idx, :);
elseif strcmpi(output, 'step'),
    idx = union(onset_idx, offset_idx);
    eventlist = eventlist(idx);
elseif strcmpi(output, 'offset'),
    event_id = eventlist(onset_idx, 3);
    eventlist = eventlist(offset_idx, :);
    eventlist(:, 2) = eventlist(:, 3);
    eventlist(:, 3) = event_id;
    eventlist(:, 1) = eventlist(:, 1) - 1;
else
    error('Invalid output parameter: %s.', output);
end
