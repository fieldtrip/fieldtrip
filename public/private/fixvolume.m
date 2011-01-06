function data = fixvolume(data)

% FIXVOLUME cleans up the volume data representation, removes old and obsoleted
% fields and ensures that it is consistent with the most recent code.
%
% Use as
%   output = fixvolume(input)
% where input is a structure representing volume data
%
% See also FT_CHECKDATA, FIXSOURCE

% remove the unwanted fields
if isfield(data, 'pos'),     data = rmfield(data, 'pos');     end
if isfield(data, 'xgrid'),   data = rmfield(data, 'xgrid');   end
if isfield(data, 'ygrid'),   data = rmfield(data, 'ygrid');   end
if isfield(data, 'zgrid'),   data = rmfield(data, 'zgrid');   end

% volume data is always 3-D and therefore does not need a dimord
if isfield(data, 'dimord'),  data = rmfield(data, 'dimord');  end

% although volume data should be 3D, better do this conversion anyway
% old data structures may use latency/frequency instead of time/freq
if isfield(data, 'frequency'),
  data.freq = data.frequency;
  data      = rmfield(data, 'frequency');
end

if isfield(data, 'latency'),
  data.time = data.latency;
  data      = rmfield(data, 'latency');
end

hastime = isfield(data, 'time');
hasfreq = isfield(data, 'freq');

if hastime
  error('the volume data representation should be 3-D and is not allowed to contain a time field');
end

if hasfreq
  error('the volume data representation should be 3-D and is not allowed to contain a freq field');
end
