function write_mrtrix_tracks (tracks, filename)

% function: write_mrtrix_tracks (tracks, filename)
%
% writes the track data stored as a cell array in the 'data' field of the
% tracks variable to the MRtrix format track file 'filename'. All other fields
% of the tracks variable will be written as text entries in the header, and are
% expected to supplied as character arrays.

assert(isfield(tracks, 'data'), ...
  'input tracks variable does not contain required ''data'' field');

assert(iscell(tracks.data), ...
  'input tracks.data variable should be a cell array');

f = fopen (filename, 'w', 'ieee-le');
assert(f ~= -1, 'error opening %s', filename);

fprintf (f, 'mrtrix tracks\ndatatype: Float32LE\ncount: %d\n', prod(size(tracks.data)));
names = fieldnames(tracks);
for i=1:size(names)
  if strcmpi (names{i}, 'data'), continue; end
  if strcmpi (names{i}, 'count'), continue; end
  if strcmpi (names{i}, 'datatype'), continue; end
  fprintf (f, '%s: %s\n', names{i}, getfield(tracks, names{i}));
end
data_offset = ftell (f) + 20;
fprintf (f, 'file: . %d\nEND\n', data_offset);

fwrite (f, zeros(data_offset-ftell(f),1), 'uint8');
for i = 1:prod(size(tracks.data))
  fwrite (f, tracks.data{i}', 'float32');
  fwrite (f, [ nan nan nan ], 'float32');
end

fwrite (f, [ inf inf inf ], 'float32');
fclose (f);
