function tracks = read_mrtrix_tracks (filename)

% function: tracks = read_mrtrix_tracks (filename)
%
% returns a structure containing the header information and data for the MRtrix 
% format track file 'filename' (i.e. files with the extension '.tck'). 
% The track data will be stored as a cell array in the 'data' field of the
% return variable.

f = fopen (filename, 'r');
assert(f ~= -1, 'error opening %s', filename);
L = fgetl(f);
if ~strncmp(L, 'mrtrix tracks', 13)
  fclose(f);
  error('%s is not in MRtrix format', filename);
end

tracks = struct();

while 1
  L = fgetl(f);
  if ~ischar(L), break, end;
  L = strtrim(L);
  if strcmp(L, 'END'), break, end;
  d = strfind (L,':');
  if isempty(d)
    disp (['invalid line in header: ''' L ''' - ignored']);
  else
    key = lower(strtrim(L(1:d(1)-1)));
    value = strtrim(L(d(1)+1:end));
    if strcmp(key, 'file')
      file = value;
    elseif strcmp(key, 'datatype')
      tracks.datatype = value;
    else 
      tracks = add_field (tracks, key, value);
    end
  end
end
fclose(f);

assert(exist('file') && isfield(tracks, 'datatype'), ...
  'critical entries missing in header - aborting');

[ file, offset ] = strtok(file);
assert(strcmp(file, '.'), ...
  'unexpected file entry (should be set to current ''.'') - aborting');

assert(~isempty(offset), 'no offset specified - aborting');
offset = str2num(char(offset));

datatype = lower(tracks.datatype);
byteorder = datatype(end-1:end);

if strcmp(byteorder, 'le')
  f = fopen (filename, 'r', 'l');
  datatype = datatype(1:end-2);
elseif strcmp(byteorder, 'be')
  f = fopen (filename, 'r', 'b');
  datatype = datatype(1:end-2);
else
  error('unexpected data type - aborting');
end

assert(f ~= -1, 'error opening %s', filename);

fseek (f, offset, -1);
data = fread(f, inf, datatype);
fclose (f);

N = floor(prod(size(data))/3);
data = reshape (data, 3, N)';
k = find (~isfinite(data(:,1)));

tracks.data = {};
pk = 1;
for n = 1:(prod(size(k))-1)
  tracks.data{end+1} = data(pk:(k(n)-1),:);
  pk = k(n)+1;
end




function image = add_field (image, key, value)
  if isfield (image, key)
    previous = getfield (image, key);
    if iscell (previous)
      image = setfield (image, key, [ previous {value} ]);
    else
      image = setfield (image, key, { previous, value });
    end
  else
    image = setfield (image, key, value);
  end
