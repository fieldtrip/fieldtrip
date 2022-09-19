function image = read_mrtrix (filename)

% function: image = read_mrtrix (filename)
%
% returns a structure containing the header information and data for the MRtrix 
% format image 'filename' (i.e. files with the extension '.mif' or '.mih').

f = fopen (filename, 'r');
assert(f ~= -1, 'error opening %s', filename);
L = fgetl(f);
if ~strncmp(L, 'mrtrix image', 12)
  fclose(f);
  error('%s is not in MRtrix format', filename);
end

transform = [];
dw_scheme = [];

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
    key = strrep(key, '-', '_'); % replace unpermitted '-' character in matlab fieldnames with '_'
    key = strrep(key, ' ', '_'); % replace unpermitted ' ' 
    value = strtrim(L(d(1)+1:end));
    if strcmp(key, 'dim')
      image.dim = str2num(char(split_strings (value, ',')))';
    elseif strcmp(key, 'vox')
      image.vox = str2num(char(split_strings (value, ',')))';
    elseif strcmp(key, 'layout')
      image.layout = value;
    elseif strcmp(key, 'datatype')
      image.datatype = value;
    elseif strcmp(key, 'transform')
      transform(end+1,:) = str2num(char(split_strings (value, ',')))';
    elseif strcmp(key, 'file')
      file = value;
    elseif strcmp(key, 'dw_scheme')
      dw_scheme(end+1,:) = str2num(char(split_strings (value, ',')))';
    else 
      image = add_field (image, key, value);
    end
  end
end
fclose(f);


if ~isempty(transform)
  image.transform = transform;
  image.transform(4,:) = [ 0 0 0 1 ];
end

if ~isempty(dw_scheme)
  image.dw_scheme = dw_scheme;
end

if ~isfield (image, 'dim') || ~exist ('file') || ...
  ~isfield (image, 'layout') || ~isfield (image, 'datatype')
  error('critical entries missing in header - not reading data');
end

layout = split_strings(image.layout, ',');
order = (abs(str2num (char(layout)))+1)';

[ file, offset ] = strtok(file);
if isempty(offset), offset = 0; else offset = str2num(char(offset)); end
[f,g] = fileparts(filename);
if strcmp(file,'.'), file = filename; else file = fullfile (f, file); end

datatype = lower(image.datatype);
byteorder = datatype(end-1:end);

if strcmp(byteorder, 'le')
  f = fopen (file, 'r', 'l');
  datatype = datatype(1:end-2);
elseif strcmp(byteorder, 'be')
  f = fopen (file, 'r', 'b');
  datatype = datatype(1:end-2);
else 
  if strcmp(datatype, 'bit')
    datatype = 'bit1';
    f = fopen(file, 'r', 'b'); 
  else
    f = fopen(file, 'r'); 
  end
end

assert(f ~= -1, 'error opening %s', filename);

fseek (f, offset, -1);
image.data = fread (f, inf, datatype);
fclose (f);

order(order)= 1:size(order,2);
image.data = reshape (image.data, image.dim(order));
image.data = ipermute (image.data, order);
for i=1:size(order,2)
  if layout{i}(1) == '-'
    image.data = flipdim(image.data, i);
  end
end




function S = split_strings (V, delim)
  S = {};
  while size(V,2) > 0
    [R, V] = strtok(V,delim);
    S{end+1} = R;
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
