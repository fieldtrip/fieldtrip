function write_mrtrix (image, filename)

% function: write_mrtrix (image, filename)
%
% write the data contained in the structure 'image' in the MRtrix 
% format image 'filename' (i.e. files with the extension '.mif' or '.mih').
%
% 'image' is either a N-dimensional array (N <= 16), or a structure containing
% the following fields:
%    image.data:            a N-dimensional array (N <= 16)
%    image.vox:             N-vector of voxel sizes (in mm) (default: { 2 }) [optional]
%    image.comments:        a cell array of strings [optional]
%    image.datatype:        the datatype specifier (default: float32) [optional]
%    image.mrtrix_version:  a character array [optional]
%    image.transform:       a 4x4 matrix [optional]
%    image.dw_scheme:       a NDWx4 matrix of gradient directions [optional]

fid = fopen (filename, 'w');
assert(fid ~= -1, 'error opening %s', filename);
fprintf (fid, 'mrtrix image\ndim: ');

if isstruct(image)
  dim = size(image.data);
else
  dim = size(image);
end

while prod(size(dim)) < 3
  dim(end+1) = 1;
end 

fprintf (fid, '%d', dim(1));
fprintf (fid, ',%d', dim(2:end));

fprintf (fid, '\nvox: ');
if isstruct (image) && isfield (image, 'vox')
  fprintf (fid, '%.3f', image.vox(1));
  fprintf (fid, ',%.3f', image.vox(2:end)); 
else
  fprintf(fid, '2');
  fprintf(fid, ',%d', 2*ones(1,size(dim,2)-1));
end

fprintf (fid, '\nlayout: +0');
fprintf (fid, ',+%d', 1:(size(dim,2)-1));

[computerType, maxSize, endian] = computer;
if isstruct (image) && isfield (image, 'datatype')
  datatype = lower(image.datatype);
  byteorder = datatype(end-1:end);

  if strcmp (byteorder, 'le')
    precision = datatype(1:end-2);
    byteorder = 'l';
  elseif strcmp(byteorder, 'be')
    precision = datatype(1:end-2);
    byteorder = 'b';
  else 
    if strcmp(datatype, 'bit')
      precision = 'bit1';
      byteorder = 'n';
    elseif strcmp (datatype, 'int8') || strcmp (datatype, 'uint8')
      precision = datatype;
      byteorder = 'n';
      if endian == 'L'
        datatype(end+1:end+3) = 'le';
      else
        datatype(end+1:end+3) = 'be';
      end 
    end
  end
else
  if endian == 'L'
    datatype = 'float32le';
  else 
    datatype = 'float32be';
  end 
  precision = 'float32';
  byteorder = 'n';
end
fprintf (fid, [ '\ndatatype: ' datatype ]);

fprintf (fid, '\nmrtrix_version: %s', 'matlab');

if isstruct (image) && isfield (image, 'transform')
  fprintf (fid, '\ntransform: %.6f', image.transform(1,1));
  fprintf (fid, ',%.6f', image.transform(1,2:4));
  fprintf (fid, '\ntransform: %6f', image.transform(2,1));
  fprintf (fid, ',%.6f', image.transform(2,2:4));
  fprintf (fid, '\ntransform: %.6f', image.transform(3,1));
  fprintf (fid, ',%.6f', image.transform(3,2:4));
end

if isstruct (image) && isfield (image, 'dw_scheme')
  for i=1:size(image.dw_scheme,1)
    fprintf (fid, '\ndw_scheme: %.6f', image.dw_scheme(i,1));
    fprintf (fid, ',%.6f', image.dw_scheme(i,2:4));
   end
end

% write out any other fields:
if isstruct (image)
  f = fieldnames(image);
  % don't worry about fields that have already been dealt with:
  f(strcmp(f,'dim')) = [];
  f(strcmp(f,'vox')) = [];
  f(strcmp(f,'layout')) = [];
  f(strcmp(f,'datatype')) = [];
  f(strcmp(f,'transform')) = [];
  f(strcmp(f,'dw_scheme')) = [];
  f(strcmp(f,'data')) = [];
  
  % write out contents of the remainging fields:
  for n = 1:numel(f)
    val = getfield (image, f{n});
    if iscell (val)
      for i=1:numel(val)
        fprintf (fid, '\n%s: %s', f{n}, val{i});
      end
    else
      fprintf (fid, '\n%s: %s', f{n}, val);
    end
  end
end
 

if strcmp(filename(end-3:end), '.mif')
  datafile = filename;
  dataoffset = ftell (fid) + 24;
  fprintf (fid, '\nfile: . %d\nEND\n                         ', dataoffset);
elseif strcmp(filename(end-3:end), '.mih')
  datafile = [ filename(end-3:end) '.dat' ];
  dataoffset = 0;
  fprintf (fid, '\nfile: %s %d\nEND\n', datafile, dataoffset);
else
  fclose(fid);
  error('unknown file suffix - aborting');
end

fclose(fid);

fid = fopen (datafile, 'r+', byteorder);
fseek (fid, dataoffset, -1);

if isstruct(image)
  fwrite (fid, image.data, precision);
else
  fwrite (fid, image, precision);
end
fclose (fid);
