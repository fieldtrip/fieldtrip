function data = read_vtk_xml(filename)

% READ_VTK_XML reads a XML-formatted vtk file containing points in 3D and
% connecting elements.
%
% this function is a trial-and-error based implementation to read xml-style
% vtk files. There is some documentation online, which seems somewhat
% incomplete, or at least not fully understood by me.

tree   = xmltree(filename);
nnodes = length(tree);
for k = 1:nnodes
  this_node = get(tree, k);
  if isfield(this_node, 'attributes')
    node_info(k) = this_node;
  end
  if isfield(this_node, 'value')
    % this is a node that contains data, probably compressed
    node_info(k).name = 'n/a';
  end
end
nodenames = {node_info.name}';

info = attributes2struct(node_info(strcmp(nodenames,'VTKFile')));
iscompressed = isfield(info, 'compressor') && strcmp(info.compressor, 'vtkZLibDataCompressor');
if ~isequal(info.type, 'PolyData')
  error('at present this function only reads PolyData files');
end
% Note: the below code is a bit trial-and-error, and not fully generic,
% e.g. it does not take into account files with multiple Cells etc. 
% Once I encounter examples that have this, I will look into this.

% now extract all DataArrays, and put them into a struct
% this part is rather generic, but does not deal with the connecting
% elements (Lines etc.) well.
sel = find(strcmp(nodenames, 'DataArray'));
for k = 1:numel(sel)
  this_node = node_info(sel(k));
  this_info = attributes2struct(this_node);
  if isequal(this_info.format, 'appended')
    data_uid  = node_info(strcmp(nodenames,'AppendedData')).contents;
  else
    data_uid  = this_node.contents;
    
    %data_uid  = node_info(this_node.contents).contents;
  end
  data_node = get(tree, data_uid);
  this_data = value2data(data_node, this_info, iscompressed);
  data.(this_info.Name) = this_data;
end

% the struct N contains the basic information about the core stuff in the
% file, there can be additional data in other DataArrays
N    = attributes2struct(node_info(strcmp(nodenames,'Piece')),true);
fn   = fieldnames(N);
for k = 1:numel(fn)
  if N.(fn{k})~=0
    % there is data of a particular type: can be Points Verts Lines Strips
    % Polys, Points are X/Y/Z locations. Other elements seem to index into 
    % the points, and describe edges
    type = fn{k}(9:end);
    this_node = node_info(strcmp(nodenames, type));
    switch type 
      case 'Points'
        this_info = attributes2struct(node_info(this_node.contents));
        if isequal(this_info.format, 'appended')
          data_uid  = node_info(strcmp(nodenames,'AppendedData')).contents;
        else
          data_uid  = node_info(this_node.contents).contents;
        end
        data_node = get(tree, data_uid);
        this_data = value2data(data_node, this_info, iscompressed);
        data.Points  = this_data;
      case 'Lines'
        % these are defined by 2 data arrays, one representing the
        % connectivity and the other the offsets
        uids = this_node.contents;
        for m = 1:numel(uids)
          this_info = attributes2struct(node_info(uids(m)));
          if isequal(this_info.format, 'appended')
            data_uid = node_info(strcmp(nodenames,'AppendedData')).contents;
          else
            data_uid = node_info(uids(m)).contents;
          end
          data_node = get(tree, data_uid);
          if isequal(this_info.Name, 'connectivity')
            connectivity = value2data(data_node, this_info, iscompressed);
          elseif isequal(this_info.Name, 'offsets')
            offsets = value2data(data_node, this_info, iscompressed);
          end
        end
        % this is a shortcut creation of the lines, that assumes that the
        % ordering of the points (connectivity) is 0:(N-1)
        offsets = double([0;offsets]);
        npoints = diff(offsets);
        
        back_and_forth = true;
        lines   = zeros(max(npoints),numel(offsets)-1)+nan;
        if back_and_forth
          % 'draw' the line back and forth, this doubles the memory
          % requirement, but will allow the patch command (used in
          % ft_plot_mesh) to draw each connected line as a proper line,
          % rather than as a polygon
          lines = [lines;lines];
        end
        
        offset  = 0;
        for m = 1:size(lines,2)
          lines(1:npoints(m),m) = offset+(1:npoints(m))';
          if back_and_forth
            lines(npoints(m)+(1:npoints(m)),m) = flip(lines(1:npoints(m),m),1);
          end
          offset = nanmax(lines(:));
        end
        data.Lines = lines';
      otherwise
        error('not yet implemented for this type')
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function info = attributes2struct(node, castflag)

% subfunction that converts xmltree style cell-array of key-value structs
% to a matlab struct

if nargin<2
  castflag = false;
end

for k = 1:numel(node.attributes)
  if ~castflag
    info.(node.attributes{k}.key) = node.attributes{k}.val;
  else
    info.(node.attributes{k}.key) = str2double(node.attributes{k}.val);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = value2data(node, info, iscompressed)

if ~iscompressed
  error('reading of non-compressed data is not yet supported');
end

if isequal(info.format, 'binary')
  % this is OK
elseif isequal(info.format, 'appended')
  offset = str2double(info.offset)+2;
  node.value = node.value(offset:end); 
else
  error('non-binary data representation is not yet supported');
end

% this is essentially a bit trial-and-error, I don't know why the
% base64decoding is needed, but it seems to work in order to make sense of
% the data.
val = base64decode(node.value);
val = val(:);

% get the first 3 numbers for the header info (assuming 4 bytes per number)
% -> it could also be occasionally 8 bytes, but this would be in the xml
% header

nbytes_per_number = 4;
nblocks    = typecast(val(1:nbytes_per_number),'uint32');
headerinfo = typecast(val(1:(3+nblocks)*nbytes_per_number),'uint32');

offset = (3+nblocks)*nbytes_per_number+1;
tmpval = val(offset:end);
out = zeros(0,1);
for k = 1:nblocks
  % the 'relocation' of the start_idx is needed to get to the right start
  % of the block (otherwise java throuws an error)
  start_idx = find(tmpval,1,'first');
  end_idx   = start_idx+headerinfo(3+k)-1;
  this_block = zlibdecode(tmpval(start_idx:end_idx));
  switch info.type
    case 'Float32'
      this_block = typecast(this_block, 'single');
    case 'Int64'
      this_block = typecast(this_block, 'int64');
    case 'Int32'
      this_block = typecast(this_block, 'int32');
  end
  out = cat(1,out, this_block(:));
  
  offset = end_idx+1;
  tmpval = tmpval(offset:end);
end

if isfield(info, 'NumberOfComponents')
  nrow = str2double(info.NumberOfComponents);
else
  nrow = 1;
end
data = reshape(out,nrow,[]).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = base64decode(input)
%BASE64DECODE Decode Base64 string to a byte array.
%
%    output = base64decode(input)
%
% The function takes a Base64 string INPUT and returns a uint8 array
% OUTPUT. JAVA must be running to use this function. The result is always
% given as a 1-by-N array, and doesn't retrieve the original dimensions.
%
% See also base64encode
error(javachk('jvm'));
if ischar(input), input = uint8(input); end
output = typecast(org.apache.commons.codec.binary.Base64.decodeBase64(input), 'uint8')';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = zlibdecode(input)
%ZLIBDECODE Decompress input bytes using ZLIB.
%
%    output = zlibdecode(input)
%
% The function takes a compressed byte array INPUT and returns inflated
% bytes OUTPUT. The INPUT is a result of GZIPENCODE function. The OUTPUT
% is always an 1-by-N uint8 array. JAVA must be enabled to use the function.
%
% See also zlibencode typecast

error(javachk('jvm'));
if ischar(input)
  warning('zlibdecode:inputTypeMismatch', ...
          'Input is char, but treated as uint8.');
  input = uint8(input);
end
if ~isa(input, 'int8') && ~isa(input, 'uint8')
    error('Input must be either int8 or uint8.');
end

buffer = java.io.ByteArrayOutputStream();
zlib = java.util.zip.InflaterOutputStream(buffer);
zlib.write(input, 0, numel(input));
zlib.close();
output = typecast(buffer.toByteArray(), 'uint8')';

