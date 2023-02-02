function [nodes, elements] = read_gmsh_binary(filename)

% READ_GMSH_BINARY reads a gmsh .msh binary file. Current support is only
% for version 2. There are some ASCII-readers floating around on the net,
% but they do not seem to work with the primary use case in FieldTrip (and
% the test data that I have available), which is SimNibs generated data.

fid = fopen(filename, 'rb');
txt = fgetl(fid);

if startsWith(txt, '$MeshFormat')
  % the file starts, as expected with the file header
  txt    = fgetl(fid);
  format = sscanf(txt, '%f');
  version  = format(1);
  isbinary = format(2)==1;
  dtype    = format(3);

  % read until the end of this chunk
  while 1
    txt = fgetl(fid);
    if startsWith(txt, '$EndMeshFormat')
      break;
    end
  end
else
  ft_error('File %s cannot be read, header information is missing', filename);
end

if floor(version)~=2
  ft_error('Currently only version 2 files can be read, the version appears to be: %f', version);
end
if ~isbinary
  ft_error('Currently only binary files can be read');
end

while 1
  txt = fgetl(fid);
  if txt==-1
    % end of file
    fclose(fid);
    break;
  end
  if startsWith(txt, '$Nodes')
    [indx, nodes] = getnodes(fid);
    nodes = nodes(indx,:);
  elseif startsWith(txt, '$Elements')
    [elements] = getelements(fid);
  else
    keyboard
  end
end

function [indx, nodes] = getnodes(fid)

% the chunk starts with an ascii-line with the number of elements
%
% the format of the nodes is 4-byte indx followed by 3x a float for the
% coordinates

txt = fgetl(fid);
N   = sscanf(txt, '%f');
ptr = ftell(fid);

fprintf('Reading position information for %d nodes\n', N);
indx = fread(fid, N, 'uint32', 24);
fseek(fid, ptr+4, 'bof');
nodes(:,1) = fread(fid, N, 'double', 20);
fseek(fid, ptr+12, 'bof');
nodes(:,2) = fread(fid, N, 'double', 20);
fseek(fid, ptr+20, 'bof');
nodes(:,3) = fread(fid, N, 'double', 20);
fseek(fid, ptr+28*N, 'bof');
txt = fgetl(fid);
assert(isequal(txt, '$EndNodes'), 'Reading of the nodes unexpectedly failed');

function [elements] = getelements(fid)

% the chunk starts with an ascii-line mentioning the overall number of elements, then 3 integers
% per sub-chunk of elements, that indicate the type, the number of elements
% in the sub-chunk, and the number of tags

txt = fgetl(fid);
N   = sscanf(txt, '%f');
ptr = ftell(fid);

fprintf('Reading element information for %d elements\n', N);
Ntotal = 0;
while N>Ntotal

  hdr  = fread(fid, 3, 'uint32');
  type = hdr(1);
  Nsub = hdr(2);
  ntag = hdr(3);

  switch type
    case 2
      % triangle
      nnode = 3;
      type_str = 'triangles';
    case 4
      % tetrahedron
      nnode = 4;
      type_str = 'tetrahedra';
    otherwise
      keyboard
  end
  
  fprintf('Reading %d %s with %d tags\n', Nsub, type_str, ntag);
  tmp = fread(fid, Nsub.*(nnode+ntag+1), 'uint32');
  tmp = reshape(tmp, nnode+ntag+1, [])';

  elements.(type_str) = tmp;
  Ntotal = Ntotal + Nsub;
end
txt = fgetl(fid);
assert(isequal(txt, '$EndElements'), 'Reading of the elements unexpectedly failed');
