function [vert, face] = read_ply(fn)

% READ_PLY reads triangles, tetraheders or hexaheders from a Stanford *.ply file
%
% Use as
%   [vert, face, prop, face_prop] = read_ply(filename)
%
% Documentation is provided on
%   http://paulbourke.net/dataformats/ply/
%   http://en.wikipedia.org/wiki/PLY_(file_format)
%
% See also WRITE_PLY, WRITE_VTK, READ_VTK

% Copyright (C) 2013, Robert Oostenveld
%
% $Id$

fid = fopen(fn, 'r');
if fid~=-1
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % the file starts with an ascii header
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  line = readline(fid);
  if ~strcmp(line, 'ply')
    fclose(fid);
    error('unexpected header line');
  end
  line = readline(fid);
  
  if ~strncmp(line, 'format', 6)
    fclose(fid);
    error('unexpected header line');
  else
    format = line(8:end);
  end
  line = readline(fid);
  
  if ~strncmp(line, 'element vertex', 14)
    fclose(fid);
    error('unexpected header line');
  else
    nvert = str2double(line(16:end));
  end
  line = readline(fid);
  
  prop = [];
  while strncmp(line, 'property', 8)
    tok = tokenize(line);
    prop(end+1).format = tok{2};
    prop(end  ).name   = tok{3};
    line = readline(fid);
  end
  
  if ~strncmp(line, 'element face', 12)
    fclose(fid);
    error('unexpected header line');
  else
    nface = str2double(line(14:end));
  end
  line = readline(fid);
  
  if ~strcmp(line, 'property list uchar int vertex_index') && ~strcmp(line, 'property list uchar int vertex_indices')
    % the wikipedia documentation specifies vertex_index, but the OPTOCAT files
    % have vertex_indices
    
    % it would not be very difficult to enhance the reader here with another
    % representation of the faces, i.e. something else than "uchar int"
    fclose(fid);
    error('unexpected header line');
  end
  line = readline(fid);
  
  while ~strcmp(line, 'end_header');
    line = readline(fid);
  end
  
  offset = ftell(fid);
  fclose(fid);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % the file continues with a section of data, which can be ascii or binary
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  switch format
    
    case 'ascii 1.0'
      
      fid = fopen(fn, 'rt');
      fseek(fid, offset, 'cof');
      dat = fscanf(fid,'%f',[numel(prop), nvert])';
      for j=1:length(prop)
        vert.(prop(j).name) = dat(:,j);
      end
      face = zeros(nface,0);
      num = zeros(nface,1);
      for i=1:nface
        % each polygon can have a different number of elements
        num(i) = fscanf(fid,'%f',1);
        face(i,1:num(i)) = fscanf(fid,'%f',num(i));
      end
      fclose(fid);
      
    case 'binary_little_endian 1.0'
      
      fid = fopen(fn, 'rb', 'l');
      fseek(fid, offset, 'cof');
      
      dat = zeros(nvert,length(prop));
      for i=1:nvert
        for j=1:length(prop)
          dat(i,j) = fread(fid, 1, prop(j).format);
        end % for each property
      end % for each vertex
      for j=1:length(prop)
        switch prop(j).format
          % the format can be one of the following: char uchar short ushort int uint float double int8 uint8 int16 uint16 int32 uint32 float32 float64
          case 'char'
            vert.(prop(j).name) = uint8(dat(:,j));
          case 'uchar'
            vert.(prop(j).name) = uint8(dat(:,j));
          case 'short'
            vert.(prop(j).name) = int16(dat(:,j));
          case 'ushort'
            vert.(prop(j).name) = uint16(dat(:,j));
          otherwise
            vert.(prop(j).name) = dat(:,j);
        end
      end
      clear dat;
      
      face = zeros(nface,0);
      num  = zeros(nface,1);
      for i=1:nface
        % each polygon can have a different number of elements
        num(i) = fread(fid, 1, 'uint8');
        face(i,1:num(i)) = fread(fid, num(i), 'int32');
      end
      fclose(fid);
      
    case 'binary_big_endian 1.0'
      % this is exactly the same as the code above, except that the file is opened as big endian
      
      fid = fopen(fn, 'rb', 'b');
      fseek(fid, offset, 'cof');
      
      dat = zeros(nvert,length(prop));
      for i=1:nvert
        for j=1:length(prop)
          dat(i,j) = fread(fid, 1, prop(j).format);
        end % for each property
      end % for each vertex
      for j=1:length(prop)
        switch prop(j).format
          % the format can be one of the following: char uchar short ushort int uint float double int8 uint8 int16 uint16 int32 uint32 float32 float64
          case 'char'
            vert.(prop(j).name) = uint8(dat(:,j));
          case 'uchar'
            vert.(prop(j).name) = uint8(dat(:,j));
          case 'short'
            vert.(prop(j).name) = int16(dat(:,j));
          case 'ushort'
            vert.(prop(j).name) = uint16(dat(:,j));
          otherwise
            vert.(prop(j).name) = dat(:,j);
        end
      end
      clear dat;
      
      face = zeros(nface,0);
      num  = zeros(nface,1);
      for i=1:nface
        % each polygon can have a different number of elements
        num(i) = fread(fid, 1, 'uint8');
        face(i,1:num(i)) = fread(fid, num(i), 'int32');
      end
      fclose(fid);
      
    otherwise
      error('unsupported format');
  end % switch
  
else
  error('unable to open file');
end

% each polygon can have a different number of elements
% mark all invalid entries with nan
if any(num<size(face,2))
  for i=1:nface
    face(i,(num(i)+1):end) = nan;
  end
end % if

if numel(face)>0
  % MATLAB indexes start at 1, inside the file they start at 0
  face = face+1;
end

end % function read_ply


function line = readline(fid)
% read the next line from the ascii header, skip all comment lines
line = fgetl(fid);
while strncmp(line, 'comment', 7);
  line = fgetl(fid);
end
end % function readline

