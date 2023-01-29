function [vertex_coords, faces, magic] = read_surf(fname, flag)
%
% [vertex_coords, faces] = read_surf(fname)
%
% or 
%
% [vertex_coords, faces] = read_surf(fname, flag)
%
% reads a the vertex coordinates and face lists from a surface file
% note that reading the faces from a quad file can take a very long
% time due to the goofy format that they are stored in. If the faces
% output variable is not specified, they will not be read so it
% should execute pretty quickly.
% 
% if a second optional boolean flag is specified to be true, this will determine whether
% the ASCII formatted coregistration information is read from the bottom of the file.
% This coregistration information is applied to the vertex_coords, to represent the 
% surface in the coordinate system of the corresponding volume. The function's default
% is not to use this coregistration information.

% changes:
%  - 20220518: also optionally read the transformation information that is present as
%    ASCII text at the bottom of the binary file. I could not find documentation
%    online about this, but based on some limited trial and error it seems that
%    the translation information (i.e. a shift of the origin), and possibly the 
%    voxel scaling can be used to represent the vertex coordinates in the coordinate
%    system of the volume from which the surface was derived. This option makes sense
%    if the surface is to be directly related to a volume, without the resulting 
%    headaches of mismatching coordinate systems.
%
% read_surf.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: fischl $
%    $Date: 2014/04/30 12:59:03 $
%    $Revision: 1.7 $
%
% Copyright (C) 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

if nargin<2
  flag = false;
end

%fid = fopen(fname, 'r') ;
%nvertices = fscanf(fid, '%d', 1);
%all = fscanf(fid, '%d %f %f %f %f\n', [5, nvertices]) ;
%curv = all(5, :)' ;

% open it as a big-endian file


%QUAD_FILE_MAGIC_NUMBER =  (-1 & 0x00ffffff) ;
%NEW_QUAD_FILE_MAGIC_NUMBER =  (-3 & 0x00ffffff) ;

TRIANGLE_FILE_MAGIC_NUMBER =  16777214 ;
QUAD_FILE_MAGIC_NUMBER =  16777215 ;
NEW_QUAD_FILE_MAGIC_NUMBER =  16777213 ;

fid = fopen(fname, 'rb', 'b') ;
if (fid < 0)
  str = sprintf('could not open surface file %s.', fname) ;
  error(str) ;
end
magic = fread3(fid) ;

if((magic == QUAD_FILE_MAGIC_NUMBER) || (magic == NEW_QUAD_FILE_MAGIC_NUMBER))
  vnum = fread3(fid) ;
  fnum = fread3(fid) ;
  vertex_coords = fread(fid, vnum*3, 'int16') ./ 100 ;
  if (nargout > 1)
    for i=1:fnum
      for n=1:4
        faces(i,n) = fread3(fid) ;
      end
    end
  end
elseif (magic == TRIANGLE_FILE_MAGIC_NUMBER)
  fgets(fid) ;
  fgets(fid) ;
  vnum = fread(fid, 1, 'int32') ;
  fnum = fread(fid, 1, 'int32') ;
  % possibly the comment line was not followed by two \n’s, this patch is suggested by J.M.Schoffelen, 20190502
  if round(abs(vnum))~=vnum || round(abs(fnum))~=fnum
    frewind(fid);
    fread3(fid);
    fgets(fid);
    vnum = fread(fid, 1, 'int32') ;
    fnum = fread(fid, 1, 'int32') ;
  end
  vertex_coords = fread(fid, vnum*3, 'float32') ;
  if (nargout > 1)
    faces = fread(fid, fnum*3, 'int32') ;
    faces = reshape(faces, 3, fnum)' ;
  end
else
  fprintf('ERROR: magic number %d unknown\n',magic);
  vertex_coords = [];
  faces = [];
  return;
end

vertex_coords = reshape(vertex_coords, 3, vnum)' ;

if flag
  % JM the below is based a bit on trial and error, but tries to extract
  % coregistration information, that may be present as plain text at the
  % bottom of the file. I could not find any documentation online about this,
  % but anecdotally this info can be used to align the mesh with the
  % specified volume. Note that only a translation seems to be required, I am
  % not sure about the voxelsize scaling
  tline = fgetl(fid);
  if ~isempty(strfind(tline, 'volume info valid'))
    while tline~=-1
      tline = fgetl(fid);
      
      % the following assumes that the text lines are of a structure: 
      % something = something else
      if ischar(tline)
        [v, rest] = strtok(tline);
        [equalsign, rest] = strtok(rest);
        if isequal(v, 'filename')
          [r, rest] = strtok(rest);
          eval(sprintf('%s = ''%s'';', v, rest));
        else
          eval(sprintf('%s = [%s];', v, rest));
        end
      end
    end
    R = diag(voxelsize); % not sure about this
    M = [R cras';0 0 0 1];
    tmp = [vertex_coords ones(size(vertex_coords,1),1)] * M';
    vertex_coords = tmp(:,1:3);
  end
end

fclose(fid);

