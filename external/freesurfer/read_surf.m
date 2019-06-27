function [vertex_coords, faces, magic] = read_surf(fname)
%
% [vertex_coords, faces] = read_surf(fname)
% reads a the vertex coordinates and face lists from a surface file
% note that reading the faces from a quad file can take a very long
% time due to the goofy format that they are stored in. If the faces
% output variable is not specified, they will not be read so it
% should execute pretty quickly.
%
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

if((magic == QUAD_FILE_MAGIC_NUMBER) | (magic == NEW_QUAD_FILE_MAGIC_NUMBER))
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
fclose(fid) ;
