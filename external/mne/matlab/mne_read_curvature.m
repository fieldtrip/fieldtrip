function [curv] = mne_read_curvature(fname)
%
% [curf] = mne_read_surface(fname)
%
% Reads a FreeSurfer curvature file
%
% fname       - The file to read
% curv        - The curvature values
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.4  2006/07/31 05:12:39  msh
%   fread3 was still called instead of mne_fread3
%
%   Revision 1.3  2006/05/22 11:01:47  msh
%   Deleted superfluous text from the comments.
%
%   Revision 1.2  2006/05/22 10:55:02  msh
%   Fixed help text.
%
%   Revision 1.1  2006/05/22 10:44:44  msh
%   Added surface and curvature reading routines.
%   Fixed error in mne_read_source_spaces: triangle normals were not normalized
%

me='MNE:mne_read_curvature';

%
%   The input file is big endian
%
fid = fopen(fname,'rb','ieee-be');

if (fid < 0)
   error(me,'Cannot open file %s', fname);
end;

magic = mne_fread3(fid) ;
NEW_VERSION_MAGIC_NUMBER = 16777215;
if (magic == NEW_VERSION_MAGIC_NUMBER)
    nvert = fread(fid, 1, 'int32');
    nface = fread(fid, 1, 'int32');
    val_per_vertex = fread(fid, 1, 'int32');
    if val_per_vertex ~= 1
        fclose(fid);
        error(me,'Values per vertex not equal to one');
    end
    curv = fread(fid, nvert, 'float') ;
    fprintf(1,'\t%d values read from a new style curvature file.\n',nvert);
else
    nvert = magic;
    nface = mne_fread3(fid);
    curv = fread(fid, nvert, 'int16') ./ 100 ;
    fprintf(1,'\t%d values read from an old style curvature file.\n',nvert);
end

fclose(fid);

return;
