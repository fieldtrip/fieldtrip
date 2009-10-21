function [shape] = read_ctf_shape(filename);

% READ_CTF_SHAPE reads headshape points and header information
% from a CTF *.shape teh accompanying *.shape_info file.
%
% Use as
%   [shape] = read_ctf_shape(filename)
% where filename should have the .shape extension 

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: read_ctf_shape.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.2  2008/11/12 16:59:51  roboos
% open explicitely as text using fopen 'rt'
%
% Revision 1.1  2003/10/08 15:28:03  roberto
% *** empty log message ***
%

shape = read_ctf_ascii([filename '_info']);

if ~strcmp(shape.MRI_Info.COORDINATES, 'HEAD')
  warning('points on head shape are NOT in headcoordinates')
end

fid = fopen(filename, 'rt');
num = fscanf(fid, '%d', 1);
shape.pnt = fscanf(fid, '%f', inf);
shape.pnt = reshape(shape.pnt, [3 num])';
fclose(fid);

