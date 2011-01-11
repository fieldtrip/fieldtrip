function ft_write_headshape(filename, bnd, varargin)

% FT_WRITE_HEADSHAPE writes a head surface, cortical sheet or
% geometrical descrition of the volume conduction model or source
% model to a file for further processing in external software.
%
% Use as
%   ft_write_headshape(filename, bnd, ...)
% or
%   ft_write_headshape(filename, pos, ...)
% where the bnd is a structure containing the vertices and triangles
% (bnd.pnt and bnd.tri), or where pnt describes the surface or source
% points.
%
% Optional input arguments should be specified as key-value pairs and
% should include
%   format		= string, see below
%
% Supported output formats are
%   'mne_tri'		MNE surface desciption in ascii format
%   'mne_pos'		MNE source grid in ascii format, described as 3D points
%
% See also FT_READ_HEADSHAPE

% Copyright (C) 2011, Lilla Magyari & Robert Oostenveld
%
% $Rev$

fileformat = keyval('format', varargin);

if ~isstruct(bnd)
  bnd.pnt = bnd;
end

fid = fopen(filename, 'wt');

switch fileformat
  case 'mne_pos'
    % convert to milimeter
    bnd = ft_convert_units(bnd, 'mm');
    n=size(bnd.pnt,1);
    for line = 1:n
      num = bnd.pnt(line,1);
      fprintf(fid, '%-1.0f ',num);
      num = bnd.pnt(line,2);
      fprintf(fid, '%-1.0f ',num);
      num = bnd.pnt(line,3);
      fprintf(fid, '%-1.0f\n',num);
    end
    
  case 'mne_tri'
    % convert to milimeter
    bnd = ft_convert_units(bnd, 'mm');
    n=size(bnd.pnt,1);
    fprintf(fid, '%-1.0f\n',n);
    for line = 1:n
      num=bnd.pnt(line,1);
      fprintf(fid,'%g ', num);
      num = bnd.pnt(line,2);
      fprintf(fid, '%g ',num);
      num = bnd.pnt(line,3);
      fprintf(fid, '%g\n',num);
    end
    n=size(bnd.tri,1);
    fprintf(fid, '%-1.0f\n',n);
    for line = 1:n
      num=bnd.tri(line,1);
      fprintf(fid,'%-1.0f ', num);
      num = bnd.tri(line,2);
      fprintf(fid, '%-1.0f ',num);
      num = bnd.tri(line,3);
      fprintf(fid, '%-1.0f\n',num);
    end
    
  case []
    error('you must specify the output format');
    
  otherwise
    error('unsupported output format "%s"');
end

fclose(fid);
