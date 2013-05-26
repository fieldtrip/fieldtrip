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
% Required input arguments should be specified as key-value pairs and
% should include
%   format		= string, see below
% Optional input arguments should be specified as key-value pairs and
% can include
%   data      = data matrix, size(1) should be number of vertices
%   dimord    = string describing the dimensions of the data, e.g. 'pos_time'
%
% Supported output formats are
%   'mne_tri'		MNE surface desciption in ascii format
%   'mne_pos'		MNE source grid in ascii format, described as 3D points
%   'off'
%   'vista'
%   'tetgen'
%   'gifti'
%   'stl'       STereoLithography file format, for use with CAD and generic 3D mesh editing programs
%   'vtk'       Visualization ToolKit file format, for use with Paraview
%   'ply'       Stanford Polygon file format, for use with Paraview or Meshlab
%   'freesurfer' Freesurfer surf-file format, using write_surf from
%   Freesurfer
%
% See also FT_READ_HEADSHAPE

% Copyright (C) 2011, Lilla Magyari & Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

fileformat = ft_getopt(varargin,'format','unknown');
data       = ft_getopt(varargin,'data',  []); % can be stored in a gifti file
dimord     = ft_getopt(varargin,'dimord',[]); % optional for data

if ~isstruct(bnd)
  bnd.pnt = bnd;
end


switch fileformat
  case 'mne_pos'
    fid = fopen(filename, 'wt');
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
    fclose(fid);
    
  case 'mne_tri'
    fid = fopen(filename, 'wt');
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
    fclose(fid);
    
  case 'off'
    write_off(filename,bnd.pnt,bnd.tri);
    
  case 'vista'
    if ft_hastoolbox('simbio',1)
      % no conversion needed (works in voxel coordinates)
      if isfield(bnd,'hex')
        write_vista_mesh(filename,bnd.pnt,bnd.hex,bnd.index); % bnd.tensor
      elseif isfield(bnd,'tet')
        write_vista_mesh(filename,bnd.pnt,bnd.tet,bnd.index);
      else
        error('unknown format')
      end
    else
      error('You need Simbio/Vista toolbox to write a .v file')
    end
    
  case 'tetgen'
    % the third argument is the element type. At the moment only type 302
    % (triangle) is supported
    surf_to_tetgen(filename, bnd.pnt, bnd.tri, 302*ones(size(bnd.tri,1),1),[],[]);
  
  case 'vtk'
    [p, f, x] = fileparts(filename);
    filename = fullfile(p, [f, '.vtk']); % ensure it has the right extension
    if isfield(bnd, 'tri')
      write_vtk(filename, bnd.pnt, bnd.tri);
    elseif isfield(bnd, 'tet')
      write_vtk(filename, bnd.pnt, bnd.tet);
    elseif isfield(bnd, 'hex')
      write_vtk(filename, bnd.pnt, bnd.hex);
    end

  case 'ply'
    [p, f, x] = fileparts(filename);
    filename = fullfile(p, [f, '.ply']); % ensure it has the right extension
    
    if isfield(bnd, 'pnt')
      vertices = bnd.pnt;
    elseif isfield(bnd, 'pos')
      vertices = bnd.pos;
    end
    
    if isfield(bnd, 'tri')
      elements = bnd.tri;
    elseif isfield(bnd, 'tet')
      elements = bnd.tet;
    elseif isfield(bnd, 'hex')
      elements = bnd.hex;
    end
    
    write_ply(filename, vertices, elements);

  case 'stl'
    nrm = normals(bnd.pnt, bnd.tri, 'triangle');
    write_stl(filename, bnd.pnt, bnd.tri, nrm);

  case 'gifti'
    ft_hastoolbox('gifti', 1);
    tmp = [];
    tmp.vertices = bnd.pnt;
    tmp.faces    = bnd.tri;
    if ~isempty(data)
      tmp.cdata = data;
    end
    tmp = gifti(tmp);
    save(tmp, filename);
  
  case 'freesurfer'
    ft_hastoolbox('freesurfer', 1);
    write_surf(filename, bnd.pnt, bnd.tri);
 
  case []
    error('you must specify the output format');
    
  otherwise
    error('unsupported output format "%s"', fileformat);
end

