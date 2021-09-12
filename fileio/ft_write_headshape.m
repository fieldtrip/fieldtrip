function ft_write_headshape(filename, mesh, varargin)

% FT_WRITE_HEADSHAPE writes a head surface, cortical sheet or
% geometrical descrition of the volume conduction model or source
% model to a file for further processing in external software.
%
% Use as
%   ft_write_headshape(filename, mesh, ...)
% or
%   ft_write_headshape(filename, pos, ...)
% where the input mesh is a structure containing the vertices and triangles
% (mesh.pos and mesh.tri), or where the input pos is a Nx3 matrix that describes
% the surface vertices.
%
% Required input arguments should be specified as key-value pairs and
% should include
%   'format'		   = string, see below
%
% Optional input arguments should be specified as key-value pairs and
% can include
%   'data'         = data matrix, size(1) should be number of vertices
%   'unit'         = string, desired geometrical units for the data, for example 'mm'
%   'coordsys'     = string, desired coordinate system for the data
%
% Supported output formats are
%   'freesurfer'  Freesurfer surf-file format, using write_surf from FreeSurfer
%   'gifti'
%   'mne_pos'		  MNE source grid in ascii format, described as 3D points
%   'mne_tri'		  MNE surface desciption in ascii format
%   'off'
%   'ply'         Stanford Polygon file format, for use with Paraview or Meshlab
%   'stl'         STereoLithography file format, for use with CAD and generic 3D mesh editing programs
%   'tetgen'
%   'vista'
%   'vtk'         Visualization ToolKit file format, for use with Paraview
%
% See also FT_READ_HEADSHAPE

% Undocumented optional option:
%   'metadata'     = struct-array, containing the fields 'name' and
%                    'value', that will be used as metadata for the vertices.
%
% In the case of 'gifti' as required format you can optionally add a struct-array
% containing metadata. This is useful of the resulting files are to be used in
% combination with workbench. In this case it makes sense to specify the
% following names (with example values in brackets):
%
%   'AnatomicalStructurePrimary'   (e.g. 'CortexLeft'),
%   'AnatomicalStructureSecondary' (e.g. 'MidLayer')

% Copyright (C) 2011-2021, Lilla Magyari & Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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

fileformat    = ft_getopt(varargin, 'format');
data          = ft_getopt(varargin, 'data');         % can be stored in a gifti file
unit          = ft_getopt(varargin, 'unit');
coordsys      = ft_getopt(varargin, 'coordsys');
metadata      = ft_getopt(varargin, 'metadata');

% ensure that vertex positions are given in pos, not in pnt
mesh = fixpos(mesh);

if ~isempty(unit)
  % convert the input to the desired units prior to writing to disk
  mesh = ft_convert_units(mesh, unit);
end

if ~isempty(coordsys)
  % convert it to the specified coordinate system, this will also interactively determine the coordinate system if required
  mesh = ft_convert_coordsys(mesh, coordsys);
end

switch fileformat
  case 'mne_pos'
    fid = fopen_or_error(filename, 'wt');
    % convert to milimeter
    mesh = ft_convert_units(mesh, 'mm');
    n = size(mesh.pnt, 1);
    for line = 1:n
      num = mesh.pnt(line, 1);
      fprintf(fid, '%-1.0f ', num);
      num = mesh.pnt(line, 2);
      fprintf(fid, '%-1.0f ', num);
      num = mesh.pnt(line, 3);
      fprintf(fid, '%-1.0f\n', num);
    end
    fclose(fid);
    
  case 'mne_tri'
    fid = fopen_or_error(filename, 'wt');
    % convert to milimeter
    mesh = ft_convert_units(mesh, 'mm');
    n = size(mesh.pos, 1);
    fprintf(fid, '%-1.0f\n', n);
    for line = 1:n
      num=mesh.pos(line, 1);
      fprintf(fid, '%g ', num);
      num = mesh.pos(line, 2);
      fprintf(fid, '%g ', num);
      num = mesh.pos(line, 3);
      fprintf(fid, '%g\n', num);
    end
    n = size(mesh.tri, 1);
    fprintf(fid, '%-1.0f\n', n);
    for line = 1:n
      num=mesh.tri(line, 1);
      fprintf(fid, '%-1.0f ', num);
      num = mesh.tri(line, 2);
      fprintf(fid, '%-1.0f ', num);
      num = mesh.tri(line, 3);
      fprintf(fid, '%-1.0f\n', num);
    end
    fclose(fid);
    
  case 'off'
    write_off(filename, mesh.pos, mesh.tri);
    
  case 'vista'
    ft_hastoolbox('simbio', 1)
    % no conversion needed (works in voxel coordinates)
    if isfield(mesh, 'hex')
      write_vista_mesh(filename, mesh.pos, mesh.hex, mesh.index); % mesh.tensor
    elseif isfield(mesh, 'tet')
      write_vista_mesh(filename, mesh.pos, mesh.tet, mesh.index);
    else
      ft_error('unknown mesh representation')
    end
    
  case 'tetgen'
    % the third argument is the element type. At the moment only type 302=triangle is supported
    surf_to_tetgen(filename, mesh.pos, mesh.tri, 302*ones(size(mesh.tri, 1), 1), [], []);
    
  case 'vtk'
    [p, f, x] = fileparts(filename);
    filename = fullfile(p, [f, '.vtk']); % ensure it has the right extension
    if isfield(mesh, 'tri')
      write_vtk(filename, mesh.pos, mesh.tri);
    elseif isfield(mesh, 'tet')
      write_vtk(filename, mesh.pos, mesh.tet);
    elseif isfield(mesh, 'hex')
      write_vtk(filename, mesh.pos, mesh.hex);
    end
    
  case {'ply', 'ply_ascii', 'ply_binary'}
    [p, f, x] = fileparts(filename);
    filename = fullfile(p, [f, '.ply']); % ensure it has the right extension
    
    if isfield(mesh, 'pos')
      vertices = mesh.pos;
    elseif isfield(mesh, 'pos')
      vertices = mesh.pos;
    end
    
    if isfield(mesh, 'tri')
      elements = mesh.tri;
    elseif isfield(mesh, 'tet')
      elements = mesh.tet;
    elseif isfield(mesh, 'hex')
      elements = mesh.hex;
    end
    
    if length(fileformat)>4
      write_ply(filename, vertices, elements, fileformat(5:end));
    else
      write_ply(filename, vertices, elements, 'ascii');
    end
    
  case 'stl'
    %nrm = normals(mesh.pos, mesh.tri, 'triangle');
    %write_stl(filename, mesh.pos, mesh.tri, nrm);
    stlwrite(filename, mesh.tri, mesh.pos);
    
  case 'gifti'
    ft_hastoolbox('gifti', 1);
    
    mesh = ft_convert_units(mesh, 'mm');  % defined in the GIFTI standard to be milimeter
    
    % start with an empty structure
    tmp          = [];
    tmp.vertices = mesh.pos;
    tmp.faces    = mesh.tri;
    if ~isempty(data)
      tmp.cdata = data;
    end
    tmp = gifti(tmp);     % construct a gifti object
    
    % check the presence of metadata
    if ~isempty(metadata)
      if isstruct(metadata)
        fnames = fieldnames(metadata);
        if any(strcmp(fnames, 'name')) && any(strcmp(fnames, 'value'))
          % this is OK, and now assume the metadata to be written at the
          % level of the vertex info
          for k = 1:numel(tmp.private.data)
            n(k, 1) = size(tmp.private.data{k}.data, 1);
          end
          ix = find(n==size(tmp.vertices, 1));
          tmp.private.data{ix}.metadata = metadata;
        else
          ft_error('the metadata structure should contain the fields ''name'' and ''value''');
        end
      else
        ft_error('metadata should be provided as a struct-array');
      end
    end
        
    save(tmp, filename);  % write the object to file

  case 'freesurfer'
    ft_hastoolbox('freesurfer', 1);
    write_surf(filename, mesh.pos, mesh.tri);
    
  case []
    ft_error('no output format specified');
    
  otherwise
    ft_error('unsupported output format "%s"', fileformat);
end
