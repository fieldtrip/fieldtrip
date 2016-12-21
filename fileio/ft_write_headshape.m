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
% (bnd.pnt and bnd.tri), or where pos is a Nx3 matrix that describes the 
% surface or source points.
%
% Required input arguments should be specified as key-value pairs and
% should include
%   'format'		  = string, see below
%
% Optional input arguments should be specified as key-value pairs and
% can include
%   'data'        = data matrix, size(1) should be number of vertices
%   'unit'        = string, e.g. 'mm'
%
% Supported output formats are
%   'mne_tri'		MNE surface desciption in ascii format
%   'mne_pos'		MNE source grid in ascii format, described as 3D points
%   'off'
%   'vista'
%   'tetgen'
%   'gifti'
%   'stl'           STereoLithography file format, for use with CAD and generic 3D mesh editing programs
%   'vtk'           Visualization ToolKit file format, for use with Paraview
%   'ply'           Stanford Polygon file format, for use with Paraview or Meshlab
%   'freesurfer'    Freesurfer surf-file format, using write_surf from FreeSurfer
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

% Copyright (C) 2011-2014, Lilla Magyari & Robert Oostenveld
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

fileformat    = ft_getopt(varargin, 'format', 'unknown');
data          = ft_getopt(varargin, 'data');         % can be stored in a gifti file
unit          = ft_getopt(varargin, 'unit');
metadata      = ft_getopt(varargin, 'metadata');

if ~isfield(bnd, 'pnt') && isfield(bnd, 'pos')
  bnd.pnt = bnd.pos;
end

if ~isstruct(bnd)
  bnd.pnt = bnd;
end

if ~isempty(unit)
  % convert to the desired units prior to writing to disk
  bnd = ft_convert_units(bnd, unit);
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
    
  case {'ply', 'ply_ascii', 'ply_binary'}
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
    
    if length(fileformat)>4
      write_ply(filename, vertices, elements, fileformat(5:end));
    else
      write_ply(filename, vertices, elements, 'ascii');
    end
    
  case 'stl'
    nrm = normals(bnd.pnt, bnd.tri, 'triangle');
    write_stl(filename, bnd.pnt, bnd.tri, nrm);
    
  case 'gifti'
    ft_hastoolbox('gifti', 1);
    
    bnd = ft_convert_units(bnd, 'mm');  % defined in the GIFTI standard to be milimeter
    
    % start with an empty structure
    tmp          = [];
    tmp.vertices = bnd.pnt;
    tmp.faces    = bnd.tri;
    if ~isempty(data)
      tmp.cdata = data;
    end
    tmp = gifti(tmp);     % construct a gifti object
    
    % check the presence of metadata
    if ~isempty(metadata)
      if isstruct(metadata)
        fnames = fieldnames(metadata);
        if any(strcmp(fnames,'name')) && any(strcmp(fnames,'value'))
          % this is OK, and now assume the metadata to be written at the
          % level of the vertex info
          for k = 1:numel(tmp.private.data)
            n(k,1) = size(tmp.private.data{k}.data,1);
          end
          ix = find(n==size(tmp.vertices,1));
          tmp.private.data{ix}.metadata = metadata;
        else
          error('the metadata structure should contain the fields ''name'' and ''value''');
        end
      else
        error('metadata should be provided as a struct-array');
      end
    end
        
    save(tmp, filename);  % write the object to file

  case 'freesurfer'
    ft_hastoolbox('freesurfer', 1);
    write_surf(filename, bnd.pnt, bnd.tri);
    
  case []
    error('you must specify the output format');
    
  otherwise
    error('unsupported output format "%s"', fileformat);
end

