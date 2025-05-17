function [headmodel] = ft_read_headmodel(filename, varargin)

% FT_READ_HEADMODEL reads a head model (or volume conduction model of the head) from
% various manufacturer specific files. Currently supported are ASA, CTF, Neuromag,
% MBFYS, MATLAB and SimNIBS. The volume conduction model is represented as a
% structure with fields that depend on the type of model.
%
% Use as
%   headmodel = ft_read_headmodel(filename, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'fileformat' = string
%
% If the fileformat is 'simnibs', an additional options can be used to specify the
% type of model that is to be returned.
%   'meshtype'   = string, 'volume' or 'surface' (default is automatic)
%
% See also FT_DATATYPE_HEADMODEL, FT_PREPARE_HEADMODEL, FT_READ_HEADMODEL,
% FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% Copyright (C) 2008-2024, Robert Oostenveld
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

% optionally get the data from the URL and make a temporary local copy
filename = fetch_url(filename);

% test whether the file exists
if ~exist(filename, 'file')
  ft_error('file ''%s'' does not exist', filename);
end

% get the options
fileformat = ft_getopt(varargin, 'fileformat', ft_filetype(filename));
meshtype   = ft_getopt(varargin, 'meshtype'); % the default is handled further down

switch fileformat
  case 'matlab'
    headmodel = loadvar(filename, 'headmodel');
    
  case 'ctf_hdm'
    headmodel = read_ctf_hdm(filename);
    headmodel.coordsys = 'ctf';
    
  case 'asa_vol'
    headmodel = read_asa_vol(filename);
    headmodel.type = 'asa';
    
  case 'mbfys_ama'
    ama = loadama(filename);
    headmodel = ama2headmodel(ama);
    
  case 'neuromag_fif'
    ft_hastoolbox('mne', 1);
    global FIFF
    bem = mne_read_bem_surfaces(filename);
    headmodel.bnd.pos = bem.rr;
    headmodel.bnd.tri = bem.tris;
    headmodel.coordsys = fif2coordsys(bem.coord_frame);
    
  case 'gmsh_binary'
    ft_error('this could be a simnibs headmodel, please specify ''fileformat'' to be ''simnibs'' if this is the case');
    
  case {'simnibs' 'simnibs_v3' 'simnibs_v4'}
    ft_hastoolbox('simnibs', 1);
    standard = standard_cond;
    mesh = ft_read_headshape(filename, 'meshtype', meshtype); % this will automatically detect ascii or binary
    
    hastri = isfield(mesh, 'tri');
    hastet = isfield(mesh, 'tet');
    hashex = isfield(mesh, 'hex');
    
    if (hastri+hastet+hashex)>1
      ft_error('multiple mesh types are not supported in a single head model');
    end
    
    if hastri
      % convert the long list of triangles into multiple surfaces
      % this also prunes unused vertices
      [bnd, bndtissue] = tri2bnd(mesh.pos, mesh.tri, mesh.tissue);
      headmodel.bnd = bnd;
      headmodel.tissuelabel = mesh.tissuelabel(bndtissue);
      
      headmodel.cond = nan(size(headmodel.tissuelabel));
      for i=1:numel(headmodel.tissuelabel)
        for j=1:numel(standard)
          if strcmpi(standard(j).type, 'COND') && ~isempty(standard(j).name) && strcmpi(standard(j).name, headmodel.tissuelabel{i})
            headmodel.cond(i) = standard(j).value;
            break
          end
        end
      end
      
    elseif hastet
      % prune unused vertices
      [headmodel.pos, headmodel.tet] = remove_unused_vertices(mesh.pos, mesh.tet);
      
      % renumber the tissues to remove unused tissue types
      tag    = mesh.tissue(:);
      utag   = unique(tag);
      tissue = zeros(size(tag));
      tissuelabel = cell(numel(utag),1);
      cond   = zeros(1,numel(utag));
      for k = 1:numel(utag)
        tissue(tag==utag(k)) = k;
        tissuelabel{k} = lower(standard(utag(k)).name);
        cond(k) = standard(utag(k)).value;
      end
      
      headmodel.tissue = tissue;
      headmodel.tissuelabel = tissuelabel;
      headmodel.cond = cond;
      
    elseif hashex
      % prune unused vertices
      [headmodel.pos, headmodel.hex] = remove_unused_vertices(mesh.pos, mesh.hex);
      
      % renumber the tissues to remove unused tissue types
      tag    = mesh.tissue(:);
      utag   = unique(tag);
      tissue = zeros(size(tag));
      tissuelabel = cell(numel(utag),1);
      cond   = zeros(1,numel(utag));
      for k = 1:numel(utag)
        tissue(tag==utag(k)) = k;
        tissuelabel{k} = lower(standard(utag(k)).name);
        cond(k) = standard(utag(k)).value;
      end
      
      headmodel.tissue = tissue;
      headmodel.tissuelabel = tissuelabel;
      headmodel.cond = cond;
      
    end % if tri, tet or hex
    
  otherwise
    ft_error('unknown fileformat for volume conductor model');
end % switch fileformat

% this will ensure that the structure is up to date, e.g. that the type is correct and that it has units
headmodel = ft_datatype_headmodel(headmodel);
headmodel = ft_determine_units(headmodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coordsys = fif2coordsys(coord_frame)
ft_hastoolbox('mne', 1);
global FIFF

switch coord_frame
  case FIFF.FIFFV_COORD_DEVICE
    coordsys = 'device';
  case FIFF.FIFFV_COORD_HPI
    coordsys = 'hpi';
  case FIFF.FIFFV_COORD_HEAD
    coordsys = 'head';
  case FIFF.FIFFV_COORD_MRI
    coordsys = 'mri';
  case FIFF.FIFFV_COORD_MRI_SLICE
    coordsys = 'mri_slice';
  case FIFF.FIFFV_COORD_MRI_DISPLAY
    coordsys = 'mri_display';
  case FIFF.FIFFV_COORD_DICOM_DEVICE
    coordsys = 'dicom_device';
  case FIFF.FIFFV_COORD_IMAGING_DEVICE
    coordsys = 'imaging_device';
  otherwise
    error('unrecognized coord_frame')
end
