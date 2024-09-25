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
    mesh = ft_read_headshape(filename); % this will automatically detect ascii or binary
    S    = standard_cond;

    if isempty(meshtype)
      % set the default
      if isfield(mesh, 'tet') || isfield(mesh, 'hex')
        meshtype = 'volume';
      else
        meshtype = 'surface';
      end
    end

    switch meshtype
      case 'volume'
        if isfield(mesh, 'tet')
          % prune the mesh
          [headmodel.pos, headmodel.tet] = remove_unused_vertices(mesh.pos, mesh.tet);

          tag    = mesh.tetrahedron_regions(:,1);
          utag   = unique(tag);
          tissue = zeros(size(tag));
          tissuelabel = cell(numel(utag),1);
          cond   = zeros(1,numel(utag));
          for k = 1:numel(utag)
            tissue(tag==utag(k)) = k;
            tissuelabel{k} = lower(S(utag(k)).name);
            cond(k) = S(utag(k)).value;
          end

          headmodel.tissue = tissue;
          headmodel.tissuelabel = tissuelabel;
          headmodel.cond = cond;
          headmodel.type = 'simnibs';
          headmodel = ft_determine_units(headmodel);
        end

        if isfield(mesh, 'hex')
          % prune the mesh
          [headmodel.pos, headmodel.hex] = remove_unused_vertices(mesh.pos, mesh.hex);

          tag    = mesh.hexahedron_regions(:,1);
          utag   = unique(tag);
          tissue = zeros(size(tag));
          tissuelabel = cell(numel(utag),1);
          cond   = zeros(1,numel(utag));
          for k = 1:numel(utag)
            tissue(tag==utag(k)) = k;
            tissuelabel{k} = lower(S(utag(k)).name);
            cond(k) = S(utag(k)).value;
          end

          headmodel.tissue = tissue;
          headmodel.tissuelabel = tissuelabel;
          headmodel.cond = cond;
          headmodel.type = 'simnibs';
          headmodel = ft_determine_units(headmodel);
        end

      case 'surface'
        % prune the mesh
        [pos, tri] = remove_unused_vertices(mesh.pos, mesh.tri);

        tag  = mesh.triangle_regions(:,1);
        utag = unique(tag);

        cond = zeros(1,numel(utag));
        tissuelabel = cell(numel(utag),1);
        for k = 1:numel(utag)

          tissuelabel{k} = lower(S(utag(k)-1000).name);
          cond(k) = S(utag(k)-1000).value;
          ft_info('Creating boundary for tissue type %s', tissuelabel{k});

          seltri = tri(tag==utag(k), :);
          usepos = unique(reshape(seltri,[],1));
          removepos = setdiff((1:size(pos,1))', usepos);
          [bnd(k,1).pos, bnd(k,1).tri] = remove_vertices(pos, tri, removepos);
        end
        headmodel.bnd  = bnd;
        headmodel.tissuelabel = tissuelabel;
        headmodel.cond = cond;
        headmodel.type = 'simnibs';

      otherwise
        ft_error('unsupported meshtype');
    end

  otherwise
    ft_error('unknown fileformat for volume conductor model');
end

% this will ensure that the structure is up to date, e.g. that the type is correct and that it has units
headmodel = ft_datatype_headmodel(headmodel);

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
