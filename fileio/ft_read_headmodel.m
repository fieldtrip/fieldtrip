function [headmodel] = ft_read_headmodel(filename, varargin)

% FT_READ_HEADMODEL reads a volume conduction model from various manufacturer
% specific files. Currently supported are ASA, CTF, Neuromag, MBFYS, Matlab
% and simnibs.
%
% Use as
%   headmodel = ft_read_headmodel(filename, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'fileformat'   string
%
% The volume conduction model is represented as a structure with fields
% that depend on the type of model.
%
% If the fileformat is 'simnibs', an additional options should specify the
% type of model that is to be returned.
%   'meshtype'     string 'volume' (default) or 'surface'
%
% See also FT_DATATYPE_HEADMODEL, FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% Copyright (C) 2008-2018 Robert Oostenveld
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
meshtype   = ft_getopt(varargin, 'meshtype', 'volume');

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
    if ~startsWith(ft_filetype(filename), 'gmsh')
      ft_error('a simnibs headmodel should have the extension .gmsh');
    end
    mesh = ft_read_headshape(filename);
    S    = standard_cond;

    switch meshtype
      case 'volume'
        if isfield(mesh, 'tet')
          npos = numel(unique(mesh.tet(:)));
          if npos~=size(mesh.pos,1)
            ft_error('pruning of vertices is not yet supported, please file PR on github if needed');
          end
          headmodel.pos = mesh.pos;
          headmodel.tet = mesh.tet;
          
          tag    = mesh.tag_tet(:,1);
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
          headmodel = ft_determine_units(headmodel);
        end
        if isfield(mesh, 'hex')
          npos = numel(unique(mesh.hex(:)));
          if npos~=size(mesh.pos,1)
            ft_error('pruning of vertices is not yet supported, please file PR on github if needed');
          end
          headmodel.pos = mesh.pos;
          headmodel.hex = mesh.hex;
          
          tag    = mesh.tag_hex(:,1);
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
          headmodel = ft_determine_units(headmodel);
        end
      case 'surface'
        npos = numel(unique(mesh.tri(:)));
        if npos~=size(mesh.pos,1)
          ft_info('pruning the vertices, not all of them are used for the triangulations');
          usepos = unique(reshape(mesh.tri,[],1));
          removepos = setdiff((1:size(mesh.pos,1))', usepos);
          [pos, tri] = remove_vertices(mesh.pos, mesh.tri, removepos);
        else
          pos = mesh.pos;
          tri = mesh.tri;
        end
        tag  = mesh.tag_tri(:,1);
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
        headmodel.cond = cond;
        headmodel.tissuelabel = tissuelabel;
      otherwise
        ft_error('unsupported meshtype requested');
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

%%%%%%%
% the below is from the simnibs repo
function S=standard_cond
% USAGE: S=standard_cond
% 
%        returns a structure with standard conductivity values
%
% A. Thielscher, A. Antunes 2013
% updated 2018; A. Thielscher, G. Saturnino


%    This program is part of the SimNIBS package.
%    Please check on www.simnibs.org how to cite our work in publications.
%
%    Copyright (C) 2013-2018 Axel Thielscher, Andre Antunes, Guilherme Saturnino
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

% WM
S(1)=sim_struct('COND');
S(1).name = 'WM';
S(1).value = 0.126;
S(1).descrip = 'brain white matter (from Wagner 2004)';

% GM
S(2)=sim_struct('COND');
S(2).name = 'GM';
S(2).value = 0.275;
S(2).descrip = 'brain gray matter (from Wagner 2004)';

% CSF
S(3)=sim_struct('COND');
S(3).name = 'CSF';
S(3).value = 1.654;
S(3).descrip = 'cerebrospinal fluid (from Wagner 2004)';

% Bone
S(4)=sim_struct('COND');
S(4).name = 'Bone';
S(4).value = 0.010;
S(4).descrip = 'average bone (from Wagner 2004)';

% Scalp
S(5)=sim_struct('COND');
S(5).name = 'Scalp';
S(5).value = 0.465;
S(5).descrip = 'average scalp (from Wagner 2004)';

% Eye balls (vitreous humour)
S(6)=sim_struct('COND');
S(6).name = 'Eye_balls';
S(6).value = 0.5;
S(6).descrip = 'vitreous humour (from Opitz, ..., Thielscher, NI, 2015)';

% Compact bone
S(7)=sim_struct('COND');
S(7).name = 'Compact_bone';
S(7).value = 0.008;
S(7).descrip = 'compact bone (from Opitz, ..., Thielscher, NI, 2015)';

% Spongy bone
S(8)=sim_struct('COND');
S(8).name = 'Spongy_bone';
S(8).value = 0.025;
S(8).descrip = 'spongy bone (from Opitz, ..., Thielscher, NI, 2015)';

% Blood
S(9)=sim_struct('COND');
S(9).name = 'Blood';
S(9).value = 0.6;
S(9).descrip = 'Blood (from Gabriel et al., 2009)';

% Muscle
S(10)=sim_struct('COND');
S(10).name = 'Muscle';
S(10).value = 0.16;
S(10).descrip = 'Muscle (from Gabriel et al., 2009)';

% Rubber
S(100)=sim_struct('COND');
S(100).name = 'Electrode_rubber';
S(100).value = 29.4;
S(100).descrip = 'for tDCS rubber electrodes (for neuroConn rubber: Wacker Elastosil R 570/60 RUSS)';

% Saline
S(500)=sim_struct('COND');
S(500).name = 'Saline';
S(500).value = 1.0;
S(500).descrip = 'for tDCS sponge electrodes';

function S = sim_struct(str)
S.name = '';
S.value = '';
S.descrip = '';