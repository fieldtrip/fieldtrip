function [headmodel] = ft_read_vol(filename, varargin)

% FT_READ_VOL reads a volume conduction model from various manufacturer
% specific files. Currently supported are ASA, CTF, Neuromag, MBFYS
% and Matlab.
%
% Use as
%   headmodel = ft_read_vol(filename, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'fileformat'   string
%
% The volume conduction model is represented as a structure with fields
% that depend on the type of model.
%
% See also FT_TRANSFORM_VOL, FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

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
if ~exist(filename)
  ft_error('file ''%s'' does not exist', filename);
end

% get the options
fileformat = ft_getopt(varargin, 'fileformat', ft_filetype(filename));

switch fileformat
  case 'matlab'
    % FIXME in the future the file should contain the variable 'headmodel' instead of vol
    headmodel = loadvar(filename, 'vol');
    
  case 'ctf_hdm'
    headmodel = read_ctf_hdm(filename);
    
  case 'asa_vol'
    headmodel = read_asa_vol(filename);
    headmodel.type = 'asa';
    
  case 'mbfys_ama'
    ama = loadama(filename);
    headmodel = ama2vol(ama);
    
  case 'neuromag_fif'
    ft_hastoolbox('mne', 1);
    global FIFF
    bem = mne_read_bem_surfaces(filename);
    headmodel.bnd.pos = bem.rr;
    headmodel.bnd.tri = bem.tris;
    headmodel.coordsys = fif2coordsys(bem.coord_frame);
    
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
