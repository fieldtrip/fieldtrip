function [sens] = ft_read_sens(filename, varargin)

% FT_READ_SENS read sensor positions from various manufacturer
% specific files. The following acsuisition system and analysis
% platform EEG and MEG file formats are currently supported:
%
%   asa_elc besa_elp besa_pos besa_sfp yokogawa_ave yokogawa_con
%   yokogawa_raw 4d 4d_pdf 4d_m4d 4d_xyz ctf_ds ctf_res4 itab_raw
%   itab_mhd netmeg neuromag_fif neuromag_mne neuromag_mne_elec
%   neuromag_mne_grad polhemus_fil polhemus_pos zebris_sfp spmeeg_mat
%   eeglab_set matlab
%
% Use as
%   grad = ft_read_sens(filename, ...)  % for gradiometers
%   elec = ft_read_sens(filename, ...)  % for electrodes
%
% Additional options should be specified in key-value pairs and can be
%   'fileformat' = string, see the list of supported file formats (the
%                  default is determined automatically)
%
% An electrode definition contain the following fields
%   elec.elecpos = Nx3 matrix with carthesian (x,y,z) coordinates of each
%                  electrode
%   elec.label   = cell-array of length N with the label of each electrode
%   elec.chanpos = Nx3 matrix with coordinates of each sensor
%
% A gradiometer definition generally consists of multiple coils per
% channel, e.g.two coils for a 1st order gradiometer in which the
% orientation of the coils is opposite. Each coil is described
% separately and a large "tra" matrix has to be given that defines
% how the forward computed field is combined over the coils to generate
% the output of each channel. The gradiometer definition constsis of
% the following fields
%   grad.coilpos = Mx3 matrix with the position of each coil
%   grad.coilori = Mx3 matrix with the orientation of each coil
%   grad.tra     = NxM matrix with the weight of each coil into each channel
%   grad.label   = cell-array of length N with the label of each of the channels
%   grad.chanpos = Nx3 matrix with the positions of each sensor
%
% See also FT_TRANSFORM_SENS, FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD,
% FT_DATATYPE_SENS

% Copyright (C) 2005-2010 Robert Oostenveld
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

% optionally get the data from the URL and make a temporary local copy
filename = fetch_url(filename);

% test whether the file exists
if ~exist(filename, 'file')
  error('file ''%s'' does not exist', filename);
end

% get the options
fileformat = ft_getopt(varargin, 'fileformat', ft_filetype(filename));

switch fileformat
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % read the content from various files that contain EEG electrode positions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'asa_elc'
    sens = read_asa_elc(filename);
    
  case 'polhemus_pos'
    sens = read_brainvision_pos(filename);
    
  case 'besa_elp'
    % the code below does not yet work
    error('unknown fileformat for electrodes or gradiometers');
    
    fid = fopen(filename);
    % the ascii file contains: type, label, angle, angle
    tmp = textscan(fid, '%s%s%f%f');
    fclose(fid);
    
    sel = strcmpi(tmp{1}, 'EEG');  % type can be EEG or POS
    sens.label = tmp{2}(sel);
    az = tmp{3}(sel) * pi/180;
    el = tmp{4}(sel) * pi/180;
    r  = ones(size(el));
    [x, y, z] = sph2cart(az, el, r);
    sens.chanpos = [x y z];
    
  case 'besa_pos'
    tmp = importdata(filename);
    if ~isnumeric(tmp)
      error('unexpected file format for fileformat=besa_pos')
    end
    [nchan,nrow] = size(tmp);
    if nrow==3
      sens.pnt = tmp;
    elseif nrow==9
      pnt1 = tmp(:,1:3);  % bottom coil
      pnt2 = tmp(:,4:6);  % top coil
      ori  = tmp(:,7:9);  % orientation of bottom coil
      sens.pnt = [pnt1; pnt2];
      sens.ori = [ori; ori];
      sens.tra = [eye(nchan) -eye(nchan)];
    else
      error('unexpected file format for fileformat=besa_pos')
    end
    [p, f, x] = fileparts(filename);
    elpfile = fullfile(p, [f '.elp']);
    elafile = fullfile(p, [f '.ela']);
    if exist(elpfile, 'file')
      warning('reading channel labels from %s', elpfile);
      % read the channel names from the accompanying ELP file
      lbl = importdata(elpfile);
      sens.label = strrep(lbl.textdata(:,2) ,'''', '');
    elseif exist(elafile, 'file')
      warning('reading channel labels from %s', elafile);
      % read the channel names from the accompanying ELA file
      lbl = importdata(elafile);
      lbl = strrep(lbl, 'MEG ', ''); % remove the channel type
      lbl = strrep(lbl, 'EEG ', ''); % remove the channel type
      sens.label = lbl;
    else
      % the file does not have channel labels in it
      warning('creating fake channel names for besa_pos');
      for i=1:nchan
        sens.label{i} = sprintf('%03d', i);
      end
    end
    
  case 'besa_sfp'
    fid        = fopen(filename);
    tmp        = textscan(fid, ' %[^ \t]%n%n%n');
    fclose(fid);
    sens.label   = tmp{1};
    sens.chanpos = [tmp{2:4}];
    sens.elecpos = sens.chanpos;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % gradiometer information is always stored in the header of the MEG dataset
    % hence uses the standard fieldtrip/fileio read_header function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  case {'ctf_ds', 'ctf_res4', 'ctf_old', 'neuromag_fif', '4d', '4d_pdf', '4d_m4d', '4d_xyz', 'yokogawa_ave', 'yokogawa_con', 'yokogawa_raw', 'itab_raw' 'itab_mhd', 'netmeg'}
    hdr = ft_read_header(filename, 'headerformat', fileformat);
    sens = hdr.grad;
    
  case 'neuromag_mne_grad'
    % the file can contain both, force reading the gradiometer info
    hdr = ft_read_header(filename,'headerformat','neuromag_mne');
    sens = hdr.grad;
    
  case 'neuromag_mne_elec'
    % the file can contain both, force reading the electrode info
    hdr = ft_read_header(filename,'headerformat','neuromag_mne');
    sens = hdr.elec;
    
  case 'neuromag_mne'
    % the file can contain both, try to be smart in determining what to return
    hdr = ft_read_header(filename,'headerformat','neuromag_mne');
    if isfield(hdr, 'elec') && isfield(hdr, 'grad')
      warning('returning electrode information, not gradiometer location');
      sens = hdr.elec;
    elseif isfield(hdr, 'elec')
      sens = hdr.elec;
    elseif isfield(hdr, 'grad')
      sens = hdr.grad;
    else
      error('cannot find electrode or gradiometer information');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is for EEG formats where electrode positions can be stored with the data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  case {'spmeeg_mat', 'eeglab_set'}
    hdr = ft_read_header(filename);
    
    if isfield(hdr, 'grad')
      sens = hdr.grad;
    elseif isfield(hdr, 'elec')
      sens = hdr.elec;
    else
      error('no electrodes or gradiometers found in the file')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % these are created at the FIL in London with a polhemus tracker
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  case 'polhemus_fil'
    [sens.fid, sens.pnt] = read_polhemus_fil(filename, 0);
    
    % the file does not have channel labels in it
    warning('no channel names in polhemus file, using numbers instead');
    for i=1:size(sens.pnt, 1)
      sens.label{i} = sprintf('%03d', i);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % matlab files can contain either electrodes or gradiometers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  case 'matlab'
    matfile = filename;   % this solves a problem with the matlab compiler v3
    ws = warning('off', 'MATLAB:load:variableNotFound');
    tmp = load(matfile, 'elec', 'grad', 'sens', 'elc');
    warning(ws);
    if isfield(tmp, 'grad')
      sens = tmp.grad;
    elseif isfield(tmp, 'elec')
      sens = tmp.elec;
    elseif isfield(tmp, 'sens')
      sens = tmp.sens;
    elseif isfield(tmp, 'elc')
      sens = tmp.elc;
    else
      error('no electrodes or gradiometers found in Matlab file');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % these are created by a Zebris tracker, at CRC in Liege at least.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  case 'zebris_sfp'
    [sens.fid, sens.chanpos, sens.fid_label, sens.label] = read_zebris(filename, 0);
    
  case '4d_el_ascii'
    fid = fopen(filename, 'rt');
    c = textscan(fid, '%s%s%f%f%f');
    l = c{:,1}; % label
    s = c{:,2}; % status, it can be 'Collected' or empty
    x = c{:,3};
    y = c{:,4};
    z = c{:,5};
    % shift the columns with one where the status is not specified
    sel = isnan(z);
    z(sel) = y(sel);
    y(sel) = x(sel);
    x(sel) = str2double(s(sel));
    s(sel) = {''};
    fclose(fid);
    % return all positions, including the ones that do not correspond to
    % electrodes per see, such as the fiducials and localizer coils
    sens          = [];
    sens.label    = l;
    sens.elecpos  = [x y z];
    
  otherwise
    error('unknown fileformat for electrodes or gradiometers');
end

% ensure that the sensor description is up-to-date
% this will also add the units to the sensor array if missing
sens = ft_datatype_sens(sens);

if ft_senstype(sens, 'eeg')
  % only keep positions and labels in case of EEG electrodes
  % FIXME what is removed here?
  dum  = sens;
  sens = [];
  sens.chanpos = dum.chanpos;
  sens.elecpos = dum.elecpos;
  sens.label   = dum.label;
end
