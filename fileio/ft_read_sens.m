function [sens] = ft_read_sens(filename, varargin)

% FT_READ_SENS read sensor positions from various manufacturer specific files. The
% following acquisition system and analysis platform file formats are currently
% supported:
%
%   asa_elc besa_elp besa_pos besa_sfp yokogawa_ave yokogawa_con yokogawa_raw 4d
%   4d_pdf 4d_m4d 4d_xyz ctf_ds ctf_res4 itab_raw itab_mhd netmeg neuromag_fif
%   neuromag_mne neuromag_mne_elec neuromag_mne_grad polhemus_fil polhemus_pos
%   zebris_sfp spmeeg_mat eeglab_set localite_pos matlab
%
% Use as
%   grad = ft_read_sens(filename, ...)  % for gradiometers
%   elec = ft_read_sens(filename, ...)  % for electrodes
%
% Additional options should be specified in key-value pairs and can be
%   'fileformat'     = string, see the list of supported file formats (the default is determined automatically)
%   'senstype'       = string, can be 'eeg' or 'meg', specifies which type of sensors to read from a fif file (default = 'eeg')
%   'coordsys'       = string, 'head' or 'dewar' (default = 'head')
%   'coilaccuracy'   = can be empty or a number (0, 1 or 2) to specify the accuracy (default = [])
%
% An electrode definition contain the following fields
%   elec.elecpos = Nx3 matrix with carthesian (x,y,z) coordinates of each
%                  electrode
%   elec.label   = cell-array of length N with the label of each electrode
%   elec.chanpos = Nx3 matrix with coordinates of each sensor
%
% A gradiometer definition generally consists of multiple coils per channel, e.g. two
% coils for a 1st order gradiometer in which the orientation of the coils is
% opposite. Each coil is described separately and a large "tra" matrix has to be
% given that defines how the forward computed field is combined over the coils to
% generate the output of each channel. The gradiometer definition constsis of the
% following fields
%   grad.coilpos = Mx3 matrix with the position of each coil
%   grad.coilori = Mx3 matrix with the orientation of each coil
%   grad.tra     = NxM matrix with the weight of each coil into each channel
%   grad.label   = cell-array of length N with the label of each of the channels
%   grad.chanpos = Nx3 matrix with the positions of each sensor
%
% See also FT_READ_HEADER, FT_TRANSFORM_SENS, FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD,
% FT_DATATYPE_SENS

% Copyright (C) 2005-2016 Robert Oostenveld
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
  error('file ''%s'' does not exist', filename);
end

% get the options
fileformat     = ft_getopt(varargin, 'fileformat', ft_filetype(filename));
senstype       = ft_getopt(varargin, 'senstype', 'eeg');  % can be eeg or meg, this is used to decide what to return if both are present in a fif file
coordsys       = ft_getopt(varargin, 'coordsys', 'head'); % this is used for ctf and neuromag_mne, it can be head or dewar
coilaccuracy   = ft_getopt(varargin, 'coilaccuracy');     % empty, or a number between 0 to 2

switch fileformat
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % read the content from various files that contain EEG electrode positions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'asa_elc'
    sens = read_asa_elc(filename);
    
  case 'polhemus_pos'
    sens = read_brainvision_pos(filename);
    
  case 'besa_elp'
    error('unknown fileformat for electrodes or gradiometers');
    % the code below does not yet work
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
    
  case 'besa_sfh'
    sfh = readBESAsfh(filename);
    sens.label   = sfh.SurfacePointsLabels(:);
    sens.elecpos = sfh.SurfacePointsCoordinates(:,1:3);
    sel = true(sfh.NrOfPoints, 1);
    for i=1:sfh.NrOfPoints
      tok = tokenize(sens.label{i}, '_');
      sens.label{i} = tok{2};
      sel(i) = ~strcmp(tok{1}, 'Fid');
    end
    sens.label   = sens.label(sel);
    sens.elecpos = sens.elecpos(sel,:);
    
  case 'besa_sfp'
    [lab, pos] = read_besa_sfp(filename);
	
    sens.label   = lab;
    sens.elecpos = pos;
		
  case {'ctf_ds', 'ctf_res4', 'ctf_old', 'neuromag_fif', 'neuromag_mne', '4d', '4d_pdf', '4d_m4d', '4d_xyz', 'yokogawa_ave', 'yokogawa_con', 'yokogawa_raw', 'itab_raw' 'itab_mhd', 'netmeg'}
    % gradiometer information is always stored in the header of the MEG dataset, hence uses the standard fieldtrip/fileio ft_read_header function
    hdr = ft_read_header(filename, 'headerformat', fileformat, 'coordsys', coordsys, 'coilaccuracy', coilaccuracy);
    % sometimes there can also be electrode position information in the header
    if isfield(hdr, 'elec') && isfield(hdr, 'grad')
      switch lower(senstype)
        case 'eeg'
          warning('both electrode and gradiometer information is present, returning the electrode information');
          sens = hdr.elec;
        case 'meg'
          warning('both electrode and gradiometer information is present, returning the gradiometer information');
          sens = hdr.grad;
      end
    elseif ~isfield(hdr, 'elec') && ~isfield(hdr, 'grad')
      error('neither electrode nor gradiometer information is present');
    elseif isfield(hdr, 'grad')
      sens = hdr.grad;
    elseif isfield(hdr, 'elec')
      sens = hdr.elec;
    end
    
  case 'neuromag_mne_grad'
    % the file can contain both, force reading the gradiometer info
    % note that this functionality overlaps with senstype=eeg/meg
    hdr = ft_read_header(filename, 'headerformat', 'neuromag_mne', 'coordsys', coordsys, 'coilaccuracy', coilaccuracy);
    sens = hdr.grad;
    
  case 'neuromag_mne_elec'
    % the file can contain both, force reading the electrode info
    % note that this functionality overlaps with senstype=eeg/meg
    hdr = ft_read_header(filename, 'headerformat', 'neuromag_mne', 'coordsys', coordsys, 'coilaccuracy', coilaccuracy);
    sens = hdr.elec;
    
  case {'spmeeg_mat', 'eeglab_set'}
    % this is for EEG formats where electrode positions can be stored with the data
    hdr = ft_read_header(filename, 'coordsys', coordsys, 'coilaccuracy', coilaccuracy);
    if isfield(hdr, 'grad')
      sens = hdr.grad;
    elseif isfield(hdr, 'elec')
      sens = hdr.elec;
    else
      error('no electrodes or gradiometers found in the file')
    end
    
  case 'polhemus_fil'
    % these are created at the FIL in London with a polhemus tracker
    [sens.fid.pnt, sens.pnt, sens.fid.label] = read_polhemus_fil(filename, 0);
    % the file does not have channel labels in it
    warning('no channel names in polhemus file, using numbers instead');
    for i=1:size(sens.pnt, 1)
      sens.label{i} = sprintf('%03d', i);
    end
    
  case 'matlab'
    % MATLAB files can contain either electrodes or gradiometers
    matfile = filename;   % this solves a problem with the MATLAB compiler v3
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
      error('no electrodes or gradiometers found in MATLAB file');
    end
    
  case 'zebris_sfp'
    % these are created by a Zebris tracker, at CRC in Liege at least.
    [sens.fid.pnt, sens.chanpos, sens.fid.label, sens.label] = read_zebris(filename, 0);
    % convert to columns
    sens.label = sens.label(:);
    sens.fid.label = sens.fid.label(:);
    
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
    if false
      % return all positions, including the ones that do not correspond to
      % electrodes per see, such as the fiducials and localizer coils
      sens          = [];
      sens.label    = l;
      sens.elecpos  = [x y z];
    else
      % split the electrodes and fiducials
      % this is consistent with zebris_sfp and with the output of ft_read_headshape
      sens            = [];
      sens.label      = l(~sel);
      sens.elecpos    = [x(~sel) y(~sel) z(~sel)];
      sens.fid.label  = l(sel);
      sens.fid.pnt    = [x(sel) y(sel) z(sel)];
    end
    
  case {'localite_pos','localite_ins'}
    if ~usejava('jvm') % Using xml2struct requires java
      fid = fopen(filename);
      
      % Read marker-file and store contents in cells of strings
      tmp = textscan(fid,'%s');
      
      fclose(fid);
      
      % Search for cells that contain coordinates
      selx = strncmp('data0',tmp{1},5);
      sely = strncmp('data1',tmp{1},5);
      selz = strncmp('data2',tmp{1},5);
      sellab = strncmp('description',tmp{1},5);
      
      % Extract cells that contain coordinates
      xtemp  = tmp{1}(selx);
      ytemp  = tmp{1}(sely);
      ztemp  = tmp{1}(selz);
      labtemp = tmp{1}(sellab);
      
      % Determine which channels are set. In localite channels that are not set
      % automatically receive coordinates [0, 0, 0] and should therefore
      % be discarded.
      settemp = tmp{1}(strncmp('set',tmp{1},3));
      selset = strncmp('set="f',settemp,6);
      
      % Remove channels that are not set
      xtemp(selset) = [];
      ytemp(selset) = [];
      ztemp(selset) = [];
      labtemp(selset) = [];
      
      % Convert cells that contain coordinates from string to double
      x = [];
      y = [];
      z = [];
      lbl = [];
      
      for i=1:numel(xtemp)
        x(i,1) = str2double(xtemp{i}(8:end-1));
        y(i,1) = str2double(ytemp{i}(8:end-1));
        z(i,1) = str2double(ztemp{i}(8:end-3));
        lbl{i,1} = labtemp{i}(14:end-1);
      end;
      
      % Create and fill sens structure
      sens = [];
      sens.elecpos = [x y z];
      sens.chanpos = sens.elecpos;
      sens.label = lbl;
    else
      tmp = xml2struct(filename);
      
      sens = [];
      
      % Loop through structure obtained from xml-file and store
      % coordinate information into sens structure of channels that are
      % set.
      for i=1:numel(tmp)
        if strcmp(tmp(i).Marker.set,'true')
          sens.elecpos(i,1) = str2double(tmp(i).Marker.ColVec3D.data0);
          sens.elecpos(i,2) = str2double(tmp(i).Marker.ColVec3D.data1);
          sens.elecpos(i,3) = str2double(tmp(i).Marker.ColVec3D.data2);
          sens.label{i} = tmp(i).Marker.description;
        end;
      end;
      
      sens.chanpos = sens.elecpos;
    end;
    
  case 'easycap_txt'
    % Read the file and store all contents in cells of strings
    fid = fopen(filename);
    tmp = textscan(fid,'%s%s%s%s');
    fclose(fid);
    
    sens = [];
    if all(cellfun(@isempty, tmp{4}))
      % it contains theta and phi
      sens.label = tmp{1}(2:end);
      theta = cellfun(@str2double, tmp{2}(2:end));
      phi   = cellfun(@str2double, tmp{3}(2:end));
      radians = @(x) pi*x/180;
      warning('assuming a head radius of 85 mm');
      x = 85*cos(radians(phi)).*sin(radians(theta));
      y = 85*sin(radians(theta)).*sin(radians(phi));
      z = 85*cos(radians(theta));
      sens.unit = 'cm';
      sens.elecpos = [x y z];
      sens.chanpos = [x y z];
    else
      % it contains X, Y, Z
      sens.label = tmp{1}(2:end);
      x = cellfun(@str2double, tmp{2}(2:end));
      y = cellfun(@str2double, tmp{3}(2:end));
      z = cellfun(@str2double, tmp{4}(2:end));
      sens.elecpos = [x y z];
      sens.chanpos = [x y z];
    end
    
  otherwise
    error('unknown fileformat for electrodes or gradiometers');
end

% ensure that the sensor description is up-to-date
% this will also add chantype and units to the sensor array if missing
sens = ft_datatype_sens(sens);
