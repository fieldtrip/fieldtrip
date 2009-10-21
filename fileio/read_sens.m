function [sens] = read_sens(filename, varargin)

% READ_SENS read sensor positions from various manufacturer specific files. 
% Currently supported are ASA, BESA, Polhemus and Matlab for EEG 
% electrodes and CTF and Neuromag for MEG gradiometers.
%
% Use as
%   grad = read_sens(filename, ...)  % for gradiometers
%   elec = read_sens(filename, ...)  % for electrodes
%
% Additional options should be specified in key-value pairs and can be
%   'fileformat'   string
%
% An electrode definition contain the following fields
%   elec.pnt     Nx3 matrix with carthesian (x,y,z) coordinates of each electrodes
%   elec.label   cell-array of length N with the label of each electrode
%
% A gradiometer definition generally consists of multiple coils per
% channel, e.g.two coils for a 1st order gradiometer in which the
% orientation of the coils is opposite. Each coil is described
% separately and a large "tra" matrix (can be sparse) has to be
% given that defines how the forward computed field is combined over
% the coils to generate the output of each channel. The gradiometer
% definition constsis of the following fields
%   grad.pnt     Mx3 matrix with the position of each coil
%   grad.ori     Mx3 matrix with the orientation of each coil
%   grad.tra     NxM matrix with the weight of each coil into each channel
%   grad.label   cell-array of length N with the label of each of the channels
%
% See also TRANSFORM_SENS, PREPARE_VOL_SENS, COMPUTE_LEADFIELD

% Copyright (C) 2005-2008, Robert Oostenveld
%
% $Log: read_sens.m,v $
% Revision 1.16  2009/10/16 07:32:08  roboos
% renamed chieti into itab for consistency with other formats
%
% Revision 1.15  2009/10/13 10:12:51  roboos
% added support for chieti_raw
%
% Revision 1.14  2009/07/02 10:34:21  vlalit
% Added eeglab_set to the list of formats where electrode locations can be found in
%  the header.
%
% Revision 1.13  2009/06/03 09:52:15  roboos
% added zebris_sfp
%
% Revision 1.12  2009/02/02 16:10:15  vlalit
% Provide the 'headertype' argument to the internal read_header call.
%
% Revision 1.11  2009/01/23 10:32:55  vlalit
% New reader for Neuromag fif format using the MNE toolbox (http://www.nmr.mgh.harvard.edu/martinos/userInfo/data/sofMNE.php)  implemented by Laurence Hunt.
%
% Revision 1.10  2008/09/18 10:38:33  vlalit
% Added 4D formats to be recognized by read_sens
%
% Revision 1.9  2008/05/22 14:33:18  vlalit
% Changes related to generalization of fiducials'  handling in SPM.
%
% Revision 1.8  2008/04/14 20:51:36  roboos
% added convert_units
%
% Revision 1.7  2008/04/11 16:17:22  roboos
% added polhemus_fil
%
% Revision 1.6  2008/03/20 13:43:14  roboos
% added support for besa_pos
%
% Revision 1.5  2008/03/18 12:34:30  roboos
% fixed bug: added varargin to input arguments, thanks to Juan
%
% Revision 1.4  2008/03/06 09:27:54  roboos
% updated documentation
%
% Revision 1.3  2008/03/05 11:06:11  roboos
% test the presence of the fileio toolbox, needed when this function is included in forwinv
%
% Revision 1.2  2008/03/05 10:54:05  roboos
% added optional argument for fileformat
% some documentation changes
%
% Revision 1.1  2008/01/28 20:10:11  roboos
% new functions based on existing fieldtrip code
%

% test whether the file exists
if ~exist(filename)
  error(sprintf('file ''%s'' does not exist', filename));
end

% get the options
fileformat = keyval('fileformat',  varargin);

% determine the filetype
if isempty(fileformat)
  fileformat = filetype(filename);
end

switch fileformat

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % read the content from various files that contain EEG electrode positions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'asa_elc'
    sens = read_asa_elc(filename);
  
  case 'polhemus_pos'
    sens = read_brainvision_pos(filename);

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
      warning(sprintf('reading channel labels from %s', elpfile));
      % read the channel names from the accompanying ELP file
      lbl = importdata(elpfile);
      sens.label = strrep(lbl.textdata(:,2) ,'''', '');
    elseif exist(elafile, 'file')
      warning(sprintf('reading channel labels from %s', elafile));
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
    sens.label = tmp{1};
    sens.pnt   = [tmp{2:4}];
   
  case 'itab_raw'
    hastoolbox('fileio');
    hdr = read_header(filename);
    sens = hdr.grad;
    
  case 'neuromag_mne'
    hastoolbox('fileio');
    hdr = read_header(filename,'headerformat','neuromag_mne');
    sens = hdr.elec;    
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % gradiometer information is always stored in the header of the MEG dataset
  % hence uses the standard fieldtrip/fileio read_header function
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  case {'ctf_ds', 'ctf_res4', 'neuromag_fif', '4d', '4d_pdf', '4d_m4d', '4d_xyz', 'yokogawa_ave', 'yokogawa_con', 'yokogawa_raw'}
    % check the availability of the required low-level toolbox
    % this is required because the read_sens function is also on itself included in the forwinv toolbox
    hastoolbox('fileio');
    hdr = read_header(filename, 'headerformat', fileformat);
    sens = hdr.grad;
    
    
  case 'neuromag_mne_grad'
    hastoolbox('fileio');
    hdr = read_header(filename,'headerformat','neuromag_mne');
    sens = hdr.grad;
    
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % This is for EEG formats where electrode positions can be stored with the data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  case {'spmeeg_mat', 'eeglab_set'}
    % check the availability of the required low-level toolbox
    % this is required because the read_sens function is also on itself included in the forwinv toolbox
    hastoolbox('fileio');
    hdr = read_header(filename);
    
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
    warning('off', 'MATLAB:load:variableNotFound');
    tmp = load(matfile, 'elec', 'grad', 'sens', 'elc');
    warning('on', 'MATLAB:load:variableNotFound');
    if isfield(tmp, 'grad')
      sens = getfield(tmp, 'grad');
    elseif isfield(tmp, 'elec')
      sens = getfield(tmp, 'elec');
    elseif isfield(tmp, 'sens')
      sens = getfield(tmp, 'sens');
    elseif isfield(tmp, 'elc')
      sens = getfield(tmp, 'elc');
    else
      error('no electrodes or gradiometers found in Matlab file');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % these are created by a Zebris tracker, at CRC in Liege at least.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  case 'zebris_sfp'
    [sens.fid, sens.pnt, sens.fid_label, sens.label] = read_zebris(filename, 0);

  otherwise
    error('unknown fileformat for electrodes or gradiometers');
end

if senstype(sens, 'eeg')
  % only keep positions and labels in case of EEG electrodes
  dum  = sens;
  sens = [];
  sens.pnt   = dum.pnt;
  sens.label = dum.label;
end

% this will add the units to the sensor array
sens = convert_units(sens);
