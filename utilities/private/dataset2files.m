function [filename, headerfile, datafile] = dataset2files(filename, format)

% DATASET2FILES manages the filenames for the dataset, headerfile, datafile and eventfile
% and tries to maintain a consistent mapping between them for each of the known fileformats
%
% Use as
%   [filename, headerfile, datafile] = dataset2files(filename, format)

% Copyright (C) 2007-2011, Robert Oostenveld
%
% $Log$

if isempty(format)
  format = ft_filetype(filename);
end

switch format
  case '4d_pdf'
    datafile   = filename;
    headerfile = [datafile '.m4d'];
    sensorfile = [datafile '.xyz'];
  case {'4d_m4d', '4d_xyz'}
    datafile   = filename(1:(end-4)); % remove the extension
    headerfile = [datafile '.m4d'];
    sensorfile = [datafile '.xyz'];
  case '4d'
    [path, file, ext] = fileparts(filename);
    datafile   = fullfile(path, [file,ext]);
    headerfile = fullfile(path, [file,ext]);
    configfile = fullfile(path, 'config');
  case {'ctf_ds', 'ctf_old'}
    % convert CTF filename into filenames
    [path, file, ext] = fileparts(filename);
    if any(strcmp(ext, {'.res4' '.meg4', '.1_meg4' '.2_meg4' '.3_meg4' '.4_meg4' '.5_meg4' '.6_meg4' '.7_meg4' '.8_meg4' '.9_meg4'}))
      filename = path;
      [path, file, ext] = fileparts(filename);
    end
    if isempty(path) && isempty(file)
      % this means that the dataset was specified as the present working directory, i.e. only with '.'
      filename = pwd;
      [path, file, ext] = fileparts(filename);
    end
    headerfile = fullfile(filename, [file '.res4']);
    datafile   = fullfile(filename, [file '.meg4']);
    if length(path)>3 && strcmp(path(end-2:end), '.ds')
      filename = path; % this is the *.ds directory
    end
  case 'ctf_meg4'
    [path, file, ext] = fileparts(filename);
    if strcmp(ext, '.ds')
      % the directory name was specified instead of the meg4 file
      path = filename;
    end
    if isempty(path)
      path = pwd;
    end
    headerfile = fullfile(path, [file '.res4']);
    datafile   = fullfile(path, [file '.meg4']);
    if length(path)>3 && strcmp(path(end-2:end), '.ds')
      filename = path; % this is the *.ds directory
    end
  case 'ctf_res4'
    [path, file, ext] = fileparts(filename);
    if strcmp(ext, '.ds')
      % the directory name was specified instead of the meg4 file
      path = filename;
    end
    if isempty(path)
      path = pwd;
    end
    headerfile = fullfile(path, [file '.res4']);
    datafile   = fullfile(path, [file '.meg4']);
    if length(path)>3 && strcmp(path(end-2:end), '.ds')
      filename = path; % this is the *.ds directory
    end
  case 'brainvision_vhdr'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.vhdr']);
    if exist(fullfile(path, [file '.eeg']))
      datafile   = fullfile(path, [file '.eeg']);
    elseif exist(fullfile(path, [file '.seg']))
      datafile   = fullfile(path, [file '.seg']);
    elseif exist(fullfile(path, [file '.dat']))
      datafile   = fullfile(path, [file '.dat']);
    end
  case 'brainvision_eeg'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.vhdr']);
    datafile   = fullfile(path, [file '.eeg']);
  case 'brainvision_seg'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.vhdr']);
    datafile   = fullfile(path, [file '.seg']);
  case 'brainvision_dat'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.vhdr']);
    datafile   = fullfile(path, [file '.dat']);
  case 'itab_raw'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.raw.mhd']);
    datafile   = fullfile(path, [file '.raw']);
  case 'fcdc_matbin'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.mat']);
    datafile   = fullfile(path, [file '.bin']);
  case 'fcdc_buffer_offline'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '/header']);
  case {'tdt_tsq' 'tdt_tev'}
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.tsq']);
    datafile   = fullfile(path, [file '.tev']);
  otherwise
    % convert filename into filenames, assume that the header and data are the same
    datafile   = filename;
    headerfile = filename;
end
