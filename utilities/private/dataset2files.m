function [filename, headerfile, datafile] = dataset2files(filename, format)

% DATASET2FILES manages the filenames for the dataset, headerfile, datafile and eventfile
% and tries to maintain a consistent mapping between them for each of the known fileformats
%
% Use as
%   [filename, headerfile, datafile] = dataset2files(filename, format)

% Copyright (C) 2007-2019, Robert Oostenveld
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

persistent previous_argin previous_argout

if iscell(filename)
  % use recursion to go over multiple files
  headerfile = cell(size(filename));
  datafile   = cell(size(filename));
  for i=1:numel(filename)
    [filename{i}, headerfile{i}, datafile{i}] = dataset2files(filename{i}, format);
  end
  return
end

if isempty(format)
  format = ft_filetype(filename);
end

current_argin = {filename, format};
if isequal(current_argin, previous_argin)
  % don't do the whole cheking again, but return the previous output from cache
  filename   = previous_argout{1};
  headerfile = previous_argout{2};
  datafile   = previous_argout{3};
  return
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
  case {'anywave_ades', 'anywave_dat'}
    [path, file, ext] = fileparts(filename);
    datafile   = fullfile(path, [file '.dat']);
    headerfile = fullfile(path, [file '.ades']);
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
  case {'ctf_meg4' 'ctf_res4' 'ctf_read_meg4' 'ctf_read_res4' 'read_ctf_meg4' 'read_ctf_res4'}
    [path, file, ext] = fileparts(filename);
    if strcmp(ext, '.ds')
      % the directory name was specified instead of the meg4/res4 file
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
    if exist(fullfile(path, [file '.eeg']), 'file')
      datafile   = fullfile(path, [file '.eeg']);
    elseif exist(fullfile(path, [file '.seg']), 'file')
      datafile   = fullfile(path, [file '.seg']);
    elseif exist(fullfile(path, [file '.dat']), 'file')
      datafile   = fullfile(path, [file '.dat']);
    else
      ft_error('cannot determine the data file that corresponds to %s', filename);
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
    if isfolder(filename)
      path = filename;
    else
      [path, file, ext] = fileparts(filename);
    end
    headerfile = fullfile(path, 'header');
    datafile   = fullfile(path, 'samples');
  case {'tdt_tsq' 'tdt_tev'}
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.tsq']);
    datafile   = fullfile(path, [file '.tev']);
  case 'egi_mff'
    if ~isfolder(filename);
      [path, file, ext] = fileparts(filename);
      headerfile = path;
      datafile   = path;
    else
      headerfile = filename;
      datafile   = filename;
    end
  case {'deymed_dat' 'deymed_ini'}
    [p, f, x] = fileparts(filename);
    headerfile = fullfile(p, [f '.ini']);
    if ~exist(headerfile, 'file')
      headerfile = fullfile(p, [f '.Ini']);
    end
    datafile = fullfile(p, [f '.dat']);
    if ~exist(datafile, 'file')
      datafile = fullfile(p, [f '.Dat']);
    end
  case 'neurosim_ds'
    % this is the directory
    filename = fullfile(filename, 'signals'); % this is the only one we care about for the continuous signals
    headerfile = filename;
    datafile   = filename;
  case 'tmsi_poly5'
    [p, f, x] = fileparts(filename);
    if strcmpi(x, '.poly5')
      headerfile = filename;
      datafile = filename;
    else
      filename = fullfile(p, f, [f '.eeg.poly5']);
      if ~exist(filename , 'file')
        filename  = fullfile(p, f, [f '.EEG.Poly5']);
      end
      headerfile = filename;
      datafile = filename;
    end
  otherwise
    % convert filename into filenames, assume that the header and data are the same
    datafile   = filename;
    headerfile = filename;
end

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
current_argout = {filename, headerfile, datafile};
previous_argin  = current_argin;
previous_argout = current_argout;
