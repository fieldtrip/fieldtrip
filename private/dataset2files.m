function cfg = dataset2files(cfg)

% Helper function that converts a dataset into headerfile and datafile
% if necessary. This is used in PREPROCESSING and DEFINETRIAL
%
% This function operates only on
%   cfg.dataset
%   cfg.datafile
%   cfg.headerfile
% and returns the updated cfg.

% Copyright (C) 2004, Robert Oostenveld
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

% start with empty fields if thery are not present
if ~isfield(cfg, 'dataset')
  cfg.dataset = [];
end
if ~isfield(cfg, 'datafile')
  cfg.datafile = [];
end
if ~isfield(cfg, 'headerfile')
  cfg.headerfile = [];
end

if ~isempty(cfg.dataset)
  if strcmp(cfg.dataset, 'gui')
    d = uigetdir;
    if d==0
      [f, p] = uigetfile;
      if f==0
        ft_error('You should select a dataset file or directory');
      else
        d = fullfile(p, f);
      end
    end
    cfg.dataset = d;
  end
  
  switch ft_filetype(cfg.dataset)
    case 'ctf_ds'
      % convert CTF dataset into filenames
      [path, file, ext] = fileparts(cfg.dataset);
      cfg.headerfile = fullfile(cfg.dataset, [file '.res4']);
      cfg.datafile   = fullfile(cfg.dataset, [file '.meg4']);
    case 'brainvision_vhdr'
      [path, file, ext] = fileparts(cfg.dataset);
      cfg.headerfile = fullfile(path, [file '.vhdr']);
      if exist(fullfile(path, [file '.eeg']))
        cfg.datafile   = fullfile(path, [file '.eeg']);
      elseif exist(fullfile(path, [file '.seg']))
        cfg.datafile   = fullfile(path, [file '.seg']);
      end
    case 'brainvision_eeg'
      [path, file, ext] = fileparts(cfg.dataset);
      cfg.headerfile = fullfile(path, [file '.vhdr']);
      cfg.datafile   = fullfile(path, [file '.eeg']);
    case 'brainvision_seg'
      [path, file, ext] = fileparts(cfg.dataset);
      cfg.headerfile = fullfile(path, [file '.vhdr']);
      cfg.datafile   = fullfile(path, [file '.seg']);
    case {'physionet_dat', 'physionet_hea'}
      [path, file, ext] = fileparts(cfg.dataset);
      cfg.headerfile = fullfile(path, [file '.hea']);
      cfg.datafile   = fullfile(path, [file '.dat']);
    otherwise
      % convert dataset into filenames, assume that the header and data are the same
      cfg.datafile   = cfg.dataset;
      cfg.headerfile = cfg.dataset;
  end
  
elseif ~isempty(cfg.datafile) && isempty(cfg.headerfile)
  % assume that  the datafile also contains the header
  cfg.headerfile = cfg.datafile;

elseif isempty(cfg.datafile) && ~isempty(cfg.headerfile)
  % assume that  the headerfile also contains the data
  cfg.datafile = cfg.headerfile;
end

