function cfg = dataset2files(cfg);

% Helper function that converts a dataset into headerfile and datafile
% if neccessary. This is used in PREPROCESSING and DEFINETRIAL
%
% This function operates only on
%   cfg.dataset
%   cfg.datafile
%   cfg.headerfile
% and returns the updated cfg.

% Copyright (C) 2004, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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
  if strcmp(cfg.dataset, 'gui');
    d = uigetdir;
    if d==0
      [f, p] = uigetfile;
      if f==0
        error('You should select a dataset file or directory');
      else
        d = fullfile(p, f);
      end
    end
    cfg.dataset = d;
  end
  
    
  
  switch filetype(cfg.dataset)
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
    otherwise
      % convert dataset into filenames, assume that the header and data are the same
      cfg.datafile   = cfg.dataset;
      cfg.headerfile = cfg.dataset;
  end
elseif ~isempty(cfg.datafile) && isempty(cfg.headerfile);
  % assume that  the datafile also contains the header
  cfg.headerfile = cfg.datafile;
elseif isempty(cfg.datafile) && ~isempty(cfg.headerfile);
  % assume that  the headerfile also contains the data
  cfg.datafile = cfg.headerfile;
end

