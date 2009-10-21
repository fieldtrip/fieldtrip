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
% $Log: dataset2files.m,v $
% Revision 1.7  2008/01/10 21:00:51  roboos
% added dataset=gui, pop up graphical interface (first asking for directory, then for file)
%
% Revision 1.6  2005/09/07 10:00:46  roboos
% fixed bug in path/file for brainvision
% added check for presence of *.eeg/seg file in case dataset=brainvision_vhdr
%
% Revision 1.5  2005/05/17 17:50:49  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.4  2005/02/16 08:05:16  roboos
% changed default behaviour so that dataset is copied into headerfile/datafile, except for ctf of brainvision
%
% Revision 1.3  2005/02/02 14:32:23  roboos
% added eep_cnt to the list of files for which dataset=headerfile+datafile
%
% Revision 1.2  2004/12/20 15:03:58  roboos
% added support for translating dataset to data+headerfile for neuroscan cnt and eeg
%
% Revision 1.1  2004/11/15 09:21:27  roboos
% simple helper function to avoid repliaction of code preprocessing and definetrial
%

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

