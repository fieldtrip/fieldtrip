function [cfg, artifact] = artifact_file(cfg);

% ARTIFACT_FILE reads rejection marks from a file
%
% Use as
%   [cfg, artifact] = arifact_file(cfg)
%   required configuration options: 
%   cfg.dataset or cfg.headerfile
%
% See also REJECTARTIFACT

% Copyright (C) 2003-2006, Robert Oostenveld
%
% $Log: artifact_file.m,v $
% Revision 1.17  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.16  2009/01/16 17:21:20  sashae
% added config tracking
%
% Revision 1.15  2008/11/25 13:15:42  estmee
% Documentation update
%
% Revision 1.14  2008/10/13 10:40:47  sashae
% added call to checkconfig
%
% Revision 1.13  2008/10/07 08:58:51  roboos
% committed the changes that Esther made recently, related to the support of data as input argument to the artifact detection functions. I hope that this does not break the functions too seriously.
%
% Revision 1.12  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.11  2008/05/13 15:37:24  roboos
% switched to using read_data/header instead of the read_fcdc_data/header wrapper functions
%
% Revision 1.10  2006/05/03 08:12:51  ingnie
% updated documentation
%
% Revision 1.9  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.8  2006/01/31 13:49:29  jansch
% included dataset2files to ensure the presence of cfg.headerfile or cfg.datafile
% whenever needed
%
% Revision 1.7  2005/12/20 08:36:47  roboos
% add the artifact Nx2 matrix to the output configuration
% changed some indentation and white space, renamed a few variables
%
% Revision 1.6  2005/06/29 12:42:00  roboos
% added version to the output configuration
%
% Revision 1.5  2005/05/17 17:50:36  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.4  2004/09/01 17:59:28  roboos
% added copyright statements to all filed
% added cfg.version to all functions that give configuration in their output
% added cfg.previous to all functions with input data containing configuration details
%
% Revision 1.3  2004/08/26 16:01:34  roboos
% added support for brainvision_marker
%
% Revision 1.2  2003/12/11 09:35:27  roberto
% fixed bug in "if isfield()"
%
% Revision 1.1  2003/12/08 12:34:27  roberto
% initial version, according to pluggable artifact rejection
%

fieldtripdefs
cfg = checkconfig(cfg, 'trackconfig', 'on');

if isfield(cfg, 'rejectfile') && ~strcmp(cfg.rejectfile, 'no')
  cfg = checkconfig(cfg, 'dataset2files', {'yes'});
  cfg = checkconfig(cfg, 'required', {'headerfile'});
  hdr = read_header(cfg.headerfile);
  if filetype(cfg.rejectfile, 'eep_rej')
    artifact = read_eep_rej(cfg.rejectfile);
  elseif filetype(cfg.rejectfile, 'brainvision_marker')
    artifact = read_brainvision_marker(cfg.rejectfile);
  else
    error(sprintf('unrecognized filetype for rejection file ''%s''', cfg.rejectfile));
  end
  % convert rejection marks from miliseconds into seconds
  artifact = artifact/1000;
  % convert rejection marks into samples
  artifact = round(artifact * hdr.Fs);
  fprintf('%d rejection marks read\n', size(artifact,1));
else
  artifact = [];
  fprintf('no rejection marks read\n');
end

% remember the details that were used here
cfg.artfctdef.file.trl      = cfg.trl;
cfg.artfctdef.file.artifact = artifact;

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: artifact_file.m,v 1.17 2009/01/20 13:01:31 sashae Exp $';
