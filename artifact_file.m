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
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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
