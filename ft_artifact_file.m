function [cfg, artifact] = ft_artifact_file(cfg);

% FT_ARTIFACT_FILE reads rejection marks from a file
%
% Use as
%   [cfg, artifact] = ft_arifact_file(cfg)
%   required configuration options: 
%   cfg.dataset or cfg.headerfile
%
% See also FT_REJECTARTIFACT

% Copyright (C) 2003-2006, Robert Oostenveld
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

fieldtripdefs
cfg = checkconfig(cfg, 'trackconfig', 'on');

if isfield(cfg, 'rejectfile') && ~strcmp(cfg.rejectfile, 'no')
  cfg = checkconfig(cfg, 'dataset2files', {'yes'});
  cfg = checkconfig(cfg, 'required', {'headerfile'});
  hdr = ft_read_header(cfg.headerfile);
  if ft_filetype(cfg.rejectfile, 'eep_rej')
    artifact = read_eep_rej(cfg.rejectfile);
  elseif ft_filetype(cfg.rejectfile, 'brainvision_marker')
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
cfg.version.id = '$Id$';
