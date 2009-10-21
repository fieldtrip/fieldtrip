function [cfg, artifact] = artifact_jump(cfg,data)

% ARTIFACT_JUMP reads the data segments of interest from file and identifies 
% SQUID jump artifacts.
%
% Use as
%   [cfg, artifact] = artifact_jump(cfg)
%   required configuration options: 
%   cfg.dataset or both cfg.headerfile and cfg.datafile
% or
%   [cfg, artifact] = artifact_jump(cfg, data)
%   forbidden configuration options: 
%   cfg.dataset, cfg.headerfile and cfg.datafile
%
% In both cases the configuration should also contain:
%   cfg.trl        = structure that defines the data segments of interest. See DEFINETRIAL
%   cfg.continuous = 'yes' or 'no' whether the file contains continuous data
%
% The data is preprocessed (again) with the following configuration parameters,
% which are optimal for identifying jump artifacts:
%   cfg.artfctdef.jump.medianfilter  = 'yes'
%   cfg.artfctdef.jump.medianfiltord = 9
%   cfg.artfctdef.jump.absdiff       = 'yes'
%
% Artifacts are identified by means of thresholding the z-transformed value
% of the preprocessed data.
%   cfg.artfctdef.jump.channel       = Nx1 cell-array with selection of channels, see CHANNELSELECTION for details
%   cfg.artfctdef.jump.cutoff        = 20      z-value at which to threshold
%   cfg.artfctdef.jump.trlpadding    = automatically determined based on the filter padding (cfg.padding)
%   cfg.artfctdef.jump.artpadding    = automatically determined based on the filter padding (cfg.padding)
%
% The output argument "artifact" is a Nx2 matrix comparable to the
% "trl" matrix of DEFINETRIAL. The first column of which specifying the
% beginsamples of an artifact period, the second column contains the
% endsamples of the artifactperiods.
%
% See also ARTIFACT_ZVALUE, REJECTARTIFACT

% Undocumented local options:
% cfg.method

% Copyright (c) 2003-2006, Jan-Mathijs Schoffelen & Robert Oostenveld
%
% $Log: artifact_jump.m,v $
% Revision 1.27  2009/10/12 14:26:47  jansch
% added default for not taking cumulated z-value across channels for artifact
% identification
%
% Revision 1.26  2009/03/10 14:25:59  roboos
% fixed bug in copying of cfg.continuous to tmpcfg, also keek data and headerformat
%
% Revision 1.25  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.24  2009/01/14 11:47:07  sashae
% changed handling of cfg.datatype
% added call to checkconfig at start and end of function
%
% Revision 1.23  2008/12/02 16:35:32  estmee
% Checkconfig cfg.datatype= forbidden
%
% Revision 1.22  2008/11/25 13:16:24  estmee
% Documentation update
%
% Revision 1.21  2008/11/18 16:20:58  estmee
% Added cfg.continuous
%
% Revision 1.20  2008/10/13 13:03:11  sashae
% added call to checkconfig (as discussed with estmee)
%
% Revision 1.19  2008/10/10 15:01:13  estmee
% Repaired determining artifacts with only cfg as input argument.
%
% Revision 1.18  2008/10/07 16:13:44  estmee
% Added data as second intput argument to artifact_jump itself and the way it calls artifact_zvalue.
%
% Revision 1.17  2008/10/07 08:58:51  roboos
% committed the changes that Esther made recently, related to the support of data as input argument to the artifact detection functions. I hope that this does not break the functions too seriously.
%
% Revision 1.16  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.15  2006/11/29 09:06:36  roboos
% renamed all cfg options with "sgn" into "channel", added backward compatibility when required
% updated documentation, mainly in the artifact detection routines
%
% Revision 1.14  2006/04/25 17:06:28  ingnie
% updated documentation
%
% Revision 1.13  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.12  2006/02/27 17:01:07  roboos
% added tmpcfg.dataset, added try-end around copying of headerfile, datafile and dataset
%
% Revision 1.11  2006/01/17 14:05:27  roboos
% do preproc absdiff instead of rectify in combination with derivative, absdiff ensures the right order
%
% Revision 1.10  2006/01/12 13:51:38  roboos
% completely new implementation, all based upon the same artifact_zvalue code
% all preprocessing is now done consistently and the various paddings have been better defined
% the functions do not have any explicit support any more for non-continuous data
% the old artifact_xxx functions from JM have been renamed to xxx_old
%

fieldtripdefs

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'trackconfig', 'on');
cfg = checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
cfg = checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});

% set default rejection parameters
if ~isfield(cfg,'artfctdef'),                      cfg.artfctdef                 = [];              end
if ~isfield(cfg.artfctdef,'jump'),                 cfg.artfctdef.jump            = [];              end
if ~isfield(cfg.artfctdef.jump,'method'),          cfg.artfctdef.jump.method     = 'zvalue';        end

% for backward compatibility
if isfield(cfg.artfctdef.jump,'sgn')
  cfg.artfctdef.jump.channel = cfg.artfctdef.jump.sgn;
  cfg.artfctdef.jump         = rmfield(cfg.artfctdef.jump, 'sgn');
end

if isfield(cfg.artfctdef.jump, 'artifact')
  fprintf('jump artifact detection has already been done, retaining artifacts\n');
  artifact = cfg.artfctdef.jump.artifact;
  return
end

if strcmp(cfg.artfctdef.jump.method, 'zvalue')
  % the following fields should be supported for backward compatibility
  dum = 0;
  if isfield(cfg.artfctdef.jump,'pretim'),
    dum = max(dum, cfg.artfctdef.jump.pretim);
    cfg.artfctdef.jump = rmfield(cfg.artfctdef.jump,'pretim');
  end
  if isfield(cfg.artfctdef.jump,'psttim'),
    dum = max(dum, cfg.artfctdef.jump.psttim);
    cfg.artfctdef.jump = rmfield(cfg.artfctdef.jump,'psttim');
  end
  if dum
    cfg.artfctdef.jump.artpadding = max(dum);
  end
  if isfield(cfg.artfctdef.jump,'padding'),
    cfg.artfctdef.jump.trlpadding   = cfg.artfctdef.jump.padding;
    cfg.artfctdef.jump = rmfield(cfg.artfctdef.jump,'padding');
  end
  % settings for preprocessing
  if ~isfield(cfg.artfctdef.jump,'medianfilter'),  cfg.artfctdef.jump.medianfilter  = 'yes';        end
  if ~isfield(cfg.artfctdef.jump,'medianfiltord'), cfg.artfctdef.jump.medianfiltord = 9;            end
  if ~isfield(cfg.artfctdef.jump,'absdiff'),       cfg.artfctdef.jump.absdiff       = 'yes';        end  % compute abs(diff(data)), whereas the order of rectify=yes in combination with derivative=yes would be diff(abs(data)) due to the ordering in preproc
  % settings for the zvalue subfunction
  if ~isfield(cfg.artfctdef.jump,'cutoff'),        cfg.artfctdef.jump.cutoff     = 20;              end
  if ~isfield(cfg.artfctdef.jump,'channel'),       cfg.artfctdef.jump.channel    = 'MEG';           end
  if ~isfield(cfg.artfctdef.jump,'cumulative'),    cfg.artfctdef.jump.cumulative = 'no';            end
  if isfield(cfg, 'padding') && cfg.padding~=0
    if ~isfield(cfg.artfctdef.jump,'trlpadding'), cfg.artfctdef.jump.trlpadding = 0.5*cfg.padding; end
    if ~isfield(cfg.artfctdef.jump,'artpadding'), cfg.artfctdef.jump.artpadding = 0.5*cfg.padding; end
    if ~isfield(cfg.artfctdef.jump,'fltpadding'), cfg.artfctdef.jump.fltpadding = 0;               end
  else
    if ~isfield(cfg.artfctdef.jump,'trlpadding'), cfg.artfctdef.jump.trlpadding = 0;               end
    if ~isfield(cfg.artfctdef.jump,'artpadding'), cfg.artfctdef.jump.artpadding = 0;               end
    if ~isfield(cfg.artfctdef.jump,'fltpadding'), cfg.artfctdef.jump.fltpadding = 0;               end
  end
  % construct a temporary configuration that can be passed onto artifact_zvalue
  tmpcfg                  = [];
  tmpcfg.trl              = cfg.trl;
  tmpcfg.artfctdef.zvalue = cfg.artfctdef.jump;
  if isfield(cfg, 'continuous'),   tmpcfg.continuous       = cfg.continuous;    end
  if isfield(cfg, 'dataformat'),   tmpcfg.dataformat       = cfg.dataformat;    end
  if isfield(cfg, 'headerformat'), tmpcfg.headerformat     = cfg.headerformat;  end
  % call the zvalue artifact detection function
  if nargin ==1
    cfg = checkconfig(cfg, 'dataset2files', {'yes'});
    cfg = checkconfig(cfg, 'required', {'headerfile', 'datafile'});  
    tmpcfg.datafile    = cfg.datafile;
    tmpcfg.headerfile  = cfg.headerfile;
    [tmpcfg, artifact] = artifact_zvalue(tmpcfg);
  elseif nargin ==2
    cfg = checkconfig(cfg, 'forbidden', {'dataset', 'headerfile', 'datafile'});
    [tmpcfg, artifact] = artifact_zvalue(tmpcfg, data);
  end
  cfg.artfctdef.jump = tmpcfg.artfctdef.zvalue;
else
  error(sprintf('jump artifact detection only works with cfg.method=''zvalue'''));
end

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 
