function [cfg, artifact] = artifact_muscle(cfg,data)

% ARTIFACT_MUSCLE reads the data segments of interest from file and
% identifies muscle artifacts.
%
% Use as
%   [cfg, artifact] = artifact_muscle(cfg)
%   required configuration options: 
%   cfg.dataset or both cfg.headerfile and cfg.datafile
% or
%   [cfg, artifact] = artifact_muscle(cfg, data)
%   forbidden configuration options: 
%   cfg.dataset, cfg.headerfile and cfg.datafile
%
% In both cases the configuration should also contain:
%   cfg.trl        = structure that defines the data segments of interest. See DEFINETRIAL
%   cfg.continuous = 'yes' or 'no' whether the file contains continuous data
%
% The data is preprocessed (again) with the following configuration parameters,
% which are optimal for identifying muscle artifacts:
%   cfg.artfctdef.muscle.bpfilter    = 'yes'
%   cfg.artfctdef.muscle.bpfreq      = [110 140]
%   cfg.artfctdef.muscle.bpfiltord   = 10
%   cfg.artfctdef.muscle.bpfilttype  = 'but'
%   cfg.artfctdef.muscle.hilbert     = 'yes'
%   cfg.artfctdef.muscle.boxcar      = 0.2
%
% Artifacts are identified by means of thresholding the z-transformed value
% of the preprocessed data.
%   cfg.artfctdef.muscle.channel     = Nx1 cell-array with selection of channels, see CHANNELSELECTION for details
%   cfg.artfctdef.muscle.cutoff      = 4       z-value at which to threshold
%   cfg.artfctdef.muscle.trlpadding  = 0.1
%   cfg.artfctdef.muscle.fltpadding  = 0.1
%   cfg.artfctdef.muscle.artpadding  = 0.1
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
% $Log: artifact_muscle.m,v $
% Revision 1.32  2009/03/10 14:25:59  roboos
% fixed bug in copying of cfg.continuous to tmpcfg, also keek data and headerformat
%
% Revision 1.31  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.30  2009/01/14 11:47:07  sashae
% changed handling of cfg.datatype
% added call to checkconfig at start and end of function
%
% Revision 1.29  2008/12/02 16:36:10  estmee
% Checkconfig cfg.datatype = forbidden
%
% Revision 1.28  2008/11/25 13:17:56  estmee
% Documentation update
%
% Revision 1.27  2008/11/18 16:21:28  estmee
% Added cfg.continuous
%
% Revision 1.26  2008/10/13 13:03:11  sashae
% added call to checkconfig (as discussed with estmee)
%
% Revision 1.25  2008/10/10 15:02:07  estmee
% Repaired determining artifacts with cfg as only input argument.
%
% Revision 1.24  2008/10/07 16:12:31  estmee
% Added data as second input argument to artifact_muscle itself and to the way is calls artifact_zvalue.
%
% Revision 1.23  2008/10/07 08:58:51  roboos
% committed the changes that Esther made recently, related to the support of data as input argument to the artifact detection functions. I hope that this does not break the functions too seriously.
%
% Revision 1.22  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.21  2006/12/04 09:30:02  roboos
% changed a comment
%
% Revision 1.20  2006/11/29 09:06:36  roboos
% renamed all cfg options with "sgn" into "channel", added backward compatibility when required
% updated documentation, mainly in the artifact detection routines
%
% Revision 1.19  2006/04/25 17:06:28  ingnie
% updated documentation
%
% Revision 1.18  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.17  2006/02/27 17:01:07  roboos
% added tmpcfg.dataset, added try-end around copying of headerfile, datafile and dataset
%
% Revision 1.16  2006/01/12 13:51:38  roboos
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
if ~isfield(cfg,'artfctdef'),                     cfg.artfctdef                    = [];        end
if ~isfield(cfg.artfctdef,'muscle'),              cfg.artfctdef.muscle             = [];        end
if ~isfield(cfg.artfctdef.muscle,'method'),       cfg.artfctdef.muscle.method      = 'zvalue';  end

% for backward compatibility
if isfield(cfg.artfctdef.muscle,'sgn')
  cfg.artfctdef.muscle.channel = cfg.artfctdef.muscle.sgn;
  cfg.artfctdef.muscle         = rmfield(cfg.artfctdef.muscle, 'sgn');
end

if isfield(cfg.artfctdef.muscle, 'artifact')
  fprintf('muscle artifact detection has already been done, retaining artifacts\n');
  artifact = cfg.artfctdef.muscle.artifact;
  return
end

if strcmp(cfg.artfctdef.muscle.method, 'zvalue')
  % the following settings should be supported for backward compatibility
  if isfield(cfg.artfctdef.muscle,'pssbnd'),
    cfg.artfctdef.muscle.bpfreq   = cfg.artfctdef.muscle.pssbnd;
    cfg.artfctdef.muscle.bpfilter = 'yes';
    cfg.artfctdef.muscle = rmfield(cfg.artfctdef.muscle,'pssbnd');
  end;
  dum = 0;
  if isfield(cfg.artfctdef.muscle,'pretim'),
    dum = max(dum, cfg.artfctdef.muscle.pretim);
    cfg.artfctdef.muscle = rmfield(cfg.artfctdef.muscle,'pretim');
  end
  if isfield(cfg.artfctdef.muscle,'psttim'),
    dum = max(dum, cfg.artfctdef.muscle.psttim);
    cfg.artfctdef.muscle = rmfield(cfg.artfctdef.muscle,'psttim');
  end
  if dum
    cfg.artfctdef.muscle.artpadding = max(dum);
  end
  if isfield(cfg.artfctdef.muscle,'padding'),
    cfg.artfctdef.muscle.trlpadding   = cfg.artfctdef.muscle.padding;
    cfg.artfctdef.muscle = rmfield(cfg.artfctdef.muscle,'padding');
  end

  % settings for preprocessing
  if ~isfield(cfg.artfctdef.muscle,'bpfilter'),   cfg.artfctdef.muscle.bpfilter    = 'yes';     end
  if ~isfield(cfg.artfctdef.muscle,'bpfreq'),     cfg.artfctdef.muscle.bpfreq      = [110 140]; end
  if ~isfield(cfg.artfctdef.muscle,'bpfiltord'),  cfg.artfctdef.muscle.bpfiltord   = 10;        end
  if ~isfield(cfg.artfctdef.muscle,'bpfilttype'), cfg.artfctdef.muscle.bpfilttype  = 'but';     end
  if ~isfield(cfg.artfctdef.muscle,'hilbert'),    cfg.artfctdef.muscle.hilbert     = 'yes';     end
  if ~isfield(cfg.artfctdef.muscle,'boxcar'),     cfg.artfctdef.muscle.boxcar      = 0.2;       end
  % settings for the zvalue subfunction
  if ~isfield(cfg.artfctdef.muscle,'channel'),    cfg.artfctdef.muscle.channel     = 'MEG';     end
  if ~isfield(cfg.artfctdef.muscle,'trlpadding'), cfg.artfctdef.muscle.trlpadding  = 0.1;       end
  if ~isfield(cfg.artfctdef.muscle,'fltpadding'), cfg.artfctdef.muscle.fltpadding  = 0.1;       end
  if ~isfield(cfg.artfctdef.muscle,'artpadding'), cfg.artfctdef.muscle.artpadding  = 0.1;       end
  if ~isfield(cfg.artfctdef.muscle,'cutoff'),     cfg.artfctdef.muscle.cutoff      = 4;         end
  % construct a temporary configuration that can be passed onto artifact_zvalue
  tmpcfg                  = [];
  tmpcfg.trl              = cfg.trl;
  tmpcfg.artfctdef.zvalue = cfg.artfctdef.muscle;
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
  cfg.artfctdef.muscle = tmpcfg.artfctdef.zvalue;
else
  error(sprintf('muscle artifact detection only works with cfg.method=''zvalue'''));
end

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 
