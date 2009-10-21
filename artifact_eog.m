function [cfg, artifact] = artifact_eog(cfg,data)

% ARTIFACT_EOG reads the data segments of interest from file and
% identifies EOG artifacts.
%
% Use as
%   [cfg, artifact] = artifact_eog(cfg)
%   required configuration options: 
%   cfg.dataset or both cfg.headerfile and cfg.datafile
% or
%   [cfg, artifact] = artifact_eog(cfg, data)
%   forbidden configuration options: 
%   cfg.dataset, cfg.headerfile and cfg.datafile
%
% In both cases the configuration should also contain:
%   cfg.trl        = structure that defines the data segments of interest. See DEFINETRIAL
%   cfg.continuous = 'yes' or 'no' whether the file contains continuous data
%
% The data is preprocessed (again) with the following configuration parameters,
% which are optimal for identifying EOG artifacts:
%   cfg.artfctdef.eog.bpfilter   = 'yes'
%   cfg.artfctdef.eog.bpfilttype = 'but'
%   cfg.artfctdef.eog.bpfreq     = [1 15]
%   cfg.artfctdef.eog.bpfiltord  = 4
%   cfg.artfctdef.eog.hilbert    = 'yes'
%
% Artifacts are identified by means of thresholding the z-transformed value
% of the preprocessed data.
%   cfg.artfctdef.eog.channel      = Nx1 cell-array with selection of channels, see CHANNELSELECTION for details
%   cfg.artfctdef.eog.cutoff       = 4       z-value at which to threshold
%   cfg.artfctdef.eog.trlpadding   = 0.5
%   cfg.artfctdef.eog.fltpadding   = 0.1
%   cfg.artfctdef.eog.artpadding   = 0.1
%
% The output argument "artifact" is a Nx2 matrix comparable to the
% "trl" matrix of DEFINETRIAL. The first column of which specifying the
% beginsamples of an artifact period, the second column contains the
% endsamples of the artifactperiods.
%
% See also ARTIFACT_ZVALUE, REJECTARTIFACT

% Undocumented local options
% cfg.method

% Copyright (c) 2003-2006, Jan-Mathijs Schoffelen & Robert Oostenveld
%
% $Log: artifact_eog.m,v $
% Revision 1.36  2009/03/10 14:25:59  roboos
% fixed bug in copying of cfg.continuous to tmpcfg, also keek data and headerformat
%
% Revision 1.35  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.34  2009/01/14 11:47:07  sashae
% changed handling of cfg.datatype
% added call to checkconfig at start and end of function
%
% Revision 1.33  2008/12/02 16:35:01  estmee
% Checkconfig cfg.datatype = forbidden
%
% Revision 1.32  2008/11/25 13:14:28  estmee
% Documentation update
%
% Revision 1.31  2008/11/18 16:20:28  estmee
% Added cfg.continuous
%
% Revision 1.30  2008/10/13 13:03:11  sashae
% added call to checkconfig (as discussed with estmee)
%
% Revision 1.29  2008/10/10 15:00:10  estmee
% Repaired determining artifacts with cfg as only input argument.
%
% Revision 1.28  2008/10/07 16:11:08  estmee
% Added data as second input argument to artifact_eog itself and to the way it calls artifact_zvalue.
%
% Revision 1.27  2008/10/07 08:58:51  roboos
% committed the changes that Esther made recently, related to the support of data as input argument to the artifact detection functions. I hope that this does not break the functions too seriously.
%
% Revision 1.26  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.25  2006/11/29 09:06:36  roboos
% renamed all cfg options with "sgn" into "channel", added backward compatibility when required
% updated documentation, mainly in the artifact detection routines
%
% Revision 1.24  2006/04/25 17:06:28  ingnie
% updated documentation
%
% Revision 1.23  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.22  2006/02/27 17:01:07  roboos
% added tmpcfg.dataset, added try-end around copying of headerfile, datafile and dataset
%
% Revision 1.21  2006/01/12 13:51:38  roboos
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
if ~isfield(cfg,'artfctdef'),                  cfg.artfctdef                 = [];       end
if ~isfield(cfg.artfctdef,'eog'),              cfg.artfctdef.eog             = [];       end
if ~isfield(cfg.artfctdef.eog,'method'),       cfg.artfctdef.eog.method      = 'zvalue'; end

% for backward compatibility
if isfield(cfg.artfctdef.eog,'sgn')
  cfg.artfctdef.eog.channel = cfg.artfctdef.eog.sgn;
  cfg.artfctdef.eog         = rmfield(cfg.artfctdef.eog, 'sgn');
end

if isfield(cfg.artfctdef.eog, 'artifact')
  fprintf('eog artifact detection has already been done, retaining artifacts\n');
  artifact = cfg.artfctdef.eog.artifact;
  return
end

if strcmp(cfg.artfctdef.eog.method, 'zvalue')
  % the following fields should be supported for backward compatibility
  if isfield(cfg.artfctdef.eog,'pssbnd'),
    cfg.artfctdef.eog.bpfreq   = cfg.artfctdef.eog.pssbnd;
    cfg.artfctdef.eog.bpfilter = 'yes';
    cfg.artfctdef.eog = rmfield(cfg.artfctdef.eog,'pssbnd');
  end;
  dum = 0;
  if isfield(cfg.artfctdef.eog,'pretim'),
    dum = max(dum, cfg.artfctdef.eog.pretim);
    cfg.artfctdef.eog = rmfield(cfg.artfctdef.eog,'pretim');
  end
  if isfield(cfg.artfctdef.eog,'psttim'),
    dum = max(dum, cfg.artfctdef.eog.psttim);
    cfg.artfctdef.eog = rmfield(cfg.artfctdef.eog,'psttim');
  end
  if dum
    cfg.artfctdef.eog.artpadding = max(dum);
  end
  if isfield(cfg.artfctdef.eog,'padding'),
    cfg.artfctdef.eog.trlpadding   = cfg.artfctdef.eog.padding;
    cfg.artfctdef.eog = rmfield(cfg.artfctdef.eog,'padding');
  end
  % settings for preprocessing
  if ~isfield(cfg.artfctdef.eog,'bpfilter'),   cfg.artfctdef.eog.bpfilter   = 'yes';     end
  if ~isfield(cfg.artfctdef.eog,'bpfilttype'), cfg.artfctdef.eog.bpfilttype = 'but';     end
  if ~isfield(cfg.artfctdef.eog,'bpfreq'),     cfg.artfctdef.eog.bpfreq     = [1 15];    end
  if ~isfield(cfg.artfctdef.eog,'bpfiltord'),  cfg.artfctdef.eog.bpfiltord  = 4;         end
  if ~isfield(cfg.artfctdef.eog,'hilbert'),    cfg.artfctdef.eog.hilbert    = 'yes';     end
  % settings for the zvalue subfunction
  if ~isfield(cfg.artfctdef.eog,'channel'),    cfg.artfctdef.eog.channel     = 'EOG';    end
  if ~isfield(cfg.artfctdef.eog,'trlpadding'), cfg.artfctdef.eog.trlpadding  = 0.5;      end
  if ~isfield(cfg.artfctdef.eog,'artpadding'), cfg.artfctdef.eog.artpadding  = 0.1;      end
  if ~isfield(cfg.artfctdef.eog,'fltpadding'), cfg.artfctdef.eog.fltpadding  = 0.1;      end
  if ~isfield(cfg.artfctdef.eog,'cutoff'),     cfg.artfctdef.eog.cutoff      = 4;        end
  % construct a temporary configuration that can be passed onto artifact_zvalue
  tmpcfg                  = [];
  tmpcfg.trl              = cfg.trl;
  tmpcfg.artfctdef.zvalue = cfg.artfctdef.eog;
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
  cfg.artfctdef.eog  = tmpcfg.artfctdef.zvalue;
else
  error(sprintf('EOG artifact detection only works with cfg.method=''zvalue'''));
end

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');
