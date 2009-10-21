function [cfg, artifact] = artifact_clip(cfg,data)

% ARTIFACT_CLIP scans the data segments of interest for channels that
% clip. A clipping artifact is detected by the signal being completely
% flat for some time.
%
% Use as
%   [cfg, artifact] = artifact_clip(cfg)
%   required configuration options:
%   cfg.dataset or both cfg.headerfile and cfg.datafile
% or
%   [cfg, artifact] = artifact_clip(cfg, data)
%   forbidden configuration options: 
%   cfg.dataset, cfg.headerfile and cfg.datafile
%
% In both cases the configuration should also contain:
%   cfg.artfctdef.clip.channel  = Nx1 cell-array with selection of channels, see CHANNELSELECTION for details
%   cfg.artfctdef.clip.pretim   = 0.000;  pre-artifact rejection-interval in seconds
%   cfg.artfctdef.clip.psttim   = 0.000;  post-artifact rejection-interval in seconds
%   cfg.artfctdef.clip.thresh   = 0.010;  minimum duration in seconds of a datasegment with consecutive identical samples to be considered as 'clipped'
%   cfg.continuous              = 'yes' or 'no' whether the file contains continuous data
%   
% See also REJECTARTIFACT

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: artifact_clip.m,v $
% Revision 1.22  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.21  2009/01/16 18:19:41  sashae
% moved some lines of code, no functional change
%
% Revision 1.20  2009/01/14 11:47:07  sashae
% changed handling of cfg.datatype
% added call to checkconfig at start and end of function
%
% Revision 1.19  2008/12/02 16:29:42  estmee
% Set default cfg.continuous/ checkconfig cfg.datatype = forbidden
%
% Revision 1.18  2008/11/25 13:12:48  estmee
% Documentation update
%
% Revision 1.17  2008/11/18 16:11:42  estmee
% Added cfg.continuous option to documentation.
%
% Revision 1.16  2008/10/13 14:37:46  estmee
% Changed the checkboundary again when calling read_data.
%
% Revision 1.15  2008/10/13 13:57:23  estmee
% Added checkboundary as input argument read_data.
%
% Revision 1.14  2008/10/13 11:39:34  sashae
% added call to checkconfig (as discussed with estmee)
%
% Revision 1.13  2008/10/07 16:16:31  estmee
% Added data as second input argument to artifact_clip.
%
% Revision 1.12  2008/10/07 08:58:51  roboos
% committed the changes that Esther made recently, related to the support of data as input argument to the artifact detection functions. I hope that this does not break the functions too seriously.
%
% Revision 1.11  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.10  2008/05/13 15:37:24  roboos
% switched to using read_data/header instead of the read_fcdc_data/header wrapper functions
%
% Revision 1.9  2006/11/29 09:06:36  roboos
% renamed all cfg options with "sgn" into "channel", added backward compatibility when required
% updated documentation, mainly in the artifact detection routines
%
% Revision 1.8  2006/06/14 12:43:48  roboos
% removed the documentation for cfg.lnfilttype, since that option is not supported by preproc
%
% Revision 1.7  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.6  2006/04/12 08:38:01  ingnie
% updated documentation
%
% Revision 1.5  2006/01/31 13:49:29  jansch
% included dataset2files to ensure the presence of cfg.headerfile or cfg.datafile
% whenever needed
%
% Revision 1.4  2005/12/20 08:36:47  roboos
% add the artifact Nx2 matrix to the output configuration
% changed some indentation and white space, renamed a few variables
%
% Revision 1.3  2005/10/06 14:28:20  roboos
% added check for case when no artifacts are found
%
% Revision 1.2  2005/09/05 10:37:03  roboos
% added pretim and psttim option for extending the detected artifacts
% added copyright and log
%

fieldtripdefs

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'trackconfig', 'on');
cfg = checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
cfg = checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});

% set default rejection parameters for clip artifacts if necessary.
if ~isfield(cfg,'artfctdef'),               cfg.artfctdef               = [];              end;
if ~isfield(cfg.artfctdef,'clip'),          cfg.artfctdef.clip          = [];              end;
if ~isfield(cfg.artfctdef.clip,'channel'),  cfg.artfctdef.clip.channel  = 'all';           end;
if ~isfield(cfg.artfctdef.clip,'thresh'),   cfg.artfctdef.clip.thresh   = 0.010;           end;
if ~isfield(cfg.artfctdef.clip,'pretim'),   cfg.artfctdef.clip.pretim   = 0.000;           end;
if ~isfield(cfg.artfctdef.clip,'psttim'),   cfg.artfctdef.clip.psttim   = 0.000;           end;

% for backward compatibility
if isfield(cfg.artfctdef.clip,'sgn')
  cfg.artfctdef.clip.channel = cfg.artfctdef.clip.sgn;
  cfg.artfctdef.clip         = rmfield(cfg.artfctdef.clip, 'sgn');
end

% start with an empty artifact list
artifact = [];

% read the header
if nargin == 1
  isfetch = 0;
  cfg = checkconfig(cfg, 'dataset2files', {'yes'});
  cfg = checkconfig(cfg, 'required', {'headerfile', 'datafile'});
  hdr = read_header(cfg.headerfile);
elseif nargin == 2
  isfetch = 1;
  cfg = checkconfig(cfg, 'forbidden', {'dataset', 'headerfile', 'datafile'});
  hdr = fetch_header(data);
end

% set default cfg.continuous
if ~isfield(cfg, 'continuous')
    if hdr.nTrials==1
      cfg.continuous = 'yes';
    else
      cfg.continuous = 'no';
    end
end

% find the channel labels present in the data and their indices
label = channelselection(cfg.artfctdef.clip.channel, hdr.label);
sgnindx = match_str(hdr.label, label);

% make a local copy for convenience
artfctdef = cfg.artfctdef.clip;

ntrl = size(cfg.trl,1);
nsgn = length(sgnindx);
for trlop=1:ntrl
  fprintf('searching for clipping artifacts in trial %d\n', trlop);
  % read the data of this trial
  if isfetch
    dat = fetch_data(data,        'header', hdr, 'begsample', cfg.trl(trlop,1), 'endsample', cfg.trl(trlop,2), 'chanindx', sgnindx);
  else
    dat = read_data(cfg.datafile, 'header', hdr, 'begsample', cfg.trl(trlop,1), 'endsample', cfg.trl(trlop,2), 'chanindx', sgnindx, 'checkboundary', strcmp(cfg.continuous, 'no'));
  end
  % apply filtering etc to the data
  datflt = preproc(dat, label, hdr.Fs, artfctdef, cfg.trl(trlop,3));
  % detect all samples that have the same value as the previous sample
  identical = (datflt(:,1:(end-1)) == datflt(:,2:end));
  % ensure that the number of samples does not change
  identical = [identical zeros(nsgn,1)];

  % determine the number of consecutively identical samples
  clip = zeros(size(dat));
  for sgnlop=1:length(sgnindx)
    up = find(diff([0 identical(sgnlop,:)], 1, 2)== 1);
    dw = find(diff([identical(sgnlop,:) 0], 1, 2)==-1);
    for k=1:length(up)
      clip(sgnlop,up(k):dw(k)) = dw(k)-up(k);
    end
  end
  % collapse over cannels
  clip = max(clip,[],1);

  % detect whether there are intervals in which the number of consecutive
  % identical samples is larger than the threshold
  thresh = (clip>=artfctdef.thresh*hdr.Fs);

  % remember the thresholded parts as artifacts
  artup = find(diff([0 thresh])== 1) + cfg.trl(trlop,1) - 1;
  artdw = find(diff([thresh 0])==-1) + cfg.trl(trlop,1) - 1;
  for k=1:length(artup)
    artifact(end+1,:) = [artup(k) artdw(k)];
  end
end

if ~isempty(artifact)
  % add the pretim and psttim to the detected artifacts
  artifact(:,1) = artifact(:,1) - artfctdef.pretim * hdr.Fs;
  artifact(:,2) = artifact(:,2) + artfctdef.psttim * hdr.Fs;
end

% remember the details that were used here
cfg.artfctdef.clip          = artfctdef;
cfg.artfctdef.clip.label    = label;
cfg.artfctdef.clip.trl      = cfg.trl;
cfg.artfctdef.clip.artifact = artifact;

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
cfg.version.id = '$Id: artifact_clip.m,v 1.22 2009/01/20 13:01:31 sashae Exp $';
