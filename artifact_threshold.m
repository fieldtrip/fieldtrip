function [cfg, artifact] = artifact_threshold(cfg,data)

% ARTIFACT_THRESHOLD scans for trials in which the range, i.e. the minimum,
% the maximum or the range (min-max difference) of the signal in any
% channel exceeds a specified threshold.
%
% use as:
%   [cfg, artifact] = artifact_threshold(cfg)
%   required configuration options: 
%   cfg.dataset or both cfg.headerfile and cfg.datafile
% or
%   [cfg, artifact] = artifact_threshold(cfg, data)
%   forbidden configuration options: 
%   cfg.dataset, cfg.headerfile and cfg.datafile
%
% In both cases the configuration should also contain:
%   cfg.continuous                    = 'yes' or 'no' whether the file contains continuous data
%
% The following configuration options can be specified
%   cfg.artfctdef.threshold.channel	  = cell-array with channel labels
%   cfg.artfctdef.threshold.bpfilter  = 'no' or 'yes'
%   cfg.artfctdef.threshold.bpfreq    = [0.3 30]
%   cfg.artfctdef.threshold.bpfiltord = 4
%
% The detection of artifacts is done according to the following settings,
% you should specify at least one of these thresholds
%   cfg.artfctdef.threshold.range     = value in uV/T, default  inf
%   cfg.artfctdef.threshold.min       = value in uV/T, default -inf
%   cfg.artfctdef.threshold.max       = value in uV/T, default  inf
%
% This function does not support partial rejections, since the whole trial
% is used to rate the minimum and maximum values. Furthermore, this
% function does not support artifact- or filterpadding.
%
% See also REJECTARTIFACT

% Copyright (c) 2003, Robert Oostenveld, SMI, FCDC
%
% $Log: artifact_threshold.m,v $
% Revision 1.29  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.28  2009/01/16 18:19:41  sashae
% moved some lines of code, no functional change
%
% Revision 1.27  2009/01/14 11:47:07  sashae
% changed handling of cfg.datatype
% added call to checkconfig at start and end of function
%
% Revision 1.26  2008/12/02 16:33:00  estmee
% Set default cfg.continuous/ checkconfig cfg.datatype= forbidden
%
% Revision 1.25  2008/11/25 13:55:31  estmee
% Documentation update
%
% Revision 1.24  2008/11/18 16:15:36  estmee
% Added cfg.continuous
%
% Revision 1.23  2008/10/13 11:39:34  sashae
% added call to checkconfig (as discussed with estmee)
%
% Revision 1.22  2008/10/07 16:15:30  estmee
% Added data as second intput argument to artifact_threshold and changed trllop in trlop several times.
%
% Revision 1.21  2008/10/07 08:58:51  roboos
% committed the changes that Esther made recently, related to the support of data as input argument to the artifact detection functions. I hope that this does not break the functions too seriously.
%
% Revision 1.20  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.19  2008/05/13 15:37:24  roboos
% switched to using read_data/header instead of the read_fcdc_data/header wrapper functions
%
% Revision 1.18  2006/06/14 12:43:52  roboos
% removed the documentation for cfg.lnfilttype, since that option is not supported by preproc
%
% Revision 1.17  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.16  2006/04/19 09:41:08  ingnie
% updated documentation
%
% Revision 1.15  2006/02/07 08:20:24  roboos
% fixed silly bug that was introduced along with dataset2files
%
% Revision 1.14  2006/01/31 13:49:29  jansch
% included dataset2files to ensure the presence of cfg.headerfile or cfg.datafile
% whenever needed
%
% Revision 1.13  2005/12/20 08:36:47  roboos
% add the artifact Nx2 matrix to the output configuration
% changed some indentation and white space, renamed a few variables
%
% Revision 1.12  2005/06/29 12:42:00  roboos
% added version to the output configuration
%
% Revision 1.11  2005/05/17 17:50:36  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.10  2005/02/02 14:29:49  roboos
% cleaned up code and help, implemented use of general preproc, changed some cfg fields (backward compatible), implemented min and max besides range, changed default from 75uV to Inf
%
% Revision 1.9  2005/01/25 13:15:41  roboos
% added check for cfg.datatype=continous and extended the call to read_fcdc_data with the boundary check for non-continous data
%
% Revision 1.8  2004/09/01 17:59:28  roboos
% added copyright statements to all filed
% added cfg.version to all functions that give configuration in their output
% added cfg.previous to all functions with input data containing configuration details
%
% Revision 1.7  2004/01/22 12:00:54  roberto
% updated documentation and help
%
% Revision 1.6  2004/01/21 22:04:12  roberto
% changed some whitespace
%
% Revision 1.5  2003/12/09 09:42:05  roberto
% aesthetic changes, not tested
%
% Revision 1.4  2003/12/08 12:30:39  roberto
% numerous changes
%
% Revision 1.3  2003/11/28 21:02:04  roberto
% added revision logging
%

fieldtripdefs

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'trackconfig', 'on');
cfg = checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
cfg = checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});

if ~isfield(cfg.artfctdef, 'threshold'), cfg.artfctdef.threshold = []; end

% copy the specific configuration for this function out of the master cfg
artfctdef = cfg.artfctdef.threshold;

% rename some cfg fields for backward compatibility
if isfield(artfctdef, 'sgn') && ~isfield(artfctdef, 'channel')
  artfctdef.channel = artfctdef.sgn;
  artfctdef         = rmfield(artfctdef, 'sgn');
end
if isfield(artfctdef, 'cutoff') && ~isfield(artfctdef, 'range')
  artfctdef.range = artfctdef.cutoff;
  artfctdef       = rmfield(artfctdef, 'cutoff');
end

% set default preprocessing parameters if necessary
if ~isfield(artfctdef, 'channel'),   artfctdef.channel   = 'all';    end
if ~isfield(artfctdef, 'bpfilter'),  artfctdef.bpfilter  = 'yes';    end
if ~isfield(artfctdef, 'bpfreq'),    artfctdef.bpfreq    = [0.3 30]; end
if ~isfield(artfctdef, 'bpfiltord'), artfctdef.bpfiltord = 4;        end

% set the default artifact detection parameters
if ~isfield(artfctdef, 'range'),    artfctdef.range = inf;           end
if ~isfield(artfctdef, 'min'),      artfctdef.min =  -inf;           end
if ~isfield(artfctdef, 'max'),      artfctdef.max =   inf;           end

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

% get the remaining settings
numtrl      = size(cfg.trl,1);
channel     = channelselection(artfctdef.channel, hdr.label);
channelindx = match_str(hdr.label,channel);
artifact    = [];

for trlop = 1:numtrl
  if isfetch
    dat = fetch_data(data,        'header', hdr, 'begsample', cfg.trl(trlop,1), 'endsample', cfg.trl(trlop,2), 'chanindx', channelindx, 'checkboundary', strcmp(cfg.continuous, 'no'));
  else
    dat = read_data(cfg.datafile, 'header', hdr, 'begsample', cfg.trl(trlop,1), 'endsample', cfg.trl(trlop,2), 'chanindx', channelindx, 'checkboundary', strcmp(cfg.continuous, 'no'));
  end
  dat = preproc(dat, channel, hdr.Fs, artfctdef, cfg.trl(trlop,3));
  % compute the min, max and range over all channels and samples
  minval   = min(dat(:));
  maxval   = max(dat(:));
  rangeval = maxval-minval;
  % test the min, max and range against the specified thresholds
  if ~isempty(artfctdef.min) && minval<artfctdef.min
    fprintf('threshold artifact scanning: trial %d from %d exceeds min-threshold\n', trlop, numtrl);
    artifact(end+1,1:2) = cfg.trl(trlop,1:2);
  elseif ~isempty(artfctdef.max) && maxval>artfctdef.max
    fprintf('threshold artifact scanning: trial %d from %d exceeds max-threshold\n', trlop, numtrl);
    artifact(end+1,1:2) = cfg.trl(trlop,1:2);
  elseif ~isempty(artfctdef.range) && rangeval>artfctdef.range
    fprintf('threshold artifact scanning: trial %d from %d exceeds range-threshold\n', trlop, numtrl);
    artifact(end+1,1:2) = cfg.trl(trlop,1:2);
  else
    fprintf('threshold artifact scanning: trial %d from %d is ok\n', trlop, numtrl);
  end
end

% remember the details that were used here
cfg.artfctdef.threshold          = artfctdef;
cfg.artfctdef.threshold.trl      = cfg.trl;         % trialdefinition prior to rejection
cfg.artfctdef.threshold.channel  = channel;         % exact channels used for detection
cfg.artfctdef.threshold.artifact = artifact;        % detected artifacts

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
cfg.version.id = '$Id: artifact_threshold.m,v 1.29 2009/01/20 13:01:31 sashae Exp $';
