function [cfg] = rejectartifact(cfg,data)

% REJECTARTIFACT removes data segments containing artifacts. It returns a
% configuration structure with a modified trial definition which can be
% used for preprocessing of only the clean data.
%
% You should start by detecting the artifacts in the data using the
% function ARTIFACT_xxx where xxx is the type of artifact. Subsequently
% REJECTARTIFACT looks at the detected artifacts and removes them from
% the trial definition or from the data.
%
% Use as
%   cfg = rejectartifact(cfg)
% with the cfg as obtained from DEFINETRIAL, or as
%   data = rejectartifact(cfg, data)
% with the data as obtained from PREPROCESSING
%
% The following configuration options are supported:
%   cfg.artfctdef.reject          = 'none', 'partial' or 'complete' (default = 'complete')
%   cfg.artfctdef.minaccepttim    = length in seconds (default = 0.1)
%   cfg.artfctdef.feedback        = 'yes' or 'no' (default = 'no')
%   cfg.artfctdef.eog.artifact    = Nx2 matrix with artifact segments, this is added to the cfg by using ARTIFACT_EOG
%   cfg.artfctdef.jump.artifact   = Nx2 matrix with artifact segments, this is added to the cfg by using ARTIFACT_JUMP
%   cfg.artfctdef.muscle.artifact = Nx2 matrix with artifact segments, this is added to the cfg by using ARTIFACT_MUSCLE
%   cfg.artfctdef.zvalue.artifact = Nx2 matrix with artifact segments, this is added to the cfg by using ARTIFACT_ZVALUE
%   cfg.artfctdef.xxx.artifact    = Nx2 matrix with artifact segments, this should be added by your own artifact detection function
%
% A trial that contains an artifact can be rejected completely or
% partially. In case of partial rejection, a minimum length of the
% resulting sub-trials can be specified.
%
% Output:
%   If cfg is used as the only input parameter, a cfg with a new trl is the output.
%   If cfg and data are both input parameters, a new raw data structure with only the clean data segments is the output.
%
% See also ARTIFACT_EOG, ARTIFACT_MUSCLE, ARTIFACT_JUMP, ARTIFACT_MANUAL, ARTIFACT_THRESHOLD, ARTIFACT_CLIP, ARTIFACT_ECG

% Undocumented local options:
% cfg.headerfile
% cfg.reject
% cfg.trl
% cfg.trlold
% cfg.version
% These old configuration options are still supported
% cfg.rejectmuscle      = 'no' or 'yes'
% cfg.rejecteog         = 'no' or 'yes'
% cfg.rejectjump        = 'no' or 'yes'
% cfg.rejectfile        = string with filename
% cfg.artfctdef.writerej = filename of rejection file
% cfg.artfctdef.type    = cell-array with strings, e.g. {'eog', 'muscle' 'jump'}

% Copyright (C) 2003-2007, Robert Oostenveld
%
% $Log: rejectartifact.m,v $
% Revision 1.48  2009/06/23 18:33:17  roboos
% use the trl from the input data and not from the config in case of nargin>1
%
% Revision 1.47  2009/03/23 21:20:10  roboos
% removed extra ;
%
% Revision 1.46  2009/01/27 17:11:42  roboos
% updated the help
%
% Revision 1.45  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.44  2009/01/14 11:29:55  sashae
% temporarily disabled previous revision
%
% Revision 1.43  2009/01/13 10:14:58  sashae
% changed handling of the output cfg: now the cfg also has cfg.previous fields,
% similar to data.cfg.previous. this way the output of definetrial and the
% artifact functions is kept separately from subsequent preprocessing steps
%
% Revision 1.42  2008/10/13 13:54:38  estmee
% Documentation is updated and added fetch_header (used when there are 2 input arguments).
%
% Revision 1.41  2008/10/13 13:05:07  sashae
% replaced call to dataset2files with checkconfig
%
% Revision 1.40  2008/10/07 16:06:27  estmee
% added error when "no trials left"
%
% Revision 1.39  2008/10/07 08:58:51  roboos
% committed the changes that Esther made recently, related to the support of data as input argument to the artifact detection functions. I hope that this does not break the functions too seriously.
%
% Revision 1.38  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.37  2008/05/13 15:37:24  roboos
% switched to using read_data/header instead of the read_fcdc_data/header wrapper functions
%
% Revision 1.36  2007/03/06 16:15:43  roboos
% keep the values of the 4th and subsequent column in the trl matrix, also when doing partial artifact rejection
%
% Revision 1.35  2007/01/04 17:08:36  roboos
% some reordering of code, no functional change
% added try-catch around the reading of the events
%
% Revision 1.34  2006/09/18 18:37:15  roboos
% updated documentation
%
% Revision 1.33  2006/06/26 15:17:33  roboos
% small change for supporting the rejectxxx backward compatible stuff, concatenate explicitely with cat function along 1st dimension, removed warnings
%
% Revision 1.32  2006/06/20 12:56:57  roboos
% updated documentation, removed unused variable in code
%
% Revision 1.31  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.30  2006/02/24 16:40:53  roboos
% added one fprintf information
%
% Revision 1.29  2006/01/17 12:28:58  roboos
% add version information to the cfg.artfctdef substructure
%
% Revision 1.28  2006/01/12 17:18:18  roboos
% fixed bug in testing which sub-structures contain artifacts
%
% Revision 1.27  2006/01/12 14:22:03  roboos
% fixed bug in negation of boolean rejectall vector
%
% Revision 1.26  2006/01/12 14:12:53  roboos
% convert rejectall vector containing artifact type into logical before proceeding with old code
%
% Revision 1.25  2006/01/12 13:48:14  roboos
% get the detected artifacts from the configuration structure instead of from the artifact_xxx function
% implemented graphical feedback using two plots (using cfg.artfctdef.feedback)
%

fieldtripdefs

if 0
  % this code snippet ensures that these functions are included in the
  % documentation as dependencies
  try, dum = artifact_ecg;       end
  try, dum = artifact_eog;       end
  try, dum = artifact_muscle;    end
  try, dum = artifact_jump;      end
  try, dum = artifact_clip;      end
  try, dum = artifact_manual;    end
  try, dum = artifact_threshold; end
end

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'trackconfig', 'on');
cfg = checkconfig(cfg, 'dataset2files', {'yes'});

% set the defaults
if ~isfield(cfg, 'artfctdef'),              cfg.artfctdef        = [];         end
if ~isfield(cfg.artfctdef,'type'),          cfg.artfctdef.type   = {};         end
if ~isfield(cfg.artfctdef,'reject'),        cfg.artfctdef.reject = 'complete'; end
if ~isfield(cfg.artfctdef,'minaccepttim'),  cfg.artfctdef.minaccepttim = 0.1;  end
if ~isfield(cfg.artfctdef,'feedback'),      cfg.artfctdef.feedback = 'no';     end

% convert from old-style to new-style configuration
if isfield(cfg,'reject')
  warning('converting from old-style artifact configuration to new-style');
  cfg.artfctdef.reject = cfg.reject;
  cfg = rmfield(cfg, 'reject');
end

% convert from old-style to new-style configuration
if isfield(cfg.artfctdef,'common')
  warning('converting from old-style artifact configuration to new-style');
  if isfield(cfg.artfctdef.common,'minaccepttim')
    cfg.artfctdef.minaccepttim = cfg.artfctdef.common.minaccepttim;
    cfg.artfctdef = rmfield(cfg.artfctdef, 'common');
  end
end

% ensure that it is a cell array
if ischar(cfg.artfctdef.type)
  cfg.artfctdef.type = {cfg.artfctdef.type};
end

% support the rejectXXX cfg settings for backward compatibility
if isfield(cfg, 'rejectmuscle')
  dum = strmatch('muscle', cfg.artfctdef.type, 'exact');
  if strcmp(cfg.rejectmuscle,'yes') && isempty(dum)
    % this overrules the other setting, add it to the type-list
    cfg.artfctdef.type = cat(1, {'muscle'}, cfg.artfctdef.type(:));
  elseif strcmp(cfg.rejectmuscle,'no') && ~isempty(dum)
    % this overrules the other setting, remove it from the type-list
    cfg.artfctdef.type(dum) = [];
  end
  cfg = rmfield(cfg, 'rejectmuscle');
end

% support the rejectXXX cfg settings for backward compatibility
if isfield(cfg, 'rejecteog')
  dum = strmatch('eog', cfg.artfctdef.type, 'exact');
  if strcmp(cfg.rejecteog,'yes') && isempty(dum)
    % this overrules the other setting, add it to the type-list
    cfg.artfctdef.type = cat(1, {'eog'}, cfg.artfctdef.type(:));
  elseif strcmp(cfg.rejecteog,'no') && ~isempty(dum)
    % this overrules the other setting, remove it from the type-list
    cfg.artfctdef.type(dum) = [];
  end
  cfg = rmfield(cfg, 'rejecteog');
end

% support the rejectXXX cfg settings for backward compatibility
if isfield(cfg, 'rejectjump')
  dum = strmatch('jump', cfg.artfctdef.type, 'exact');
  if strcmp(cfg.rejectjump,'yes') && isempty(dum)
    % this overrules the other setting, add it to the type-list
    cfg.artfctdef.type = cat(1, {'jump'}, cfg.artfctdef.type(:));
  elseif strcmp(cfg.rejectjump,'no') && ~isempty(dum)
    % this overrules the other setting, remove it from the type-list
    cfg.artfctdef.type(dum) = [];
  end
  cfg = rmfield(cfg, 'rejectjump');
end

% support the rejectXXX cfg settings for backward compatibility
if isfield(cfg, 'rejectfile')
  % this is slightly different to the ones above, since rejectfile is either 'no' or contains the filename
  dum = strmatch('file', cfg.artfctdef.type, 'exact');
  if ~strcmp(cfg.rejectfile,'no') && isempty(dum)
    % this overrules the other setting, add it to the type-list
    cfg.artfctdef.type = cat(1, {'file'}, cfg.artfctdef.type(:));
  elseif strcmp(cfg.rejectfile,'no') && ~isempty(dum)
    % this overrules the other setting, remove it from the type-list
    cfg.artfctdef.type(dum) = [];
  end
end

if nargin>1
  trl = findcfg(data.cfg, 'trl');
elseif isfield(cfg, 'trl')
  trl = cfg.trl;
end

% ensure that there are trials that can be scanned for artifacts and/or rejected
if isempty(trl)
  error('no trials were selected, cannot perform artifact detection/rejection');
end

% prevent double occurences of artifact types, ensure that the order remains the same
[dum, i] = unique(cfg.artfctdef.type);
cfg.artfctdef.type = cfg.artfctdef.type(sort(i));
% ensure that it is a row vector
cfg.artfctdef.type = cfg.artfctdef.type(:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call the appropriate function for each of the artifact types
% this will produce a Nx2 matrix with the begin and end sample of artifacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for type=1:length(cfg.artfctdef.type)
  fprintf('evaluating artifact_%s\n', cfg.artfctdef.type{type});
  % each call to artifact_xxx adds cfg.artfctdef.xxx.artifact
  cfg = feval(sprintf('artifact_%s', cfg.artfctdef.type{type}), cfg);
end

% collect the artifacts that have been detected from cfg.artfctdef.xxx.artifact
dum = fieldnames(cfg.artfctdef);
sel = [];
artifact = {};
for i=1:length(dum)
  sel(i) = issubfield(cfg.artfctdef, dum{i}) && issubfield(cfg.artfctdef, [dum{i} '.artifact']);
  if sel(i)
    artifact{end+1} = getsubfield(cfg.artfctdef, [dum{i} '.artifact']);
    num = size(artifact{end}, 1);
    if isempty(num)
      num = 0;
    end
    fprintf('detected %3d %s artifacts\n', num, dum{i});
  end
end
% update the configuration to reflect the artifacts types that were scanned
cfg.artfctdef.type = dum(find(sel));

% combine all trials into a single boolean vector
trialall = zeros(1,max(trl(:,2)));
for j=1:size(trl,1)
  trialall(trl(j,1):trl(j,2)) = 1;
end

% combine all artifacts into a single boolean vector
rejectall = zeros(1,max(trl(:,2)));
for i=1:length(cfg.artfctdef.type)
  dum = artifact{i};
  for j=1:size(dum,1)
    rejectall(dum(j,1):dum(j,2)) = i;  % the artifact type is coded here
  end
end

% ensure that both vectors are the same length
if length(trialall)>length(rejectall)
  rejectall(length(trialall)) = 0;
elseif length(trialall)<length(rejectall)
  trialall(length(rejectall)) = 0;
end

% make header, needed only for sampling frequency
if nargin ==1
    hdr = read_header(cfg.headerfile);
elseif nargin ==2
    hdr = fetch_header(data);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% give visual feedback on the trial definition and the artifacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(cfg.artfctdef.feedback, 'yes')
  COLOR = 'grcmykgrcmykgrcmykgrcmyk';
  % use the trial definition present in the local configuration
  trl  = trl;
  % compute the time axis that corresponds with each trial
  time = {};
  for i=1:size(trl,1)
    time{i} = offset2time(trl(i,3), hdr.Fs, trl(i,2)-trl(i,1)+1);
  end

  figure
  title('linear display of the continuous data')
  xlabel('sample number');
  hold on
  x = 1:length(trialall);
  y = ones(size(trialall)) * 0.53;
  y(~trialall) = nan;
  plot(x, y, 'b');
  for i=1:length(cfg.artfctdef.type)
    x = 1:length(rejectall);
    y = ones(size(rejectall)) * (0.5 - i*0.03);
    y(~(rejectall==i)) = nan;
    plot(x, y, COLOR(i));
  end
  set(gca, 'ytick', []);
  dum = (trialall~=0) | (rejectall~=0);
  [dum2, dum] = find(dum);
  smpbeg = dum(1)   - hdr.Fs;
  smpend = dum(end) + hdr.Fs;
  axis([smpbeg smpend 0 1]);
  legend({'defined trials', cfg.artfctdef.type{:}});

  figure
  title('individual trials after alignment')
  xlabel('time (s)');
  ylabel('trial number');
  hold on
  for i=1:size(trl,1)
    x = [time{i}(1) time{i}(end)];
    y = [i i];
    plot(x, y, 'b')
    for j=1:length(cfg.artfctdef.type)
      x = time{i};
      y = rejectall(trl(i,1):trl(i,2));
      y(y~=j) = nan;
      y(y==j) = i + j*0.03;
      plot(x, y, COLOR(j));
    end
    timebeg(i) = time{i}(1);
    timeend(i) = time{i}(end);
  end
  axis([min(timebeg)-0.1 max(timeend)+0.1 0.5 size(trl,1)+0.5]);
  axis ij
  legend({'defined trials', cfg.artfctdef.type{:}});
end % feedback

% convert to logical, NOTE: this is required for the following code
rejectall = (rejectall~=0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write the rejection to an EEP format file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(cfg.artfctdef, 'writerej') && ~isempty(cfg.artfctdef.writerej)
  fid = fopen(cfg.artfctdef.writerej, 'w');
  if fid<0
    error('could not open rejection file for writing');
  else
    % determine the begin and end of each rejected period (in samples)
    rejectonset = find(diff([0 rejectall])== 1);
    rejectofset = find(diff([rejectall 0])==-1);
    % determine the begin and end of each rejected period (in seconds)
    rejectonset = (rejectonset-1)/hdr.Fs;
    rejectofset = (rejectofset-1)/hdr.Fs;
    for rejlop=1:length(rejectonset)
      fprintf(fid, '%f-%f\n', rejectonset(rejlop), rejectofset(rejlop));
    end
    fclose(fid);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove the trials that (partially) coincide with a rejection mark
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(cfg.artfctdef.reject, 'partial') || strcmp(cfg.artfctdef.reject, 'complete')
  trl = trl;
  trialok = [];
  count_complete_reject = 0;
  count_partial_reject  = 0;
  for trial=1:size(trl,1)
    % cut out the part of the rejection-axis corresponding with this trial
    rejecttrial = rejectall(trl(trial,1):trl(trial,2));
    if all(not(rejecttrial))
      % the whole trial is good
      trialok = [trialok; trl(trial,:)];
    elseif all(rejecttrial)
      % the whole trial is bad
      count_complete_reject = count_complete_reject + 1;
    elseif any(rejecttrial) && strcmp(cfg.artfctdef.reject, 'complete')
      % some part of the trial is bad, reject the whole trial
      count_complete_reject = count_complete_reject + 1;
      continue;
    elseif any(rejecttrial) && strcmp(cfg.artfctdef.reject, 'partial')
      % some part of the trial is bad, reject only the bad part
      trialnew = [];
      rejecttrial = [0 not(rejecttrial) 0];
      % the first column is the begin sample, the second the end sample and the third is the offset
      trialnew(:,1) = find(diff(rejecttrial(:))== 1) + trl(trial,1) - 1;
      trialnew(:,2) = find(diff(rejecttrial(:))==-1) + trl(trial,1) - 2;
      trialnew(:,3) = find(diff(rejecttrial(:))== 1) - 1 + trl(trial,3);
      % some people use additional columns in the trl matrix to store trigger values and/or reaction times
      % these should remain linked to the original trial, i.e. they should be copied for each new fragment
      for i=4:size(trl,2)
        trialnew(:,i) = trl(trial,i);
      end
      minacceptnumsmp = round(cfg.artfctdef.minaccepttim .* hdr.Fs);
      trialnew(find(trialnew(:,2)-trialnew(:,1)<minacceptnumsmp),:) = [];
      count_partial_reject = count_partial_reject + 1;
      trialok = [trialok; trialnew];
    end
  end

  fprintf('rejected  %3d trials completely\n', count_complete_reject);
  fprintf('rejected  %3d trials partially\n', count_partial_reject);
  fprintf('resulting %3d trials\n', size(trialok,1));
  cfg.trlold = trl;      % return the original trial definition in the configuration
  cfg.trl    = trialok;  % return the cleaned trial definition in the configuration

else
  fprintf('not rejecting any data, only marking the artifacts\n');
end

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add version information to the artfctdef substructure
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: rejectartifact.m,v 1.48 2009/06/23 18:33:17 roboos Exp $';

% % remember the exact configuration details in the output
% cfgtmp = cfg;
% cfg = [];
% try cfg.trl        = cfgtmp.trl;        end
% try cfg.dataset    = cfgtmp.dataset;    end
% try cfg.datafile   = cfgtmp.datafile;   end
% try cfg.headerfile = cfgtmp.headerfile; end
% try cfg.continuous = cfgtmp.continuous; end
% cfg.previous = cfgtmp;

% apply the updated trial definition on the data
if nargin>1
    if isempty(cfg.trl)
        error('No trials left after artifact rejection.')
    else
        tmpcfg     = [];
        tmpcfg.trl = cfg.trl;
        data       = redefinetrial(tmpcfg,data);
        % remember the configuration details, this overwrites the stored configuration of redefinetrial
        data.cfg = cfg;
        % return the data instead of the cfg
        cfg = data;
    end
end

