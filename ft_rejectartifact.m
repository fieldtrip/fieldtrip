function [cfg] = ft_rejectartifact(cfg, data)

% FT_REJECTARTIFACT removes data segments containing artifacts. It returns a
% configuration structure with a modified trial definition which can be used for
% preprocessing of only the clean data.
%
% You should start by detecting the artifacts in the data using the function
% FT_ARTIFACT_xxx where xxx is the type of artifact. Subsequently FT_REJECTARTIFACT
% looks at the detected artifacts and removes them from the trial definition or from
% the data. In case you wish to replace bad parts by nans, you have to specify data
% as an input parameter.
%
% Use as
%   cfg = ft_rejectartifact(cfg)
% with the cfg as obtained from FT_DEFINETRIAL, or as
%   data = ft_rejectartifact(cfg, data)
% with the data as obtained from FT_PREPROCESSING
%
% The following configuration options are supported:
%   cfg.artfctdef.reject          = 'none', 'partial', 'complete', 'nan', or 'value' (default = 'complete')
%   cfg.artfctdef.minaccepttim    = when using partial rejection, minimum length
%                                   in seconds of remaining trial (default = 0.1)
%   cfg.artfctdef.crittoilim      = when using complete rejection, reject trial only when artifacts occur within
%                                   this time window (default = whole trial). This only works with in-memory data,
%                                   since trial time axes are unknown for data on disk.
%   cfg.artfctdef.feedback        = 'yes' or 'no' (default = 'no')
%   cfg.artfctdef.invert          = 'yes' or 'no' (default = 'no')
%   cfg.artfctdef.value           = scalar value to replace the data in the artifact segments (default = nan)
%   cfg.artfctdef.eog.artifact    = Nx2 matrix with artifact segments, this is added to the cfg by using FT_ARTIFACT_EOG
%   cfg.artfctdef.jump.artifact   = Nx2 matrix with artifact segments, this is added to the cfg by using FT_ARTIFACT_JUMP
%   cfg.artfctdef.muscle.artifact = Nx2 matrix with artifact segments, this is added to the cfg by using FT_ARTIFACT_MUSCLE
%   cfg.artfctdef.zvalue.artifact = Nx2 matrix with artifact segments, this is added to the cfg by using FT_ARTIFACT_ZVALUE
%   cfg.artfctdef.visual.artifact = Nx2 matrix with artifact segments, this is added to the cfg by using FT_DATABROWSER
%   cfg.artfctdef.xxx.artifact    = Nx2 matrix with artifact segments, this could be added by your own artifact detection function
%
% A trial that contains an artifact can be rejected completely or partially. In case
% of partial rejection, a minimum length of the resulting sub-trials can be
% specified using minaccepttim.
%
% Output:
%   If cfg is used as the only input parameter, the output is a cfg structure with an updated trl.
%   If cfg and data are both input parameters, the output is an updated raw data structure with only the clean data segments.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_ARTIFACT_EOG, FT_ARTIFACT_MUSCLE, FT_ARTIFACT_JUMP, FT_ARTIFACT_THRESHOLD,
% FT_ARTIFACT_CLIP, FT_ARTIFACT_ECG, FT_DATABROWSER, FT_REJECTVISUAL

% Undocumented local options:
% cfg.headerfile
% cfg.reject
% cfg.trl
% cfg.trlold
% cfg.version
%
% These old configuration options are still supported
% cfg.rejectmuscle      = 'no' or 'yes'
% cfg.rejecteog         = 'no' or 'yes'
% cfg.rejectjump        = 'no' or 'yes'
% cfg.rejectfile        = string with filename
% cfg.artfctdef.writerej = filename of rejection file
% cfg.artfctdef.type    = cell-array with strings, e.g. {'eog', 'muscle' 'jump'}

% Copyright (C) 2003-2018, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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

% FIXME this function contains a lot of lines of code that pertain to backward
% compatibility support that dates back to 2004/2005. It would be good to strip
% that code and only keep the relevant parts

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% ft_checkdata is done further down

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');

% set the defaults
cfg.artfctdef              = ft_getopt(cfg, 'artfctdef');
cfg.artfctdef.type         = ft_getopt(cfg.artfctdef, 'type', {});
cfg.artfctdef.reject       = ft_getopt(cfg.artfctdef, 'reject', 'complete');
cfg.artfctdef.minaccepttim = ft_getopt(cfg.artfctdef, 'minaccepttim', 0.1);
cfg.artfctdef.crittoilim   = ft_getopt(cfg.artfctdef, 'crittoilim', []);
cfg.artfctdef.feedback     = ft_getopt(cfg.artfctdef, 'feedback', 'no');
cfg.artfctdef.invert       = ft_getopt(cfg.artfctdef, 'invert', 'no');
cfg.artfctdef.value        = ft_getopt(cfg.artfctdef, 'value', nan);

% convert from old-style to new-style configuration
if isfield(cfg,'reject')
  ft_warning('converting from old-style artifact configuration to new-style');
  cfg.artfctdef.reject = cfg.reject;
  cfg = rmfield(cfg, 'reject');
end

% convert from old-style to new-style configuration
if isfield(cfg.artfctdef,'common')
  ft_warning('converting from old-style artifact configuration to new-style');
  if isfield(cfg.artfctdef.common,'minaccepttim')
    cfg.artfctdef.minaccepttim = cfg.artfctdef.common.minaccepttim;
    cfg.artfctdef = rmfield(cfg.artfctdef, 'common');
  end
end

% ensure that it is a cell-array
if ischar(cfg.artfctdef.type)
  cfg.artfctdef.type = {cfg.artfctdef.type};
end

% support the rejectXXX cfg settings for backward compatibility
if isfield(cfg, 'rejectmuscle')
  dum = find(strcmp('muscle', cfg.artfctdef.type));
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
  dum = find(strcmp('eog', cfg.artfctdef.type));
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
  dum = find(strcmp('jump', cfg.artfctdef.type));
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
  dum = find(strcmp('file', cfg.artfctdef.type));
  if ~strcmp(cfg.rejectfile,'no') && isempty(dum)
    % this overrules the other setting, add it to the type-list
    cfg.artfctdef.type = cat(1, {'file'}, cfg.artfctdef.type(:));
  elseif strcmp(cfg.rejectfile,'no') && ~isempty(dum)
    % this overrules the other setting, remove it from the type-list
    cfg.artfctdef.type(dum) = [];
  end
end

% the data can be specified as input variable or through cfg.inputfile
hasdata = exist('data', 'var');

if hasdata
  % check if the input data is valid for this function
  data = ft_checkdata(data, 'hassampleinfo', 'yes');
  
  if isfield(data, 'sampleinfo')
    trl = zeros(numel(data.trial), 3);
    trl(:,[1 2]) = data.sampleinfo;
    
    % recreate offset vector (artifact functions depend on this)
    % TODO: the artifact rejection stuff should be rewritten to avoid
    % needing this workaround
    for ntrl = 1:numel(data.trial)
      trl(ntrl,3) = time2offset(data.time{ntrl}, data.fsample);
    end
    
    if isfield(data, 'trialinfo')
      if istable(data.trialinfo)
        % convert table into normal array, keep the column labels
        VariableNames = data.trialinfo.Properties.VariableNames;
        data.trialinfo = table2array(data.trialinfo);
      end
      trl = [trl data.trialinfo];
    end
  else
    trl = [];
  end
  
elseif isfield(cfg, 'trl')
  trl = cfg.trl;
end

% ensure the crittoilim input argument is valid
if ~isempty(cfg.artfctdef.crittoilim)
  
  if (size(cfg.artfctdef.crittoilim,2) ~= 2 ...
      || (size(cfg.artfctdef.crittoilim,1) ~= size(trl,1) ...
      && size(cfg.artfctdef.crittoilim,1) ~= 1))
    ft_error('if specified, cfg.artfctdef.crittoilim should be a 1x2 or Nx2 vector');
  end
  
  % if specified as 1x2 vector, expand into Nx2
  if (size(cfg.artfctdef.crittoilim,1) == 1)
    cfg.artfctdef.crittoilim = repmat(cfg.artfctdef.crittoilim,size(trl,1),1);
  end
  
  checkCritToi = 1; % flag for convenience
else
  checkCritToi = 0;
end

% ensure that there are trials that can be scanned for artifacts and/or rejected
if isempty(trl)
  ft_error('no trials were selected, cannot perform artifact detection/rejection');
end

% prevent double occurences of artifact types, ensure that the order remains the same
[dum, i] = unique(cfg.artfctdef.type);
cfg.artfctdef.type = cfg.artfctdef.type(sort(i));
% ensure that it is a row vector
cfg.artfctdef.type = cfg.artfctdef.type(:)';

% If bad parts are to be filled with nans, make sure data is available
if ~hasdata && (strcmp(cfg.artfctdef.reject, 'nan') || strcmp(cfg.artfctdef.reject, 'value'))
  ft_error('If bad parts are to be filled with nans or another value, the input data has to be specified');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call the appropriate function for each of the artifact types
% this will produce a Nx2 matrix with the begin and end sample of artifacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for type=1:length(cfg.artfctdef.type)
  funhandle = ft_getuserfun(cfg.artfctdef.type{type}, 'artifact');
  fprintf('evaluating %s\n', func2str(funhandle));
  % each call to artifact_xxx adds cfg.artfctdef.xxx.artifact
  if hasdata
    cfg = feval(funhandle, cfg, data);
  else
    cfg = feval(funhandle, cfg);
  end
end

% collect the artifacts that have been detected from cfg.artfctdef.xxx.artifact
dum = fieldnames(cfg.artfctdef);
sel = false(size(dum));
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
cfg.artfctdef.type = dum(sel);

% combine all trials into a single boolean vector
trialall = convert_event(trl(:,1:2), 'boolvec');

% combine all artifacts into a single vector
rejectall = zeros(1, max(trl(:,2)));
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
if hasdata
  hdr = ft_fetch_header(data);
else
  if isfield(cfg, 'headerformat')
    hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
  else
    hdr = ft_read_header(cfg.headerfile);
  end
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

% convert to logical, this is required for the subsequent code
rejectall = (rejectall~=0);

% invert the artifact selection
if istrue(cfg.artfctdef.invert)
  fprintf('inverting selection of clean/artifactual data\n');
  rejectall = ~rejectall;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write the rejection to an EEP format file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(cfg.artfctdef, 'writerej') && ~isempty(cfg.artfctdef.writerej)
  fid = fopen(cfg.artfctdef.writerej, 'w');
  if fid<0
    ft_error('could not open rejection file for writing');
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
if any(strcmp(cfg.artfctdef.reject, {'partial', 'complete', 'nan', 'value'}))
  trialok = [];
  
  count_complete_reject = 0;
  count_partial_reject  = 0;
  count_nan             = 0;
  count_value           = 0;
  count_outsidecrit     = 0;
  
  trlCompletelyRemovedInd = [];
  trlPartiallyRemovedInd  = [];
  
  for trial=1:size(trl,1)
    % cut out the part of the rejection-axis corresponding with this trial
    rejecttrial = rejectall(trl(trial,1):trl(trial,2));
    
    if all(not(rejecttrial))
      % the whole trial is good
      trialok = [trialok; trl(trial,:)];
      
    elseif all(rejecttrial) && strcmp(cfg.artfctdef.reject, 'nan')
      % the whole trial is bad, but it is requested to be replaced with nans
      data.trial{trial}(:,rejecttrial) = nan;
      count_nan = count_nan + 1;
      trialok = [trialok; trl(trial,:)]; % Mark the trial as good as nothing will be removed
      
    elseif all(rejecttrial) && strcmp(cfg.artfctdef.reject, 'value')
      % the whole trial is bad, but it is requested to be replaced with a specific value
      data.trial{trial}(:,rejecttrial) = cfg.artfctdef.value;
      count_value = count_value + 1;
      trialok = [trialok; trl(trial,:)]; % Mark the trial as good as nothing will be removed

    elseif all(rejecttrial)
      % the whole trial is bad
      count_complete_reject = count_complete_reject + 1;
      trlCompletelyRemovedInd = [trlCompletelyRemovedInd trial];
      
    elseif any(rejecttrial) && strcmp(cfg.artfctdef.reject, 'complete')
      % some part of the trial is bad, check if within crittoilim?
      if (checkCritToi)
        critInd = (data.time{trial} >= cfg.artfctdef.crittoilim(trial,1) ...
          & data.time{trial} <= cfg.artfctdef.crittoilim(trial,2));
        if (any(critInd & rejecttrial))
          count_complete_reject = count_complete_reject + 1;
          trlCompletelyRemovedInd = [trlCompletelyRemovedInd trial];
          continue;
        else
          trialok = [trialok; trl(trial,:)];
          count_outsidecrit = count_outsidecrit + 1;
        end
      else % no crittoilim checking required
        count_complete_reject = count_complete_reject + 1;
        trlCompletelyRemovedInd = [trlCompletelyRemovedInd trial];
        continue;
      end
      
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
      trialnew((trialnew(:,2)-trialnew(:,1))<minacceptnumsmp,:) = [];
      count_partial_reject = count_partial_reject + 1;
      trialok = [trialok; trialnew];
      trlPartiallyRemovedInd = [trlPartiallyRemovedInd trial];
      
    elseif any(rejecttrial) && strcmp(cfg.artfctdef.reject, 'nan')
      % Some part of the trial is bad, replace bad part with nans
      data.trial{trial}(:,rejecttrial) = nan;
      count_nan = count_nan + 1;
      trialok = [trialok; trl(trial,:)]; % Mark the trial as good as nothing will be removed
   
    elseif any(rejecttrial) && strcmp(cfg.artfctdef.reject, 'value')
      % Some part of the trial is bad, replace bad part with specified value
      data.trial{trial}(:,rejecttrial) = cfg.artfctdef.value;
      count_value = count_value + 1;
      trialok = [trialok; trl(trial,:)]; % Mark the trial as good as nothing will be removed

    end
  end % for each trial
  
  fprintf('rejected  %3d trials completely\n', count_complete_reject);
  fprintf('rejected  %3d trials partially\n', count_partial_reject);
  fprintf('filled parts of  %3d trials with nans\n', count_nan);
  fprintf('filled parts of  %3d trials with the specified value\n', count_value);
  if (checkCritToi)
    fprintf('retained  %3d trials with artifacts outside critical window\n', count_outsidecrit);
  end
  fprintf('resulting %3d trials\n', size(trialok,1));
  cfg.trlold = trl;      % return the original trial definition in the configuration
  cfg.trl    = trialok;  % return the cleaned trial definition in the configuration
  
  if strcmp(cfg.artfctdef.feedback, 'yes')
    fprintf('the following trials were completely removed: ');
    for k = trlCompletelyRemovedInd
      fprintf('%d ', k);
    end
    fprintf('\nthe following trials were partially removed: ');
    for k = trlPartiallyRemovedInd
      fprintf('%d ', k);
    end
    fprintf('\n');
  end
  
else
  fprintf('not rejecting any data, only marking the artifacts\n');
end

if isempty(cfg.trl)
  ft_error('No trials left after artifact rejection.')
else
  if hasdata && ~strcmp(cfg.artfctdef.reject, 'nan') % Skip this step to avoid removing parts that should be filled with nans
    % apply the updated trial definition on the data
    tmpcfg      = keepfields(cfg, {'trl', 'showcallinfo'});
    data        = removefields(data, {'trialinfo'});
    data        = ft_redefinetrial(tmpcfg, data);
    % restore the provenance information
    [cfg, data] = rollback_provenance(cfg, data);
  end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
if hasdata
  ft_postamble previous data
  ft_postamble history data
  ft_postamble savevar data
  % return the data, the output variable is called cfg instead of data
  cfg = data;
end
