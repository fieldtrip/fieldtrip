function [data] = ft_appenddata(cfg, varargin)

% FT_APPENDDATA combines multiple datasets that have been preprocessed separately
% into a single large dataset.
%
% Use as
%   data = ft_appenddata(cfg, data1, data2, data3, ...)
% where the configuration can be empty.
%
% If the input datasets all have the same channels, the trials will be
% concatenated. This is useful for example if you have different
% experimental conditions, which, besides analyzing them separately, for
% some reason you also want to analyze together. The function will check
% for consistency in the order of the channels. If the order is inconsistent
% the channel order of the output will be according to the channel order of
% the first data structure in the input.
%
% If the input datasets have different channels, but the same number of
% trials, the channels will be concatenated within each trial. This is
% useful for example if the data that you want to analyze contains both
% MEG and EMG channels which require different preprocessing options.
%
% Occasionally, the data needs to be concatenated in the trial dimension while
% there's a slight discrepancy in the channels in the input data (e.g. missing
% channels in one of the data structures). The function will then return a data
% structure containing only the channels which are present in all inputs.
%
% See also FT_PREPROCESSING

% Undocumented local options:
%   cfg.inputfile  = one can specifiy preanalysed saved data as input
%                     The data should be provided in a cell array
%   cfg.outputfile = one can specify output as file to save to disk

% Copyright (C) 2005-2008, Robert Oostenveld
% Copyright (C) 2009, Jan-Mathijs Schoffelen
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

% set the defaults
if ~isfield(cfg, 'inputfile'),    cfg.inputfile  = [];          end
if ~isfield(cfg, 'outputfile'),   cfg.outputfile = [];          end

hasdata = nargin>1;
if ~isempty(cfg.inputfile) % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    for i=1:numel(cfg.inputfile)
      varargin{i} = loadvar(cfg.inputfile{i}, 'data'); % read datasets from array inputfile
      Ndata       = numel(cfg.inputfile); % use Ndata as if separate datafiles were specified
    end
  end
else Ndata = nargin-1; % use old implementation of Ndata if no inputfile is specified
  if Ndata<2
    error('you must give at least two datasets to append');
  end
end

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = checkdata(varargin{i}, 'datatype', 'raw', 'feedback', 'no', 'hastrialdef', 'yes');
end

% determine the dimensions of the data
Nchan  = zeros(1,Ndata);
Ntrial = zeros(1,Ndata);
label  = {};
for i=1:Ndata
  Nchan(i)  = length(varargin{i}.label);
  Ntrial(i) = length(varargin{i}.trial);
  fprintf('input dataset %d, %d channels, %d trials\n', i, Nchan(i), Ntrial(i));
  label = cat(1, label(:), varargin{i}.label(:));
end

% try to locate the trial definition (trl) in the nested configuration and
% check whether the input data contains trialinfo
hastrialinfo = 0;
for i=1:Ndata
  if isfield(varargin{i}, 'cfg')
    trl{i} = findcfg(varargin{i}.cfg, 'trl');
  else
    trl{i} = [];
  end
  if isempty(trl{i})
    % a trial definition is expected in each continuous data set
    warning(sprintf('could not locate the trial definition ''trl'' in data structure %d', i));
  end
  hastrialinfo = isfield(varargin{i}, 'trialinfo') + hastrialinfo;
end
hastrialinfo = hastrialinfo==Ndata;

hassampleinfo = 0;
for i=1:Ndata
  if isfield(varargin{i}, 'sampleinfo')
     sampleinfo{i} = varargin{i}.sampleinfo;
  else
     sampleinfo{i} = [];
  end
  if isempty(sampleinfo{i})
    % a sample definition is expected in each data set
    warning(sprintf('no ''sampleinfo'' field in data structure %d', i));
  end
  hassampleinfo = isfield(varargin{i}, 'sampleinfo') + hassampleinfo;
end
hassampleinfo = hassampleinfo==Ndata;

% check the consistency of the labels across the input-structures
[alllabel, indx1, indx2] = unique(label, 'first');
order    = zeros(length(alllabel),Ndata);
%for i=1:length(alllabel)
%  for j=1:Ndata
%    tmp = strmatch(alllabel{i}, varargin{j}.label, 'exact');
%    if ~isempty(tmp)
%      order(i,j) = tmp;
%    end
%  end
%end

% replace the nested for-loops with something faster
for j=1:Ndata
  tmplabel = varargin{j}.label;
  [ix,iy]  = match_str(alllabel, tmplabel);
  order(ix,j) = iy;
end

% check consistency of sensor positions across inputs
for j=1:Ndata
  haselec(j) = isfield(varargin{j}, 'elec');
  hasgrad(j) = isfield(varargin{j}, 'grad');
end
removesens = 0;
if all(haselec==1) || all(hasgrad==1),
  for j=1:Ndata
    if haselec, sens{j} = getfield(varargin{j}, 'elec'); end
    if hasgrad, sens{j} = getfield(varargin{j}, 'grad'); end
    if j>1,
      if numel(sens{j}.pnt) ~= numel(sens{1}.pnt) || any(sens{j}.pnt(:) ~= sens{1}.pnt(:)),
        removesens = 1;
        warning('sensor information does not seem to be consistent across the input arguments');
        break;
      end
    end
  end
elseif haselec(1)==1 || hasgrad(1)==1,
  % apparently the first input has a grad or elec, but not all the other
  % ones. because the output data will be initialized to be the first input
  % data structure, this should be removed
  removesens = 1;   
end

% check whether the data are obtained from the same datafile
origfile1      = findcfg(varargin{1}.cfg, 'datafile');
removesampleinfo = 0;
for j=2:Ndata
    if ~strcmp(origfile1, findcfg(varargin{j}.cfg, 'datafile')),
        removesampleinfo = 1;
        warning('input data comes from different datafiles');
        break;
    end
end

catlabel   = all(sum(order~=0,2)==1);
cattrial   = any(sum(order~=0,2)==Ndata);
shuflabel  = cattrial && ~all(all(order-repmat(order(:,1),[1 Ndata])==0));
prunelabel = cattrial && sum(sum(order~=0,2)==Ndata)<length(alllabel);

if shuflabel,
  fprintf('the channel order in the input-structures is not consistent, reordering\n');
  if prunelabel,
    fprintf('not all input-structures contain the same channels, pruning the input prior to concatenating over trials\n');
    selall    = find(sum(order~=0,2)==Ndata);
    alllabel  = alllabel(selall);
    order     = order(selall,:);
  end
  for i=1:Ndata
    varargin{i}.label = varargin{i}.label(order(:,i));
    for j=1:length(varargin{i}.trial)
      varargin{i}.trial{j} = varargin{i}.trial{j}(order(:,i),:);
    end
  end
end

% FIXME create the output from scratch and don't use the first varargin
% (both for cattrial and catlabel
if cattrial && catlabel
  error('cannot determine how the data should be concatenated');
  % FIXME think whether this can ever happen
  
elseif cattrial
  % concatenate the trials
  fprintf('concatenating the trials over all datasets\n');
  data = varargin{1};
  data.trial  = {};
  data.time   = {};
  data.sampleinfo = [];
  if hastrialinfo, data.trialinfo = []; end;
  for i=1:Ndata
    data.trial    = cat(2, data.trial,  varargin{i}.trial(:)');
    data.time     = cat(2, data.time,   varargin{i}.time(:)');
    % check if all datasets to merge have the sampleinfo field
    if hassampleinfo
      data.sampleinfo = cat(1, data.sampleinfo, varargin{i}.sampleinfo);
    else
      if isempty(sampleinfo{i})
        varargin{i}.sampleinfo = [];
      end
      data.sampleinfo = cat(1, data.sampleinfo, varargin{i}.sampleinfo);
    end
    if hastrialinfo, data.trialinfo = cat(1, data.trialinfo, varargin{i}.trialinfo); end;
    % FIXME is not entirely robust if the different inputs have different
    % number of columns in trialinfo
  end
  % also concatenate the trial specification
  cfg.trl = cat(1, trl{:});
  
elseif catlabel
  % concatenate the channels in each trial
  fprintf('concatenating the channels within each trial\n');
  data = varargin{1};
  if ~all(diff(Ntrial)==0)
    error('not all datasets have the same number of trials')
  else
    Ntrial = Ntrial(1);
  end
  Nch(1) = numel(data.label);
  for i=2:Ndata
    Nch(i,1)   = numel(varargin{i}.label);
    data.label = cat(1, data.label(:), varargin{i}.label(:));
  end
  
  for j=1:Ntrial
    %pre-allocate memory for this trial
    data.trial{j} = [data.trial{j}; zeros(sum(Nch(2:end)), size(data.trial{j},2))];
    
    %fill this trial with data
    endchan = Nch(1);
    for i=2:Ndata
      if ~all(data.time{j}==varargin{i}.time{j})
        error('there is a difference in the time axes of the input data');
      end
      begchan = endchan+1;
      endchan = endchan+Nch(i);
      data.trial{j}(begchan:endchan,:) = varargin{i}.trial{j};
    end
  end
  
else
  % labels are inconsistent, cannot determine how to concatenate the data
  error('cannot determine how the data should be concatenated');
end

% unshuffle the channels again to match the order of the first input data-structure
if shuflabel
  [srt,reorder] = sort(order(order(:,1)~=0,1));
  
  fprintf('reordering the channels\n');
  for i=1:length(data.trial)
    data.trial{i} = data.trial{i}(reorder,:);
  end
  data.label = data.label(reorder);
end

if removesens
  fprintf('removing sensor information from output\n');
  if haselec, data = rmfield(data, 'elec'); end
  if hasgrad, data = rmfield(data, 'grad'); end
end

if removesampleinfo
    fprintf('removing trial definition from output\n');
    data            = rmfield(data, 'sampleinfo');
    %cfg.trl(:, 1:2) = nan;
    if isfield(cfg, 'trl'), cfg = rmfield(cfg, 'trl'); end
end

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
% remember the configuration details of the input data
cfg.previous = [];
for i=1:Ndata
  try, cfg.previous{i} = varargin{i}.cfg; end
end

% remember the exact configuration details in the output
data.cfg = cfg;

fprintf('output dataset, %d channels, %d trials\n', length(data.label), length(data.trial));

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', data); % use the variable name "data" in the output file
end
