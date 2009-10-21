function [data] = appenddata(cfg, varargin);

% APPENDDATA combines multiple datasets that have been preprocessed separately 
% into a single large dataset.
%
% Use as
%   data = appenddata(cfg, data1, data2, data3, ...)
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
% See also PREPROCESSING

% undocumented options:
%   none

% Copyright (C) 2005-2008, Robert Oostenveld
% Copyright (C) 2009, Jan-Mathijs Schoffelen
%
% $Log: appenddata.m,v $
% Revision 1.17  2009/04/30 14:42:20  jansch
% included explicit check on the channel order in the inputs. if the order is
% not the same the channels will be reordered according to the first input.
% included possibility to concatenate in the trial dimension when the number of
% channels is not entirely consistent (missing channels). updated documentation
%
% Revision 1.16  2008/09/22 20:17:42  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.15  2007/11/13 15:25:07  roboos
% removed handling of spike data, since that can be done in appendspike
%
% Revision 1.14  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.13  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.12  2007/01/09 09:50:39  roboos
% some changes in whitespace
%
% Revision 1.11  2007/01/04 17:06:33  roboos
% implemented support for appending spike data (i.e. timestamps) to a continuously sampled data set (i.e. one obtained from PREPROCESSING)
%
% Revision 1.10  2006/07/24 11:29:29  roboos
% use private/findcfg function for locating the trl and event in the nested (previous) cfgs
%
% Revision 1.9  2006/05/02 19:17:04  roboos
% search for th trl in cfg and cfg.previous etc., and also append it
%
% Revision 1.8  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.7  2006/01/30 12:15:01  roboos
% ensure consistent function declaration: "function [x] = funname()"
%
% Revision 1.6  2005/08/05 09:16:22  roboos
% removed the obsolete data.offset
%
% Revision 1.5  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.4  2005/03/14 11:17:29  jansch
% fixed small bug in handling of dimensions in the case of concatenating trials
%
% Revision 1.3  2005/02/16 17:07:13  roboos
% fixed error that occurred when input data label was not a column-cell-array
%
% Revision 1.2  2005/01/27 12:06:33  roboos
% fixed bug in checking of identical time axes
%
% Revision 1.1  2005/01/17 14:39:43  roboos
% new implementation
%

fieldtripdefs

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = checkdata(varargin{i}, 'datatype', 'raw', 'feedback', 'no');
end

% set the defaults
cfg = [];

Ndata = nargin-1;
if Ndata<2
  error('you must give at least two datasets to append');
end

% determine the dimensions of the data
Nchan  = zeros(1,Ndata);
Ntrial = zeros(1,Ndata);
label  = {};
for i=1:Ndata
  Nchan(i) = length(varargin{i}.label);
  Ntrial(i) = length(varargin{i}.trial);
  fprintf('input dataset %d, %d channels, %d trials\n', i, Nchan(i), Ntrial(i));
  label = [label(:); varargin{i}.label(:)];
end

% try to locate the trial definition (trl) in the nested configuration
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
end

% check the consistency of the labels across the input-structures
[alllabel, indx1, indx2] = unique(label, 'first');
order    = zeros(length(alllabel),Ndata);
for i=1:length(alllabel)
  for j=1:Ndata
    tmp = strmatch(alllabel{i}, varargin{j}.label, 'exact');
    if ~isempty(tmp)
      order(i,j) = tmp;
    end
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

if cattrial && catlabel
  error('cannot determine how the data should be concatenated');
  %FIXME think whether this can ever happen
elseif cattrial
  % concatenate the trials
  fprintf('concatenating the trials over all datasets\n');
  data = varargin{1};
  data.trial  = {};
  data.time   = {};
  for i=1:Ndata
    data.trial  = cat(2, data.trial,  varargin{i}.trial(:)');
    data.time   = cat(2, data.time,   varargin{i}.time(:)');
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
  for i=2:Ndata
    for j=1:Ntrial
      if ~all(data.time{j}==varargin{i}.time{j})
        error('there is a difference in the time axes of the input data');
      end
      data.trial{j} = [data.trial{j}; varargin{i}.trial{j}];
    end
    data.label = [data.label(:); varargin{i}.label(:)];
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


% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: appenddata.m,v 1.17 2009/04/30 14:42:20 jansch Exp $';
% remember the configuration details of the input data
cfg.previous = [];
for i=1:Ndata
  try, cfg.previous{i} = varargin{i}.cfg; end
end
% remember the exact configuration details in the output 
data.cfg = cfg;

fprintf('output dataset, %d channels, %d trials\n', length(data.label), length(data.trial));

