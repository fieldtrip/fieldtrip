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
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure. The data structure in the input file should be a
% cell array for this particular function.
%
% See also FT_PREPROCESSING, FT_APPENDFREQ

% Copyright (C) 2005-2008, Robert Oostenveld
% Copyright (C) 2009-2011, Jan-Mathijs Schoffelen
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar varargin
ft_preamble provenance varargin
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'no');
end

% determine the dimensions of the data
Ndata = length(varargin);

if Ndata<2
  error('you must give at least two datasets to append');
end

% ensure consistent input data
for i=2:Ndata
  if isfield(varargin{1}, 'topo')
    assert(isequaln(varargin{1}.topo, varargin{i}.topo), 'the input has inconsistent topo fields')
  end
  if isfield(varargin{1}, 'topolabel')
    assert(isequaln(varargin{1}.topolabel, varargin{i}.topolabel), 'the input has inconsistent topolabel fields')
  end
  if isfield(varargin{1}, 'unmixing')
    assert(isequaln(varargin{1}.unmixing, varargin{i}.unmixing), 'the input has inconsistent unmixing fields')
  end
end

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
% this is DEPRECATED - don't look in cfg-tree for stuff anymore
% hastrialinfo = 0;
% trl = cell(1, Ndata);
% for i=1:Ndata
%   if isfield(varargin{i}, 'cfg')
%     trl{i} = ft_findcfg(varargin{i}.cfg, 'trl');
%   else
%     trl{i} = [];
%   end
%   if isempty(trl{i})
%     % a trial definition is expected in each continuous data set
%     warning('could not locate the trial definition ''trl'' in data structure %d', i);
%   end
%   hastrialinfo = isfield(varargin{i}, 'trialinfo') + hastrialinfo;
% end
% hastrialinfo = hastrialinfo==Ndata;

hastrialinfo = 0;
hassampleinfo = 0;
sampleinfo = cell(1, Ndata);
for i=1:Ndata
  if isfield(varargin{i}, 'sampleinfo')
    sampleinfo{i} = varargin{i}.sampleinfo;
  else
    sampleinfo{i} = [];
  end

  % the function should behave properly even if no sampleinfo is present,
  % hence the warning seems inappropriate (ES, 24-apr-2014)
%   if isempty(sampleinfo{i})
%     % a sample definition is expected in each data set
%     warning('no ''sampleinfo'' field in data structure %d', i);
%   end

  hassampleinfo = isfield(varargin{i}, 'sampleinfo') + hassampleinfo;
  hastrialinfo = isfield(varargin{i}, 'trialinfo') + hastrialinfo;
end
hassampleinfo = hassampleinfo==Ndata;
hastrialinfo = hastrialinfo==Ndata;

% check the consistency of the labels across the input-structures
alllabel = unique(label, 'first');
order    = zeros(length(alllabel),Ndata);
for j=1:Ndata
  tmplabel = varargin{j}.label;
  [ix,iy]  = match_str(alllabel, tmplabel);
  order(ix,j) = iy;
end

% check consistency of sensor positions across inputs
haselec = 1;
hasgrad = 1;
for j=1:Ndata
  haselec = isfield(varargin{j}, 'elec') && haselec;
  hasgrad = isfield(varargin{j}, 'grad') && hasgrad;
  hasopto = isfield(varargin{j}, 'opto') && hasopto;
end

removesens = 0;
if haselec || hasgrad || hasopto
  sens = cell(1, Ndata);
  for j=1:Ndata
    if haselec, sens{j} = varargin{j}.elec; end
    if hasgrad, sens{j} = varargin{j}.grad; end
    if hasopto, sens{j} = varargin{j}.opto; end
    if j>1
      if ~isequaln(sens{j}, sens{1})
        removesens = 1;
        warning('sensor information does not seem to be consistent across the input arguments');
        break;
      end
    end
  end
end

% check whether the data are obtained from the same datafile in case either
% (1) we have sampleinfos and they are not identical or (2) we don't have
% sampleinfos
removesampleinfo = 0;
removetrialinfo  = 0;
try
  origfile1 = ft_findcfg(varargin{1}.cfg, 'datafile');
  for j=2:Ndata
    hassampleinfos = isfield(varargin{1}, 'sampleinfo') &&...
      isfield(varargin{j}, 'sampleinfo');

    if ((hassampleinfos &&...
        ~isequal(varargin{1}.sampleinfo, varargin{j}.sampleinfo)) ||...
        ~hassampleinfos) &&...
        ~isempty(origfile1) && ~strcmp(origfile1, ft_findcfg(varargin{j}.cfg, 'datafile'))
        removesampleinfo = 1;
        warning('input data comes from different datafiles; removing sampleinfo field');
        break;
    end
  end
catch err
  if strcmp(err.identifier, 'MATLAB:nonExistentField')
    % this means no data.cfg is present; should not be treated as a fatal error
    fprintf('cannot determine from which datafiles the data is taken\n');
  else
    % not sure which error, probably a bigger problem
    throw(err);
  end
end

catlabel   = all(sum(order~=0,2)==1);
cattrial   = any(sum(order~=0,2)==Ndata);
shuflabel  = cattrial && ~all(all(order-repmat(order(:,1),[1 Ndata])==0));
prunelabel = cattrial && sum(sum(order~=0,2)==Ndata)<length(alllabel);

if shuflabel
  fprintf('the channel order in the input-structures is not consistent, reordering\n');
  if prunelabel
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

elseif cattrial
  fprintf('concatenating the trials over all datasets\n');

  data = [];
  data.label  = varargin{1}.label;
  data.trial  = {};
  data.time   = {};
  if hassampleinfo, data.sampleinfo = []; end
  if hastrialinfo,  data.trialinfo  = []; end;

  for i=1:Ndata
    data.trial    = cat(2, data.trial,  varargin{i}.trial(:)');
    data.time     = cat(2, data.time,   varargin{i}.time(:)');
    % check if all datasets to merge have the sampleinfo field
    if hassampleinfo, data.sampleinfo = cat(1, data.sampleinfo, varargin{i}.sampleinfo); end
    if hastrialinfo,  data.trialinfo  = cat(1, data.trialinfo,  varargin{i}.trialinfo);  end
    % FIXME is not entirely robust if the different inputs have different number of columns in trialinfo
  end

elseif catlabel
  fprintf('concatenating the channels within each trial\n');

  if ~all(diff(Ntrial)==0)
    error('not all datasets have the same number of trials');
  else
    Ntrial = Ntrial(1);
  end

  data = [];
  data.label = varargin{1}.label;
  data.trial = varargin{1}.trial;
  data.time  = varargin{1}.time;
  if hassampleinfo, data.sampleinfo=varargin{i}.sampleinfo; end
  if hastrialinfo,  data.trialinfo =varargin{i}.trialinfo;  end

  for i=2:Ndata
    % concatenate the labels
    data.label = cat(1, data.label(:), varargin{i}.label(:));

    % check whether the trialinfo and sampleinfo fields are consistent
    if hassampleinfo && ~isequaln(data.sampleinfo, varargin{i}.sampleinfo)
      removesampleinfo = 1;
    end
    if hastrialinfo && ~isequaln(data.trialinfo, varargin{i}.trialinfo)
      removetrialinfo = 1;
    end
  end

  if ~isfield(data, 'fsample')
    fsample = 1/mean(diff(data.time{1}));
  else
    fsample = data.fsample;
  end

  for j=1:Ntrial
    %pre-allocate memory for this trial
    data.trial{j} = [data.trial{j}; zeros(sum(Nchan(2:end)), size(data.trial{j},2))];

    %fill this trial with data
    endchan = Nchan(1);
    %allow some jitter for irregular sample frequencies
    tolerance = 0.01*(1/fsample);
    for i=2:Ndata
      if ~all(data.time{j}-varargin{i}.time{j}<tolerance)
        error('there is a difference in the time axes of the input data');
      end
      begchan = endchan+1;
      endchan = endchan+Nchan(i);
      data.trial{j}(begchan:endchan,:) = varargin{i}.trial{j};
    end
  end

else
  % labels are inconsistent, cannot determine how to concatenate the data
  error('cannot determine how the data should be concatenated');
end

% some fields from the input should be copied over in the output
copyfield = {'grad', 'elec', 'topo', 'topolabel', 'unmixing', 'fsample'};
for i=1:length(copyfield)
  if isfield(varargin{1}, copyfield{i})
    data.(copyfield{i}) = varargin{1}.(copyfield{i});
  end
end

% unshuffle the channels again to match the order of the first input data-structure
if shuflabel
  fprintf('reordering the channels back to the original input order\n');
  [dum,reorder] = sort(order(order(:,1)~=0,1));
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
  fprintf('removing sampleinfo field from output\n');
  if isfield(data, 'sampleinfo'), data = rmfield(data, 'sampleinfo'); end
end

if removetrialinfo
  fprintf('removing trialinfo field from output\n');
  if isfield(data, 'trialinfo'), data = rmfield(data, 'trialinfo'); end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous varargin
ft_postamble history data
ft_postamble savevar data
