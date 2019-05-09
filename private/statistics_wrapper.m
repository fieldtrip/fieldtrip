function [stat, cfg] = statistics_wrapper(cfg, varargin)

% STATISTICS_WRAPPER performs the selection of the biological data for
% timelock, frequency or source data and sets up the design vector or
% matrix.
%
% The specific configuration options for selecting timelock, frequency
% of source data are described in FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS and
% FT_SOURCESTATISTICS, respectively.
%
% After selecting the data, control is passed on to a data-independent
% statistical subfunction that computes test statistics plus their associated
% significance probabilities and critical values under some null-hypothesis. The statistical
% subfunction that is called is FT_STATISTICS_xxx, where cfg.method='xxx'. At
% this moment, we have implemented two statistical subfunctions:
% FT_STATISTICS_ANALYTIC, which calculates analytic significance probabilities and critical
% values (exact or asymptotic), and FT_STATISTICS_MONTECARLO, which calculates
% Monte-Carlo approximations of the significance probabilities and critical values.
%
% The specific configuration options for the statistical test are
% described in FT_STATISTICS_xxx.

% This function depends on PREPARE_TIMEFREQ_DATA which has the following options:
% cfg.avgoverchan
% cfg.avgoverfreq
% cfg.avgovertime
% cfg.channel
% cfg.channelcmb
% cfg.datarepresentation (set in STATISTICS_WRAPPER cfg.datarepresentation = 'concatenated')
% cfg.frequency
% cfg.latency
% cfg.precision
% cfg.previous (set in STATISTICS_WRAPPER cfg.previous = [])
% cfg.version (id and name set in STATISTICS_WRAPPER)

% Copyright (C) 2005-2006, Robert Oostenveld
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

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed',     {'approach',   'method'});
cfg = ft_checkconfig(cfg, 'required',    {'method'});
cfg = ft_checkconfig(cfg, 'forbidden',   {'transform'});

% set the defaults
cfg.latency     = ft_getopt(cfg, 'latency', 'all');
cfg.frequency   = ft_getopt(cfg, 'frequency', 'all');
cfg.roi         = ft_getopt(cfg, 'roi',       []);
cfg.avgovertime = ft_getopt(cfg, 'avgovertime', 'no');
cfg.avgoverfreq = ft_getopt(cfg, 'avgoverfreq', 'no');
cfg.avgoverroi  = ft_getopt(cfg, 'avgoverroi',  'no');

% determine the type of the input and hence the output data
if ~exist('OCTAVE_VERSION', 'builtin')
  [s, i] = dbstack;
  if length(s)>1
    [caller_path, caller_name, caller_ext] = fileparts(s(2).name);
  else
    caller_path = '';
    caller_name = '';
    caller_ext  = '';
  end
  % evalin('caller', 'mfilename') does not work for Matlab 6.1 and 6.5
  istimelock = strcmp(caller_name,'ft_timelockstatistics');
  isfreq     = strcmp(caller_name,'ft_freqstatistics');
  issource   = strcmp(caller_name,'ft_sourcestatistics');
else
  % cannot determine the calling function in Octave, try looking at the
  % data instead
  istimelock  = isfield(varargin{1},'time') && ~isfield(varargin{1},'freq') && isfield(varargin{1},'avg');
  isfreq      = isfield(varargin{1},'time') && isfield(varargin{1},'freq');
  issource    = isfield(varargin{1},'pos') || isfield(varargin{1},'transform');
end

if (istimelock+isfreq+issource)~=1
  ft_error('Could not determine the type of the input data');
end

if istimelock || isfreq,
  % these defaults only apply to channel level data
  cfg.channel     = ft_getopt(cfg, 'channel',     'all');
  cfg.avgoverchan = ft_getopt(cfg, 'avgoverchan', 'no');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect the biological data (the dependent parameter)
%  and the experimental design (the independent parameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if issource

  % test that all source inputs have the same dimensions and are spatially aligned
  for i=2:length(varargin)
    if isfield(varargin{1}, 'dim') && (length(varargin{i}.dim)~=length(varargin{1}.dim) || ~all(varargin{i}.dim==varargin{1}.dim))
      ft_error('dimensions of the source reconstructions do not match, use FT_VOLUMENORMALISE first');
    end
    if isfield(varargin{1}, 'pos') && (length(varargin{i}.pos(:))~=length(varargin{1}.pos(:)) || ~all(varargin{i}.pos(:)==varargin{1}.pos(:)))
      ft_error('grid locations of the source reconstructions do not match, use FT_VOLUMENORMALISE first');
    end
    if isfield(varargin{1}, 'transform') && ~all(varargin{i}.transform(:)==varargin{1}.transform(:))
      ft_error('spatial coordinates of the source reconstructions do not match, use FT_VOLUMENORMALISE first');
    end
  end

  Nsource = length(varargin);
  Nvoxel  = length(varargin{1}.inside) + length(varargin{1}.outside);

  if ~isempty(cfg.roi)
    if ischar(cfg.roi)
      cfg.roi = {cfg.roi};
    end
    % the source representation should specify the position of each voxel in MNI coordinates
    x = varargin{1}.pos(:,1);  % this is from left (negative) to right (positive)
    % determine the mask to restrict the subsequent analysis
    % process each of the ROIs, and optionally also left and/or right separately
    roimask  = {};
    roilabel = {};
    for i=1:length(cfg.roi)
      if islogical(cfg.roi{i})
        tmp = cfg.roi{i};
      else
        tmpcfg.roi = cfg.roi{i};
        tmpcfg.inputcoord = cfg.inputcoord;
        tmpcfg.atlas = cfg.atlas;
        tmp = ft_volumelookup(tmpcfg, varargin{1});
      end
      if strcmp(cfg.avgoverroi, 'no') && ~isfield(cfg, 'hemisphere')
        % no reason to deal with separated left/right hemispheres
        cfg.hemisphere = 'combined';
      end

      if     strcmp(cfg.hemisphere, 'left')
        tmp(x>=0)    = 0;  % exclude the right hemisphere
        roimask{end+1}  = tmp;
        roilabel{end+1} = ['Left '  cfg.roi{i}];

      elseif strcmp(cfg.hemisphere, 'right')
        tmp(x<=0)    = 0;  % exclude the right hemisphere
        roimask{end+1}  = tmp;
        roilabel{end+1} = ['Right ' cfg.roi{i}];

      elseif strcmp(cfg.hemisphere, 'both')
        % deal separately with the voxels on the left and right side of the brain
        tmpL = tmp; tmpL(x>=0) = 0;  % exclude the right hemisphere
        tmpR = tmp; tmpR(x<=0) = 0;  % exclude the left hemisphere
        roimask{end+1}  = tmpL;
        roimask{end+1}  = tmpR;
        roilabel{end+1} = ['Left '  cfg.roi{i}];
        roilabel{end+1} = ['Right ' cfg.roi{i}];
        clear tmpL tmpR

      elseif strcmp(cfg.hemisphere, 'combined')
        % all voxels of the ROI can be combined
        roimask{end+1}  = tmp;
        if ischar(cfg.roi{i})
          roilabel{end+1} = cfg.roi{i};
        else
          roilabel{end+1} = ['ROI ' num2str(i)];
        end

      else
        ft_error('incorrect specification of cfg.hemisphere');
      end
      clear tmp
    end % for each roi

    % note that avgoverroi=yes is implemented differently at a later stage
    % avgoverroi=no is implemented using the inside/outside mask
    if strcmp(cfg.avgoverroi, 'no')
      for i=2:length(roimask)
        % combine them all in the first mask
        roimask{1} = roimask{1} | roimask{i};
      end
      roimask = roimask{1};  % only keep the combined mask
      % the source representation should have an inside and outside vector containing indices
      sel = find(~roimask);
      varargin{1}.inside  = setdiff(varargin{1}.inside, sel);
      varargin{1}.outside = union(varargin{1}.outside, sel);
      clear roimask roilabel
    end % if avgoverroi=no
  end

  % get the source parameter on which the statistic should be evaluated
  if strcmp(cfg.parameter, 'mom') && isfield(varargin{1}, 'avg') && isfield(varargin{1}.avg, 'csdlabel') && isfield(varargin{1}, 'cumtapcnt')
    [dat, cfg] = get_source_pcc_mom(cfg, varargin{:});
  elseif strcmp(cfg.parameter, 'mom') && isfield(varargin{1}, 'avg') && ~isfield(varargin{1}.avg, 'csdlabel')
    [dat, cfg] = get_source_lcmv_mom(cfg, varargin{:});  
  elseif isfield(varargin{1}, 'trial')
    [dat, cfg] = get_source_trial(cfg, varargin{:});
  else
    [dat, cfg] = get_source_avg(cfg, varargin{:});
  end
  cfg.dimord = 'voxel';

  % note that avgoverroi=no is implemented differently at an earlier stage
  if strcmp(cfg.avgoverroi, 'yes')
    tmp = zeros(length(roimask), size(dat,2));
    for i=1:length(roimask)
      % the data only reflects those points that are inside the brain,
      % the atlas-based mask reflects points inside and outside the brain
      roi = roimask{i}(varargin{1}.inside);
      tmp(i,:) = mean(dat(roi,:), 1);
    end
    % replace the original data with the average over each ROI
    dat = tmp;
    clear tmp roi roimask
    % remember the ROIs
    cfg.dimord = 'roi';
  end

  % check whether the original input data contains a dim, which would allow
  % for 3D reshaping and clustering with bwlabeln
  if isfield(varargin{1}, 'transform') || ((isfield(varargin{1}, 'dim') && prod(varargin{1}.dim)==size(varargin{1}.pos,1)))
    cfg.connectivity = 'bwlabeln';
  end
  
elseif isfreq || istimelock
  % get the ERF/TFR data by means of PREPARE_TIMEFREQ_DATA
  cfg.datarepresentation = 'concatenated';
  [cfg, data] = prepare_timefreq_data(cfg, varargin{:});
  cfg = rmfield(cfg, 'datarepresentation');

  dim = size(data.biol);
  if length(dim)<3
    % seems to be singleton frequency and time dimension
    dim(3)=1;
    dim(4)=1;
  elseif length(dim)<4
    % seems to be singleton time dimension
    dim(4)=1;
  end
  cfg.dimord = 'chan_freq_time';

  % the dimension of the original data (excluding the replication dimension) has to be known for clustering
  cfg.dim = dim(2:end);
  % all dimensions have to be concatenated except the replication dimension and the data has to be transposed
  dat = transpose(reshape(data.biol, dim(1), prod(dim(2:end))));

  % remove to save some memory
  data.biol = [];

  % verify that neighbours are present in case needed (not needed when
  % averaging over channels)
  if ~(isfield(cfg, 'avgoverchan') && istrue(cfg.avgoverchan)) &&...
      ~isfield(cfg,'neighbours') && isfield(cfg, 'correctm') && strcmp(cfg.correctm, 'cluster')
    ft_error('You need to specify a neighbourstructure');
    %cfg.neighbours = ft_neighbourselection(cfg,varargin{1});
  end

end

% get the design from the information in cfg and data.
if ~isfield(cfg,'design')
  ft_warning('Please think about how you would create cfg.design.  Soon the call to prepare_design will be deprecated')
  cfg.design = data.design;
  [cfg] = prepare_design(cfg);
end

if size(cfg.design,2)~=size(dat,2)
  cfg.design = transpose(cfg.design);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the statistic, using the data-independent statistical subfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine the function handle to the intermediate-level statistics function
if exist(['ft_statistics_' cfg.method], 'file')
  statmethod = str2func(['ft_statistics_' cfg.method]);
else
  ft_error('could not find the corresponding function for cfg.method="%s"\n', cfg.method);
end
fprintf('using "%s" for the statistical testing\n', func2str(statmethod));

% check that the design completely describes the data
if size(dat,2) ~= size(cfg.design,2)
  ft_error('the size of the design matrix (%d) does not match the number of observations in the data (%d)', size(cfg.design,2), size(dat,2));
end

% determine the number of output arguments
try
  % the nargout function in Matlab 6.5 and older does not work on function handles
  num = nargout(statmethod);
catch
  num = 1;
end

design = cfg.design;
cfg    = rmfield(cfg,'design'); % to not confuse lower level functions with both cfg.design and design input

% perform the statistical test 
if strcmp(func2str(statmethod),'ft_statistics_montecarlo') % because ft_statistics_montecarlo (or to be precise, clusterstat) requires to know whether it is getting source data, 
                                                        % the following (ugly) work around is necessary                                             
  if num>1
    [stat, cfg] = statmethod(cfg, dat, design);
  else
    [stat] = statmethod(cfg, dat, design);
  end
else
  if num>1
    [stat, cfg] = statmethod(cfg, dat, design);
  else
    [stat] = statmethod(cfg, dat, design);
  end
end

if isstruct(stat)
  % the statistical output contains multiple elements, e.g. F-value, beta-weights and probability
  statfield = fieldnames(stat);
else
  % only the probability was returned as a single matrix, reformat into a structure
  dum = stat; stat = []; % this prevents a Matlab warning that appears from release 7.0.4 onwards
  stat.prob = dum;
  statfield = fieldnames(stat);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add descriptive information to the output and rehape into the input format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if issource
   if ~isfield(varargin{1},'dim') % FIX ME; this is added temporarily (20110427) to cope with ft_sourceanalysis output not having a dim field since r3273
        varargin{1}.dim = [Nvoxel 1];
   end
  if isempty(cfg.roi) || strcmp(cfg.avgoverroi, 'no')
    % remember the definition of the volume, assume that they are identical for all input arguments
    try stat.dim       = varargin{1}.dim;        end
    try stat.xgrid     = varargin{1}.xgrid;      end
    try stat.ygrid     = varargin{1}.ygrid;      end
    try stat.zgrid     = varargin{1}.zgrid;      end
    try stat.inside    = varargin{1}.inside;     end
    try stat.outside   = varargin{1}.outside;    end
    try stat.pos       = varargin{1}.pos;        end
    try stat.transform = varargin{1}.transform;  end
    try stat.freq      = varargin{1}.freq;       end
    try stat.time      = varargin{1}.time;       end
  else
    stat.inside  = 1:length(roilabel);
    stat.outside = [];
    stat.label   = roilabel(:);
  end
  for i=1:length(statfield)
    tmp = getsubfield(stat, statfield{i});
    if isfield(varargin{1}, 'inside') && numel(tmp)==length(varargin{1}.inside)
      % the statistic was only computed on voxels that are inside the brain
      % sort the inside and outside voxels back into their original place
      if islogical(tmp)
        tmp(varargin{1}.inside)  = tmp;
        tmp(varargin{1}.outside) = false;
      else
        tmp(varargin{1}.inside)  = tmp;
        tmp(varargin{1}.outside) = nan;
      end
      extradim = 1;
    elseif isfield(varargin{1}, 'inside') && isfield(varargin{1}, 'freq') && isfield(varargin{1}, 'time') && numel(tmp)==length(varargin{1}.inside)*length(varargin{1}.freq)*length(varargin{1}.time)
      % the statistic is a higher dimensional matrix (here as a function of freq) computed only on the
      % inside voxels
      newtmp = zeros(length(varargin{1}.inside)+length(varargin{1}.outside), length(varargin{1}.freq), length(varargin{1}.time));
      if islogical(tmp)
        newtmp(varargin{1}.inside, :, :)  = reshape(tmp, length(varargin{1}.inside), length(varargin{1}.freq), []);
        newtmp(varargin{1}.outside, :, :) = false;
      else
        newtmp(varargin{1}.inside, :, :)  = reshape(tmp, length(varargin{1}.inside), length(varargin{1}.freq), []);
        newtmp(varargin{1}.outside, :, :) = nan;
      end
      tmp = newtmp; clear newtmp;
      extradim = length(varargin{1}.freq);
    elseif isfield(varargin{1}, 'inside') && isfield(varargin{1}, 'freq') && numel(tmp)==length(varargin{1}.inside)*length(varargin{1}.freq)
      % the statistic is a higher dimensional matrix (here as a function of freq) computed only on the
      % inside voxels
      newtmp = zeros(length(varargin{1}.inside)+length(varargin{1}.outside), length(varargin{1}.freq));
      if islogical(tmp)
        newtmp(varargin{1}.inside, :)  = reshape(tmp, length(varargin{1}.inside), []);
        newtmp(varargin{1}.outside, :) = false;
      else
        newtmp(varargin{1}.inside, :)  = reshape(tmp, length(varargin{1}.inside), []);
        newtmp(varargin{1}.outside, :) = nan;
      end
      tmp = newtmp; clear newtmp;
      extradim = length(varargin{1}.freq);
    elseif isfield(varargin{1}, 'inside') && isfield(varargin{1}, 'time') && numel(tmp)==length(varargin{1}.inside)*length(varargin{1}.time)
      % the statistic is a higher dimensional matrix (here as a function of time) computed only on the
      % inside voxels
      newtmp = zeros(length(varargin{1}.inside)+length(varargin{1}.outside), length(varargin{1}.time));
      if islogical(tmp)
        newtmp(varargin{1}.inside, :)  = reshape(tmp, length(varargin{1}.inside), []);
        newtmp(varargin{1}.outside, :) = false;
      else
        newtmp(varargin{1}.inside, :)  = reshape(tmp, length(varargin{1}.inside), []);
        newtmp(varargin{1}.outside, :) = nan;
      end
      tmp = newtmp; clear newtmp;
      extradim = length(varargin{1}.time);
    else
      extradim = 1;
    end
    if numel(tmp)==prod(varargin{1}.dim)
      % reshape the statistical volumes into the original format
      stat = setsubfield(stat, statfield{i}, reshape(tmp, varargin{1}.dim));
    else
      stat = setsubfield(stat, statfield{i}, tmp);
    end
  end
else
  haschan    = isfield(data, 'label');    % this one remains relevant, even after averaging over channels
  haschancmb = isfield(data, 'labelcmb'); % this one remains relevant, even after averaging over channels
  hasfreq = strcmp(cfg.avgoverfreq, 'no') && ~any(isnan(data.freq));
  hastime = strcmp(cfg.avgovertime, 'no') && ~any(isnan(data.time));
  stat.dimord = '';

  if haschan
    stat.dimord = [stat.dimord 'chan_'];
    stat.label  = data.label;
    chandim = dim(2);
  elseif haschancmb
    stat.dimord   = [stat.dimord 'chancmb_'];
    stat.labelcmb = data.labelcmb;
    chandim = dim(2);
  end

  if hasfreq
    stat.dimord = [stat.dimord 'freq_'];
    stat.freq   = data.freq;
    freqdim = dim(3);
  end

  if hastime
    stat.dimord = [stat.dimord 'time_'];
    stat.time   = data.time;
    timedim = dim(4);
  end

  if ~isempty(stat.dimord)
    % remove the last '_'
    stat.dimord = stat.dimord(1:(end-1));
  end

  for i=1:length(statfield)
    try
      % reshape the fields that have the same dimension as the input data
      if     strcmp(stat.dimord, 'chan')
        stat = setfield(stat, statfield{i}, reshape(getfield(stat, statfield{i}), [chandim 1]));
      elseif strcmp(stat.dimord, 'chan_time')
        stat = setfield(stat, statfield{i}, reshape(getfield(stat, statfield{i}), [chandim timedim]));
      elseif strcmp(stat.dimord, 'chan_freq')
        stat = setfield(stat, statfield{i}, reshape(getfield(stat, statfield{i}), [chandim freqdim]));
      elseif strcmp(stat.dimord, 'chan_freq_time')
        stat = setfield(stat, statfield{i}, reshape(getfield(stat, statfield{i}), [chandim freqdim timedim]));
      elseif strcmp(stat.dimord, 'chancmb_time')
        stat = setfield(stat, statfield{i}, reshape(getfield(stat, statfield{i}), [chandim timedim]));
      elseif strcmp(stat.dimord, 'chancmb_freq')
        stat = setfield(stat, statfield{i}, reshape(getfield(stat, statfield{i}), [chandim freqdim]));
      elseif strcmp(stat.dimord, 'chancmb_freq_time')
        stat = setfield(stat, statfield{i}, reshape(getfield(stat, statfield{i}), [chandim freqdim timedim]));
      elseif strcmp(stat.dimord, 'time')
        stat = setfield(stat, statfield{i}, reshape(getfield(stat, statfield{i}), [1 timedim]));
      elseif strcmp(stat.dimord, 'freq')
        stat = setfield(stat, statfield{i}, reshape(getfield(stat, statfield{i}), [1 freqdim]));
      elseif strcmp(stat.dimord, 'freq_time')
        stat = setfield(stat, statfield{i}, reshape(getfield(stat, statfield{i}), [freqdim timedim]));
      end
    end
  end
end

return % statistics_wrapper main()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for extracting the data of interest
% data resembles PCC beamed source reconstruction, multiple trials are coded in mom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dat, cfg] = get_source_pcc_mom(cfg, varargin)
Nsource = length(varargin);
Nvoxel  = length(varargin{1}.inside) + length(varargin{1}.outside);
Ninside = length(varargin{1}.inside);
dim     = varargin{1}.dim;
for i=1:Nsource
  ntrltap = sum(varargin{i}.cumtapcnt);
  dat{i}  = zeros(Ninside, ntrltap);
  for j=1:Ninside
    k = varargin{1}.inside(j);
    if ischar(varargin{i}.avg.csdlabel{1})
      % this represents the 'old-style' csdlabel, one list for all dipoles
      dipsel  = find(strcmp(varargin{i}.avg.csdlabel, 'scandip'));
    else
      % this represents the 'new-style' csdlabel, one list for each of the dipoles
      dipsel  = find(strcmp(varargin{i}.avg.csdlabel{k}, 'scandip'));
    end
    dat{i}(j,:) = reshape(varargin{i}.avg.mom{k}(dipsel,:), 1, ntrltap);
  end
end
% concatenate the data matrices of the individual input arguments
dat = cat(2, dat{:});
% remember the dimension of the source data
cfg.dim = dim;
% remember which voxels are inside the brain
cfg.inside = varargin{1}.inside;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for extracting the data of interest
% data resemples LCMV beamed source reconstruction, mom contains timecourse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dat, cfg] = get_source_lcmv_mom(cfg, varargin)
Nsource = length(varargin);
Nvoxel  = length(varargin{1}.inside) + length(varargin{1}.outside);
Ntime   = length(varargin{1}.avg.mom{varargin{1}.inside(1)});
Ninside = length(varargin{1}.inside);
dim     = [varargin{1}.dim Ntime];
dat     = zeros(Ninside*Ntime, Nsource);
for i=1:Nsource
  % collect the 4D data of this input argument
  tmp = nan(Ninside, Ntime);
  for j=1:Ninside
    k = varargin{1}.inside(j);
    tmp(j,:) = reshape(varargin{i}.avg.mom{k}, 1, dim(4));
  end
  dat(:,i) = tmp(:);
end
% remember the dimension of the source data
cfg.dim = dim;
% remember which voxels are inside the brain
cfg.inside = varargin{1}.inside;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for extracting the data of interest
% data contains single-trial or single-subject source reconstructions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dat, cfg] = get_source_trial(cfg, varargin)
Nsource = length(varargin);
Nvoxel  = length(varargin{1}.inside) + length(varargin{1}.outside);

for i=1:Nsource
  Ntrial(i) = length(varargin{i}.trial);
end
k = 1;
for i=1:Nsource
  for j=1:Ntrial(i)
    tmp = getsubfield(varargin{i}.trial(j), cfg.parameter);
    if ~iscell(tmp),
      %dim = size(tmp);
      dim = [Nvoxel 1];
    else
      dim = [Nvoxel size(tmp{varargin{i}.inside(1)})];
    end
    if i==1 && j==1 && numel(tmp)~=Nvoxel,
      ft_warning('the input-data contains more entries than the number of voxels in the volume, the data will be concatenated');
      dat    = zeros(prod(dim), sum(Ntrial)); %FIXME this is old code should be removed
    elseif i==1 && j==1 && iscell(tmp),
      ft_warning('the input-data contains more entries than the number of voxels in the volume, the data will be concatenated');
      dat    = zeros(Nvoxel*numel(tmp{varargin{i}.inside(1)}), sum(Ntrial));
    elseif i==1 && j==1,
      dat = zeros(Nvoxel, sum(Ntrial));
    end
    if ~iscell(tmp),
      dat(:,k) = tmp(:);
    else
      Ninside   = length(varargin{i}.inside);
      %tmpvec    = (varargin{i}.inside-1)*prod(size(tmp{varargin{i}.inside(1)}));
      tmpvec    = varargin{i}.inside;
      insidevec = [];
      for m = 1:numel(tmp{varargin{i}.inside(1)})
        insidevec = [insidevec; tmpvec(:)+(m-1)*Nvoxel];
      end
      insidevec = insidevec(:)';
      tmpdat  = reshape(permute(cat(3,tmp{varargin{1}.inside}), [3 1 2]), [Ninside numel(tmp{varargin{1}.inside(1)})]);
      dat(insidevec, k) = tmpdat(:);
    end
    % add original dimensionality of the data to the configuration, is required for clustering
    %FIXME: this was obviously wrong, because often trial-data is one-dimensional, so no dimensional information is present
    %cfg.dim = dim(2:end);
    k = k+1;
  end
end
if isfield(varargin{1}, 'inside')
  fprintf('only selecting voxels inside the brain for statistics (%.1f%%)\n', 100*length(varargin{1}.inside)/prod(dim));
  for j=prod(dim(2:end)):-1:1
    dat((j-1).*dim(1) + varargin{1}.outside, :) = [];
  end
end
% remember the dimension of the source data
if isfield(varargin{1}, 'dim')
  cfg.dim = varargin{1}.dim;
else
  cfg.dim = size(varargin{1}.trial(1).(cfg.parameter));
  cfg.dim(1) = numel(varargin{1}.inside);
  cfg.insideorig = varargin{1}.inside;
  cfg.inside     = varargin{1}.inside;
end

%if ~isfield(cfg, 'dim')
  %warning('for clustering on trial-based data you explicitly have to specify cfg.dim');
%end
% remember which voxels are inside the brain
cfg.inside = varargin{1}.inside;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for extracting the data of interest
% get the average source reconstructions, the repetitions are in multiple input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dat, cfg] = get_source_avg(cfg, varargin)
Nsource = length(varargin);
Nvoxel  = length(varargin{1}.inside) + length(varargin{1}.outside);
if isfield(varargin{1}, 'dim')
  dim    = varargin{1}.dim;
  inside = false(prod(dim),1);
else
  dim    = numel(varargin{1}.inside);
  inside = false(size(varargin{1}.pos,1),1);
end
inside(varargin{1}.inside) = true;

tmp     = getsubfield(varargin{1}, cfg.parameter);
if size(tmp,2)>1
  if isfield(varargin{1}, 'freq')
    Nvoxel = Nvoxel*numel(varargin{1}.freq);
    dim    = [dim numel(varargin{1}.freq)];
    inside = repmat(inside, [1 numel(varargin{1}.freq)]);
  end
  if isfield(varargin{1}, 'time')
    Nvoxel = Nvoxel*numel(varargin{1}.time);
    dim    = [dim numel(varargin{1}.time)];
    if isfield(varargin{1},'freq')
      inside = repmat(inside, [1 1 numel(varargin{1}.time)]);
    else
      inside = repmat(inside, [1 numel(varargin{1}.time)]);
    end
  end
end
dat = zeros(Nvoxel, Nsource);
for i=1:Nsource
  tmp = getsubfield(varargin{i}, cfg.parameter);
  dat(:,i) = tmp(:);
end
insideorig = find(inside(:,1));
inside = find(inside(:));
if isfield(varargin{1}, 'inside')
  if isfield(varargin{1}, 'pos')
    npos = size(varargin{1}.pos,1);
  elseif isfield(varargin{1}, 'dim')
    npos = prod(varargin{1}.dim);
  else
    npos = length(varargin{1}.inside); % uninformative, at least prevents crash
  end
  fprintf('only selecting voxels inside the brain for statistics (%.1f%%)\n', 100*length(varargin{1}.inside)/npos);
  dat = dat(inside,:);
end
% remember the dimension of the source data
cfg.dim = dim;
% remember which voxels are inside the brain
cfg.inside     = inside(:);
cfg.insideorig = insideorig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for creating a design matrix
% should be called in the code above, or in prepare_design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cfg] = get_source_design_pcc_mom(cfg, varargin)
% should be implemented

function [cfg] = get_source_design_lcmv_mom(cfg, varargin)
% should be implemented

function [cfg] = get_source_design_trial(cfg, varargin)
% should be implemented

function [cfg] = get_source_design_avg(cfg, varargin)
% should be implemented
