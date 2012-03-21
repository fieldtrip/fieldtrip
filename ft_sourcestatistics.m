function [stat] = ft_sourcestatistics(cfg, varargin)

% FT_SOURCESTATISTICS computes the probability for a given null-hypothesis using
% a parametric statistical test or using a non-parametric randomization test.
% 
% Use as
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% where the input data is the result from FT_SOURCEANALYSIS, FT_SOURCEDESCRIPTIVES
% or FT_SOURCEGRANDAVERAGE.  The source structures should be spatially alligned
% to each other and should have the same positions for the source grid.
%
% The configuration should contain the following option for data selection
%   cfg.parameter  = string, describing the functional data to be processed, e.g. 'pow', 'nai' or 'coh'
%
% Furthermore, the configuration should contain:
%   cfg.method       = different methods for calculating the probability of the null-hypothesis,
%                    'montecarlo'    uses a non-parametric randomization test to get a Monte-Carlo estimate of the probability,
%                    'analytic'      uses a parametric test that results in analytic probability,
%                    'stats'         (soon deprecated) uses a parametric test from the Matlab statistics toolbox,
%                    'parametric'    uses the Matlab statistics toolbox (very similar to 'stats'),
%                    'randomization' uses randomization of the data prior to source reconstruction,
%                    'randcluster'   uses randomization of the data prior to source reconstruction 
%                                    in combination with spatial clusters.
%
% You can restrict the statistical analysis to regions of interest (ROIs)
% or to the average value inside ROIs using the following options:
%   cfg.atlas        = filename of the atlas
%   cfg.roi          = string or cell of strings, region(s) of interest from anatomical atlas
%   cfg.avgoverroi   = 'yes' or 'no' (default = 'no')
%   cfg.hemisphere   = 'left', 'right', 'both', 'combined', specifying this is
%                      required when averaging over regions
%   cfg.inputcoord   = 'mni' or 'tal', the coordinate system in which your source 
%                      reconstruction is expressed
%
% The other cfg options depend on the method that you select. You
% should read the help of the respective subfunction STATISTICS_XXX
% for the corresponding configuration options and for a detailed
% explanation of each method.
%
% See also FT_SOURCEANALYSIS, FT_SOURCEDESCRIPTIVES, FT_SOURCEGRANDAVERAGE

% Copyright (C) 2005-2008, Robert Oostenveld
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig
ft_preamble loadvar varargin

% this wrapper should be compatible with the already existing statistical
% functions that only work for source input data

if ~isfield(cfg, 'implementation'), cfg.implementation = 'old'; end

cfg = ft_checkconfig(cfg, 'forbidden',   {'trials'});

if strcmp(cfg.implementation, 'old'),
  
  %--------------------------------
  % use the original implementation

  % check if the input data is valid for this function
  for i=1:length(varargin)
    if isfield(cfg, 'roi') && ~isempty(cfg.roi)
      varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'source', 'feedback', 'no', 'inside', 'index');
    else
      varargin{i} = ft_checkdata(varargin{i}, 'datatype', {'source', 'volume'}, 'feedback', 'no', 'inside', 'index');
    end
  end
  
  if strcmp(cfg.method, 'parametric')
    % use the source-specific statistical subfunction
    stat = sourcestatistics_parametric(cfg, varargin{:});
  elseif strcmp(cfg.method, 'randomization')
    % use the source-specific statistical subfunction
    stat = sourcestatistics_randomization(cfg, varargin{:});
  elseif strcmp(cfg.method, 'randcluster')
    % use the source-specific statistical subfunction
    stat = sourcestatistics_randcluster(cfg, varargin{:});
  else
    [stat, cfg] = statistics_wrapper(cfg, varargin{:});
  end

  
elseif strcmp(cfg.implementation, 'new')
  
  %---------------------------
  % use the new implementation
  issource = ft_datatype(varargin{1}, 'source');
  isvolume = ft_datatype(varargin{1}, 'volume'); 
 
  % check if the input data is valid for this function
  for i=1:length(varargin)
    if isfield(cfg, 'roi') && ~isempty(cfg.roi)
      % FIXME implement roi-based statistics for the new implementation
      % (code is copied over from the old implementation but not yet tested
      error('roi based sourcestatistics is not yet implemented for the new implementation');
      varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'source', 'feedback', 'no', 'inside', 'index');
    else
      varargin{i} = ft_checkdata(varargin{i}, 'datatype', {'source', 'volume'}, 'feedback', 'no', 'inside', 'index', 'sourcerepresentation', 'new');
      if strcmp(cfg.parameter, 'pow') && ~isfield(varargin{i}, 'pow'),
        varargin{i} = ft_checkdata(varargin{i}, 'sourcerepresentation', 'new', 'haspow', 'yes');
      end
    end
  end
  
  if any(strcmp(cfg.method, {'parametric' 'randomization' 'randcluster'}))
    % FIXME only supported for old-style source representation
    for i = 1:numel(varargin)
      varargin{i} = ft_checkdata(varargin{i}, 'sourcerepresentation', 'old');
    end
    
    if exist(['statistics_',cfg.method]),
      statmethod = str2func(['statistics_' cfg.method]);
    else
      error(sprintf('could not find the corresponding function for cfg.method="%s"\n', cfg.method));
    end
    stat = statmethod(cfg, varargin{:});

  else
    
    % convert representation of input data to new style
    for i = 1:numel(varargin)
      varargin{i} = ft_checkdata(varargin{i}, 'sourcerepresentation', 'new');
    end

    % check the input configuration
    cfg = ft_checkconfig(cfg, 'renamed',     {'approach',   'method'});
    cfg = ft_checkconfig(cfg, 'required',    {'method', 'parameter'});
    cfg = ft_checkconfig(cfg, 'forbidden',   {'transform'});

    % set the defaults
    if ~isfield(cfg, 'channel'),     cfg.channel = 'all';            end
    if ~isfield(cfg, 'latency'),     cfg.latency = 'all';            end
    if ~isfield(cfg, 'frequency'),   cfg.frequency = 'all';          end
    if ~isfield(cfg, 'roi'),         cfg.roi = [];                   end
    if ~isfield(cfg, 'avgoverchan'), cfg.avgoverchan = 'no';         end
    if ~isfield(cfg, 'avgovertime'), cfg.avgovertime = 'no';         end
    if ~isfield(cfg, 'avgoverfreq'), cfg.avgoverfreq = 'no';         end
    if ~isfield(cfg, 'avgoverroi'),  cfg.avgoverroi = 'no';          end

    % test that all source inputs have the same dimensions and are spatially aligned
    for i=2:length(varargin)
      if isfield(varargin{1}, 'dim') && (numel(varargin{i}.dim)~=length(varargin{1}.dim) || ~all(varargin{i}.dim==varargin{1}.dim))
        error('dimensions of the source reconstructions do not match, use VOLUMENORMALISE first');
      end
      if isfield(varargin{1}, 'pos') && (numel(varargin{i}.pos(:))~=numel(varargin{1}.pos(:)) || ~all(varargin{i}.pos(:)==varargin{1}.pos(:)))
        error('grid locations of the source reconstructions do not match, use VOLUMENORMALISE first');
      end
    end 
  
    Nsource = length(varargin);
    Nvoxel  = length(varargin{1}.inside) + length(varargin{1}.outside);

    %FIXME ft_selectdata should be used for the subselection
    %FIXME ft_selectdata has to be adjusted to work with new style source data
    %if isfield(varargin{1}, 'freq') && ~strcmp(cfg.frequency, 'all'),
    %  for i=1:length(varargin)
    %    varargin{i} = ft_selectdata(varargin{i}, 'foilim', cfg.frequency, ...
    %                             'avgoverfreq', cfg.avgoverfreq);
    %  end
    %end

    % this part contains the functionality of the old statistics_wrapper
    % with source data in the input
    if ~isempty(cfg.roi)
      if ischar(cfg.roi)
        cfg.roi = {cfg.roi};
      end
      % the source representation should specify the position of each voxel in MNI coordinates
      x = varargin{1}.pos(:,1);  % this is from left (negative) to right (positive)
      % determine the mask to restrict the subsequent analysis
      % process each of the ROIs, and optionally also left and/or right seperately
      roimask  = {};
      roilabel = {};
      for i=1:length(cfg.roi)
        tmpcfg.roi = cfg.roi{i};
        tmpcfg.inputcoord = cfg.inputcoord;
        tmpcfg.atlas = cfg.atlas;
        tmp = ft_volumelookup(tmpcfg, varargin{1});
        if strcmp(cfg.avgoverroi, 'no') && ~isfield(cfg, 'hemisphere')
          % no reason to deal with seperated left/right hemispheres
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
          % deal seperately with the voxels on the left and right side of the brain
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
          roilabel{end+1} = cfg.roi{i};

        else
          error('incorrect specification of cfg.hemisphere');
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
    end % if ~isempty cfg.roi
    
    % get the required source level data  
    [dat, cfg] = getfunctional(cfg, varargin{:});
    
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
  end

  %get the design from the information in the cfg and data
  if ~isfield(cfg, 'design'),
    cfg.design = data.design;
    cfg        = prepare_design(cfg);
  end
  
  if size(cfg.design, 2)~=size(dat, 2)
    cfg.design = transpose(cfg.design);
  end
  
  % determine the function handle to the intermediate-level statistics function
  if exist(['statistics_' cfg.method])
    statmethod = str2func(['statistics_' cfg.method]);
  else
    error(sprintf('could not find the corresponding function for cfg.method="%s"\n', cfg.method));
  end
  fprintf('using "%s" for the statistical testing\n', func2str(statmethod));
  
  % check that the design completely describes the data
  if size(dat,2) ~= size(cfg.design,2)
    error('the size of the design matrix does not match the number of observations in the data');
  end
  
  % determine the number of output arguments
  try
    % the nargout function in Matlab 6.5 and older does not work on function handles
    num = nargout(statmethod);
  catch
    num = 1;
  end
  
  % perform the statistical test 
  if strcmp(func2str(statmethod),'statistics_montecarlo') 
    % because statistics_montecarlo (or to be precise, clusterstat)
    % requires to know whether it is getting source data, 
    % the following (ugly) work around is necessary                                             
    if num>1
      [stat, cfg] = statmethod(cfg, dat, cfg.design, 'issource', 1);
    else
      [stat] = statmethod(cfg, dat, cfg.design, 'issource', 1);
    end
  else
    if num>1
      [stat, cfg] = statmethod(cfg, dat, cfg.design);
    else
      [stat] = statmethod(cfg, dat, cfg.design);
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
  
  % add descriptive information to the output and rehape into the input format
  if isempty(cfg.roi) || strcmp(cfg.avgoverroi, 'no')
    % remember the definition of the volume, assume that they are identical for all input arguments
    try, stat.dim       = varargin{1}.dim;        end
    try, stat.xgrid     = varargin{1}.xgrid;      end
    try, stat.ygrid     = varargin{1}.ygrid;      end
    try, stat.zgrid     = varargin{1}.zgrid;      end
    try, stat.inside    = varargin{1}.inside;     end
    try, stat.outside   = varargin{1}.outside;    end
    try, stat.pos       = varargin{1}.pos;        end
    try, stat.transform = varargin{1}.transform;  end
  else
    stat.inside  = 1:length(roilabel);
    stat.outside = [];
    stat.label   = roilabel(:);
  end
    
  % additional descriptive fields
  hasfreq = strcmp(cfg.avgoverfreq, 'no') && isfield(varargin{1},'freq');
  hastime = strcmp(cfg.avgovertime, 'no') && isfield(varargin{1},'time');
  stat.dimord = 'pos_';

  if hasfreq, 
    stat.dimord = [stat.dimord, 'freq_']; 
    stat.freq   = varargin{1}.freq;
    nfreq       = numel(varargin{1}.freq);
  else
    nfreq       = 1;
  end
  if hastime,
    stat.dimord = [stat.dimord, 'time_'];
    stat.time   = varargin{1}.time;
    ntime       = numel(varargin{1}.time);
  else
    ntime       = 1;
  end
  stat.dimord = stat.dimord(1:end-1);

  if issource, 
    if hasfreq,
       newdim = [size(stat.pos,1) nfreq ntime];
    else
       newdim = [size(stat.pos,1) ntime];
    end
  elseif isvolume,
    if hasfreq,
      newdim = [stat.dim nfreq ntime];
    else
      newdim = [stat.dim ntime];
    end
  end
  
  for i=1:length(statfield)
    tmp   = getsubfield(stat, statfield{i});
    tmp2  = [];
    ntmp  = numel(tmp);
    if hasfreq, 
      tmpdim = [ntmp/(nfreq*ntime) nfreq ntime];
    else
      tmpdim = [ntmp/ntime ntime];
    end
    if isfield(varargin{1}, 'inside') && numel(tmp)==nfreq*ntime*length(varargin{1}.inside)
      % the statistic was only computed on voxels that are inside the brain
      % sort the inside and outside voxels back into their original place
      if islogical(tmp)
        if hasfreq,
          tmp2 = logical(zeros(prod(varargin{1}.dim),nfreq,ntime));
          tmp2(varargin{1}.inside,  1:nfreq, 1:ntime) = reshape(tmp, [ntmp/(nfreq*ntime) nfreq ntime]);
        else
          tmp2 = logical(zeros(prod(varargin{1}.dim),nfreq,ntime));
          tmp2(varargin{1}.inside,  1:nfreq, 1:ntime) = reshape(tmp, [ntmp/ntime ntime]);
        end
      else
        if hasfreq,
          tmp2 = zeros(prod(varargin{1}.dim),nfreq,ntime)+nan;
          tmp2(varargin{1}.inside,  1:nfreq, 1:ntime) = reshape(tmp, [ntmp/(nfreq*ntime) nfreq ntime]);
        else
          tmp2 = zeros(prod(varargin{1}.dim),nfreq,ntime)+nan;
          tmp2(varargin{1}.inside,  1:nfreq, 1:ntime) = reshape(tmp, [ntmp/ntime ntime]);
        end
      end
    end
    if numel(tmp2)==prod(newdim)
      % reshape the statistical volumes into the original format
      stat = setsubfield(stat, statfield{i}, reshape(tmp2, newdim));
    end
  end

else
  error('cfg.implementation can be only old or new');
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous varargin
ft_postamble history stat
ft_postamble savevar stat


%-----------------------------------------------------
%subfunction to extract functional data from the input
%and convert it into a 2D representation
function [dat, cfg] = getfunctional(cfg, varargin)

%FIXME think of how this generalizes to volumetric data (boolean inside etc)
Nsource = numel(varargin);
Nvox    = size(varargin{1}.pos, 1);
inside  = varargin{1}.inside;
Ninside = numel(inside);

dimord = varargin{1}.([cfg.parameter,'dimord']);
dimtok = tokenize(dimord, '_');

%check whether the requested parameter is represented as cell-array
x1 = strfind(dimord, '{');
x2 = strfind(dimord, '}');
if ~isempty(x1) && ~isempty(x2)
  cellparam = 1;
  cellsiz   = size(varargin{1}.(cfg.parameter));
  %this only explicitly keeps the first singleton dimension
  %the last to be removed
  cellsiz(find(cellsiz(2:end)==1)+1) = [];
  cellsiz(cellsiz==Nvox)             = Ninside;
else
  cellparam = 0;
end

%check whether there are single observations/subjects in the data
rptdim = ~cellfun(@isempty, strfind(dimtok, 'rpt')) | ~cellfun(@isempty, strfind(dimtok, 'subj'));
hasrpt = sum(rptdim);

if hasrpt && Nsource>1,
  error('only a single input with multiple observations or multiple inputs with a single observation are supported');
end

if hasrpt,
  if cellparam,
    tmpsiz  = [cellsiz size(varargin{1}.(cfg.parameter){inside(1)})];
    tmp     = zeros(tmpsiz);
    %FIXME what about volumetric data?
    for k = 1:Ninside
      tmp(k,:,:,:,:,:) = varargin{1}.(cfg.parameter){inside(k)};
    end
    %tmp    = cell2mat(varargin{1}.(cfg.parameter)(inside,:,:,:,:));
  else
    if find(rptdim)==1,
      tmp = varargin{1}.(cfg.parameter)(:,inside,:,:,:);
    else
      tmp = varargin{1}.(cfg.parameter)(inside,:,:,:,:);
    end
  end
  
  if numel(rptdim)==1,
    rptdim = [rptdim 0];
  end
  %put the repetition dimension to the last dimension
  tmp = permute(tmp, [find(rptdim==0) find(rptdim==1)]);

  %reshape the data to 2D
  siz = size(tmp);
  dat = reshape(tmp, [prod(siz(1:end-1)) siz(end)]);
  
else
  
  for k = 1:Nsource
    %check for cell-array representation in the input data
    %FIXME this assumes positions to be in the first dimension always
    %FIXME what about volumetric dadta
    if cellparam
      tmp = cell2mat(varargin{k}.(cfg.parameter)(inside,:,:,:,:));
    else
      tmp = varargin{k}.(cfg.parameter)(inside,:,:,:,:);
    end

    %allocate memory    
    if k==1, dat = zeros(numel(tmp), Nsource); end

    %reshape the data
    siz      = size(tmp);
    dat(:,k) = tmp(:);
  end
end
%cfg.dim     = varargin{1}.dim;
%cfg.inside  = varargin{1}.inside; %FIXME take the intersection between all inputs
%FIXME don't do the previous lines in order to take the unfolded inside
%across the dimensions in the input and to get the 4D dimensionality
%correct
cfg.dimord  = 'voxel';
cfg.origdim = [cfg.dim siz(2:end-1)];
