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
%                    'stats'         (soon deprecated) uses a parametric test from the MATLAB statistics toolbox,
%
% The other cfg options depend on the method that you select. You
% should read the help of the respective subfunction FT_STATISTICS_XXX
% for the corresponding configuration options and for a detailed
% explanation of each method.
%
% See also FT_SOURCEANALYSIS, FT_SOURCEDESCRIPTIVES, FT_SOURCEGRANDAVERAGE

% Deprecated cfg.method options:
%                    'parametric'    uses the MATLAB statistics toolbox (very similar to 'stats'),
%                    'randomization' uses randomization of the data prior to source reconstruction,
%                    'randcluster'   uses randomization of the data prior to source reconstruction
%                                    in combination with spatial clusters.

% FIXME the following needs to be reimplemented
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

% Copyright (C) 2005-2014, Robert Oostenveld
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


% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'required',    {'method', 'design'});
cfg = ft_checkconfig(cfg, 'renamed',     {'approach',   'method'});


%%%%%%%%%

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'required',    {'method', 'design'});
cfg = ft_checkconfig(cfg, 'renamed',     {'approach',   'method'});
cfg = ft_checkconfig(cfg, 'forbidden',   {'transform'});
cfg = ft_checkconfig(cfg, 'forbidden',   {'trials'}); % this used to be present until 24 Dec 2014, but was deemed too confusing by Robert

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'source', 'feedback', 'no');
  varargin{i} = ft_datatype_source(varargin{i}, 'version', 'upcoming');
end

% set the defaults
cfg.parameter   = ft_getopt(cfg, 'parameter'); % default is set below
cfg.correctm    = ft_getopt(cfg, 'correctm');
cfg.latency     = ft_getopt(cfg, 'latency',     'all');
cfg.avgovertime = ft_getopt(cfg, 'avgovertime', 'no');
cfg.frequency   = ft_getopt(cfg, 'frequency',   'all');
cfg.avgoverfreq = ft_getopt(cfg, 'avgoverfreq', 'no');
cfg.parameter   = ft_getopt(cfg, 'parameter', 'pow');

if strncmp(cfg.parameter, 'avg.', 4)
  cfg.parameter = cfg.parameter(5:end); % remove the 'avg.' part
end
for i=1:length(varargin)
  assert(isfield(varargin{i}, cfg.parameter), 'data does not contain parameter "%s"', cfg.parameter);
end

% ensure that the data in all inputs has the same channels, time-axis, etc.
tmpcfg = keepfields(cfg, {'frequency', 'avgoverfreq', 'latency', 'avgovertime', 'showcallinfo'});
[varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
% restore the provenance information
[cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});

dimord = getdimord(varargin{1}, cfg.parameter);
dimtok = tokenize(dimord, '_');
dimsiz = getdimsiz(varargin{1}, cfg.parameter);
dimsiz(end+1:length(dimtok)) = 1; % there can be additional trailing singleton dimensions
rptdim = find( strcmp(dimtok, 'subj') |  strcmp(dimtok, 'rpt') |  strcmp(dimtok, 'rpttap'));
datdim = find(~strcmp(dimtok, 'subj') & ~strcmp(dimtok, 'rpt') & ~strcmp(dimtok, 'rpttap'));
datsiz = dimsiz(datdim);

% pass these fields to the low-level functions, they should be removed further down
cfg.dimord = sprintf('%s_', dimtok{datdim});
cfg.dimord = cfg.dimord(1:end-1); % remove trailing _
cfg.dim    = dimsiz(datdim);

if isfield(varargin{1}, 'dim')
  % the positions can be mapped onto a 3-D volume
  if prod(cfg.dim)==prod(varargin{1}.dim)
    % one value per grid position, i.e. [nx ny nz]
    cfg.origdim = varargin{1}.dim;
  else
    % multiple value per grid position
    cfg.origdim = [varargin{1}.dim cfg.dim(2:end)];
  end
end

if isfield(varargin{1}, 'inside'), 
  cfg.inside = varargin{1}.inside;
else
  cfg.inside = true(size(varargin{1}.pos,1),1);
end

% also create an inside vector for the reshaped data, which is needed
% when there are multiple values per grid position
cfg.originside = cfg.inside;
cfg.inside     = repmat(cfg.inside, prod(cfg.dim)./numel(cfg.inside), 1);


if numel(cfg.dim)==1
  cfg.dim(2) = 1;  % add a trailing singleton dimensions
end

if isempty(rptdim)
  % repetitions are across multiple inputs
  dat = nan(prod(dimsiz), length(varargin));
  for i=1:length(varargin)
    tmp = varargin{i}.(cfg.parameter);
    dat(:,i) = tmp(:);
  end
else
  % repetitions are within inputs
  dat = cell(size(varargin));
  for i=1:length(varargin)
    tmp = varargin{i}.(cfg.parameter);
    if rptdim~=1
      % move the repetitions to the first dimension
      tmp = permute(tmp, [rptdim datdim]);
    end
    dat{i} = reshape(tmp, size(tmp,1), []);
  end
  dat = cat(1, dat{:});   % repetitions along 1st dimension
  dat = dat';             % repetitions along 2nd dimension
end
dat = dat(cfg.inside,:);

if size(cfg.design,2)~=size(dat,2)
  cfg.design = transpose(cfg.design);
end

design = cfg.design;

% determine the function handle to the intermediate-level statistics function
if exist(['ft_statistics_' cfg.method], 'file')
  statmethod = str2func(['ft_statistics_' cfg.method]);
else
  error('could not find the corresponding function for cfg.method="%s"\n', cfg.method);
end
fprintf('using "%s" for the statistical testing\n', func2str(statmethod));

% check that the design completely describes the data
if size(dat,2) ~= size(cfg.design,2)
  error('the length of the design matrix (%d) does not match the number of observations in the data (%d)', size(cfg.design,2), size(dat,2));
end

% determine the number of output arguments
try
  % the nargout function in MATLAB 6.5 and older does not work on function handles
  num = nargout(statmethod);
catch
  num = 1;
end

% perform the statistical test
if num>1
  [stat, cfg] = statmethod(cfg, dat, design);
  cfg         = rollback_provenance(cfg); % ensure that changes to the cfg are passed back to the right level
else
  [stat] = statmethod(cfg, dat, design);
end

if ~isstruct(stat)
  % only the probability was returned as a single matrix, reformat into a structure
  stat = struct('prob', stat);
end

% the statistical output contains multiple elements, e.g. F-value, beta-weights and probability
fn = fieldnames(stat);

for i=1:length(fn)
  if numel(stat.(fn{i}))==sum(cfg.inside)
    % get the data back onto the inside positions
    tmp = nan+zeros(numel(cfg.inside),1);
    tmp(cfg.inside) = stat.(fn{i});
    stat.(fn{i})    = tmp;
  end  
  if numel(stat.(fn{i}))==prod(datsiz)
    % reformat into the same dimensions as the input data
    stat.(fn{i}) = reshape(stat.(fn{i}), [datsiz 1]);
  end
end

% describe the dimensions of the output data
stat.dimord = cfg.dimord;

% copy the descripive fields into the output
stat = copyfields(varargin{1}, stat, {'freq', 'time', 'pos', 'dim', 'transform', 'tri', 'inside'});

% these were only present to inform the low-level functions
cfg = removefields(cfg, {'dim', 'dimord', 'tri', 'inside'});

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   varargin
ft_postamble provenance stat
ft_postamble history    stat
ft_postamble savevar    stat
