function [stat, cfg] = ft_statistics_crossvalidate(cfg, dat, design)

% FT_STATISTICS_CROSSVALIDATE performs cross-validation using a prespecified
% multivariate analysis given by cfg.mva
%
% Use as
%   stat = ft_timelockstatistics(cfg, data1, data2, data3, ...)
%   stat = ft_freqstatistics    (cfg, data1, data2, data3, ...)
%   stat = ft_sourcestatistics  (cfg, data1, data2, data3, ...)
%
% Options:
%   cfg.mva           = a multivariate analysis (default = {dml.standardizer dml.svm})
%   cfg.statistic     = a cell-array of statistics to report (default = {'accuracy' 'binomial'})
%   cfg.nfolds        = number of cross-validation folds (default = 5)
%   cfg.resample      = true/false; upsample less occurring classes during
%                       training and downsample often occurring classes
%                       during testing (default = false)
%
% Returns:
%   stat.statistic    = the statistics to report
%   stat.model        = the models associated with this multivariate analysis
%
% See also FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS, FT_SOURCESTATISTICS

% Copyright (c) 2007-2011, F.C. Donders Centre, Marcel van Gerven
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

% do a sanity check on the input data
assert(isnumeric(dat),    'this function requires numeric data as input, you probably want to use FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS instead');
assert(isnumeric(design), 'this function requires numeric data as input, you probably want to use FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS instead');

cfg.mva       = ft_getopt(cfg, 'mva');
cfg.statistic = ft_getopt(cfg, 'statistic', {'accuracy', 'binomial'});
cfg.nfolds    = ft_getopt(cfg, 'nfolds',   5);
cfg.resample  = ft_getopt(cfg, 'resample', false);
cfg.cv        = ft_getopt(cfg, 'cv', []);
cfg.cv.type   = ft_getopt(cfg.cv, 'type', 'nfold');


% specify classification procedure or ensure it's the correct object
if isempty(cfg.mva)
  cfg.mva = dml.analysis({ dml.standardizer('verbose',true) ...
                           dml.svm('verbose',true)});
elseif ~isa(cfg.mva,'dml.analysis')
  cfg.mva = dml.analysis(cfg.mva);
end

cv_options = {'mva', cfg.mva, 'type', cfg.cv.type, 'resample', cfg.resample, 'compact', true, 'verbose', true};
if strcmp(cfg.cv.type, 'nfold')
  cv_options = cat(2, cv_options, {'folds', cfg.nfolds});
end
cv = dml.crossvalidator(cv_options{:});

if any(isinf(dat(:)))
  ft_warning('Inf encountered; replacing by zeros');
  dat(isinf(dat(:))) = 0;
end

if any(isnan(dat(:)))
  ft_warning('Nan encountered; replacing by zeros');
  dat(isnan(dat(:))) = 0;
end

% perform everything
cv = cv.train(dat',design');

% extract the statistic of interest
s = cv.statistic(cfg.statistic);
for i=1:length(cfg.statistic)
 stat.statistic.(cfg.statistic{i}) = s{i};
end

% get the model averaged over folds
stat.model = cv.model;

fn = fieldnames(stat.model{1});
if any(ismember(fn,  {'weights', 'primal'})),
  selfn = find(ismember(fn, {'weights', 'primal'}));
  
  % the mean subtraction is needed only once, but speeds up the covariance
  % computation
  dat = bsxfun(@minus, dat, nanmean(dat,2)); 
  dat_transp = dat.';
  for j=1:numel(selfn)
    % create the 'encoding' matrix from the weights, as per Haufe 2014.
    %covdat = cov(dat');
    for i=1:length(stat.model)
      i
      W = stat.model{i}.(fn{selfn});
      
      sW   = size(W);
      sdat = size(dat);
      if sW(2)==sdat(1) && sW(1)~=sdat(1)
        W = transpose(W);
      end
      
      M    = dat'*W;
      covM = cov(M);
      WcovM = (W/covM)./(size(dat,2)-1); % with the correction term for the covariance computation
      
      %stat.model{i}.(sprintf('%sinv',fn{selfn})) = covdat*W/covM;
      stat.model{i}.(sprintf('%sinv',fn{selfn})) = dat*(dat_transp*WcovM);
      
    end
  end
end
fn = fieldnames(stat.model{1}); % update the fieldnames, because some might have been added

fn = fieldnames(stat.model{1}); % may now also contain weightsinv
for i=1:length(stat.model)
  for k=1:length(fn)
    if numel(stat.model{i}.(fn{k}))==prod(cfg.dim)
      stat.model{i}.(fn{k}) = reshape(stat.model{i}.(fn{k}),cfg.dim);
    end
  end

end

% required
stat.trial = [];

% add some stuff to the cfg
cfg.cv = cv;
