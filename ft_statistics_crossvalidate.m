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

% Copyright (c) 2007-2011, Marcel van Gerven, F.C. Donders Centre
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

cfg.mva       = ft_getopt(cfg, 'mva');
cfg.statistic = ft_getopt(cfg, 'statistic', {'accuracy', 'binomial'});
cfg.nfolds    = ft_getopt(cfg, 'nfolds',   5);
cfg.resample  = ft_getopt(cfg, 'resample', false);

% specify classification procedure or ensure it's the correct object
if isempty(cfg.mva)
  cfg.mva = dml.analysis({ dml.standardizer('verbose',true) ...
                           dml.svm('verbose',true)});
elseif ~isa(cfg.mva,'dml.analysis')
  cfg.mva = dml.analysis(cfg.mva);
end

cv = dml.crossvalidator('mva', cfg.mva, 'type', 'nfold', 'folds', cfg.nfolds,...
  'resample', cfg.resample, 'compact', true, 'verbose', true);

if any(isinf(dat(:)))
  warning('Inf encountered; replacing by zeros');
  dat(isinf(dat(:))) = 0;
end

if any(isnan(dat(:)))
  warning('Nan encountered; replacing by zeros');
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
if any(strcmp(fn, 'weights')),
  % create the 'encoding' matrix from the weights, as per Haufe 2014.
  covdat = cov(dat');
  for i=1:length(stat.model)
    W = stat.model{i}.weights;
    M = dat'*W;
    covM = cov(M);
    stat.model{i}.weightsinv = covdat*W*inv(covM);
  end
end

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
