function stat = statistics_crossvalidate(cfg, dat, design)

% STATISTICS_CROSSVALIDATE performs cross-validation using a prespecified
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
%
% Returns:
%   stat.statistic    = the statistics to report
%   stat.pvalue       = p-value for the specified significance test
%   stat.model<n>     = the n-th model associated with this multivariate analysis
%

% Copyright (c) 2007-2011, Marcel van Gerven, F.C. Donders Centre
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

% specify classification procedure
  
if ~isfield(cfg,'mva')
  cfg.mva = dml.analysis({ ...
    dml.standardizer('verbose',true) ...
    dml.svm('verbose',true) ...
    });
else
  if ~isa(cfg.mva,'dml.analysis')
    cfg.mva = dml.analysis(cfg.mva);
  end
end

if ~isfield(cfg,'statistic'),
  cfg.statistic = {'accuracy' 'binomial'};
end

if ~isfield(cfg,'nfolds'), cfg.nfolds = 5; end

cv = dml.crossvalidator('mva',cfg.mva,'type','nfold','folds',cfg.nfolds,'compact',true,'verbose',true);
  
if isfield(cfg, 'trainfolds')
  cv.trainfolds = cfg.trainfolds;
end
if isfield(cfg, 'testfolds')
  cv.testfolds = cfg.testfolds;
end

% perform everything!
cv = cv.train(dat',design');

% the statistic of interest
s = cv.statistic(cfg.statistic);
for i=1:length(cfg.statistic)
 stat.statistic.(cfg.statistic{i}) = s{i};
end

% get the model averaged over folds
stat.model = cv.model; 

fn = fieldnames(stat.model{1});
for i=1:length(stat.model)
  
  for k=1:length(fn)
    if prod(size(stat.model{i}.(fn{k})))==prod(cfg.dim)
      stat.model{i}.(fn{k}) = squeeze(reshape(stat.model{i}.(fn{k}),cfg.dim));
    end
  end
     
end
  
% required
stat.trial = [];
