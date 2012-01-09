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
%   cfg.metric        = the performance metric to report (default = 'accuracy')
%   cfg.sigtest       = the significance test to report (default = 'mcnemar')
%   cfg.cv            = crossvalidator object
% overloads the following
%   cfg.mva           = a multivariate analysis (default = {standardizer svmmethod})
%   cfg.nfolds        = number of folds (default = 10)
%   cfg.compact       = whether or not to save the mva procedure (true)
%   cfg.model         = whether or not to save the average model (true)
%
% Returns:
%   stat.performance  = computed using the specified metric
%   stat.pvalue       = p-value for the specified significance test
%   stat.cv           = the trained crossvalidator
%   stat.model<n>     = the n-th model associated with this multivariate analysis
%
% See also CROSSVALIDATE, MVA

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
if isfield(cfg,'cv')
  cv = cfg.cv;
else
  
  if ~isfield(cfg,'mva')
    cfg.mva = ft_mv_analysis({ ...
      ft_mv_standardizer('verbose',true) ...
      ft_mv_svm('verbose',true) ...
      });
  else
    if ~isa(cfg.mva,'ft_mv_analysis')
      cfg.mva = ft_mv_analysis(cfg.mva);
    end
  end
  
  if ~isfield(cfg,'nfolds'), cfg.nfolds = 10; end
  if ~isfield(cfg,'compact'), cfg.compact = true; end
  
  cv = ft_mv_crossvalidator('mva',cfg.mva,'nfolds',cfg.nfolds,'compact',cfg.compact,'verbose',true);
  
end

if ~isfield(cfg,'metric'),
  cv.metric = 'accuracy';
else
  cv.metric = cfg.metric;
end

if ~isfield(cfg,'sigtest'),
  cv.sigtest = 'binomial';
else
  cv.sigtest = cfg.sigtest;
end

if isfield(cfg, 'trainfolds')
  cv.trainfolds = cfg.trainfolds;
end
if isfield(cfg, 'testfolds')
  cv.testfolds = cfg.testfolds;
end

% check for transfer learning; this is implemented by cfg.dataset,
% indicating the dataset number for each element in the design matrix

if isfield(cfg,'dataset') && ~isempty(cfg.dataset)
  
  % split up the datasets
  
  tmpdat = dat';
  tmpdesign = design';
  
  n = max(cfg.dataset);
  
  dat = cell(1,n);
  design = cell(1,n);
  for c=1:n
    dat{c}    = tmpdat(cfg.dataset == c,:);
    design{c} = tmpdesign(cfg.dataset == c,:);
  end
  
else
  
  dat = dat';
  design = design';
  
end

% perform everything ;o)
cv = cv.train(dat,design);

% the statistic of interest
stat.performance = cv.performance();

% null-hypothesis rejected?
stat.pvalue = cv.significance();

if ~cfg.compact
  stat.X = dat;
  stat.Y = design;
end

% get the models

m = cv.model;
desc = cv.description;

for c=1:size(m,1) % iterate over parameter types
  
  if size(m,2) > 1 % transfer learning
    
    for j=1:size(m,2)
      
      if numel(m{c,j}) == prod(cfg.dim)
        stat.(sprintf('model%d_%d',c,j)) = reshape(m{c,j},cfg.dim);
      end
    end
    
  else
    if numel(m{c}) == prod(cfg.dim)
      stat.(sprintf('model%d',c)) = reshape(m{c},cfg.dim);
    end
  end
end

for c=1:size(desc,1) % iterate over parameter types
  stat.(sprintf('desc%d',c)) = desc{c};
end

% save crossvalidator object
stat.cv = cv;

% required
stat.trial = [];
