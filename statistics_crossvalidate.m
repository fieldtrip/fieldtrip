function stat = statistics_crossvalidate(cfg, dat, design)

% STATISTICS_CROSSVALIDATE performs cross-validation using a prespecified
% multivariate analysis given by cfg.mva
%
% Options:
% cfg.metric        = the metric to report (default = 'accuracy')
% cfg.cv            = crossvalidator object
%  overloads the following
%   cfg.mva         = a multivariate analysis (default = {standardizer svmmethod})
%   cfg.nfolds      = number of folds (default = 10)
%   cfg.compact     = whether or not to save the mva procedure (true)
%   cfg.model       = whether or not to save the average model (true)
%
% Returns:
% stat.prob         = computed using the specified metric
% stat.significance = w.r.t. the null-hypothesis that we do no better than a
%                     majority classifier
% stat.cv           = the trained crossvalidator
%
% See also CROSSVALIDATE, MVA
%
% Requires: multivariate analysis toolbox

% Copyright (c) 2007, Marcel van Gerven
% F.C. Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

% specify classification procedure
if isfield(cfg,'cv')
  cv = cfg.cv;
else

  if ~isfield(cfg,'mva')
     cfg.mva = mva({ ...
        standardizer('verbose',true) ...
        svmmethod('verbose',true) ...
        });
   else
     if ~isa(cfg.mva,'mva')
       cfg.mva = mva(cfg.mva,'verbose',true);
     end
   end
   
   if ~isfield(cfg,'nfolds'), cfg.nfolds = 10; end
   if ~isfield(cfg,'compact'), cfg.compact = true; end
   if ~isfield(cfg,'model'), cfg.model = true; end

   for j=1:length(cfg.mva.mvmethods)
     cfg.mva.mvmethods{j}.indims = cfg.dim;
   end
   
   cv = crossvalidator('procedure',cfg.mva,'nfolds',cfg.nfolds,'compact',cfg.compact,'model',cfg.model,'verbose',true);

end

% specify metric returned in stat.prob
if ~isfield(cfg,'metric'), cfg.metric = 'accuracy'; end

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
cv = cv.validate(dat,design);

% the statistic of interest
res = cv.evaluate('metric',cfg.metric);
stat.prob = res;

% is the statistic significant?
stat.significance = cv.significance();

% get the models
if ~cfg.compact || cfg.model

  if cfg.model
    m = cv.model;
    desc = cv.desc;
  else
    [m,desc] = cv.getmodel();
  end
  
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
  
end

% save crossvalidator object
stat.cv = cv;
