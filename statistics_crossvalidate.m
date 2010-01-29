function stat = statistics_crossvalidate(cfg, dat, design)
% STATISTICS_CROSSVALIDATE performs cross-validation using a prespecified
% classifier as given by cfg.classifier.
%
% Options:
% cfg.metric        = the metric to report (default = 'accuracy')
% cfg.cv            = crossvalidator object
%  overloads the following
%   cfg.clfproc     = a classification procedure (default = {standardizer svmmethod})
%   cfg.nfolds      = number of folds (default = 10)
%   cfg.compact     = whether or not to save the classification procedure (true)
%   cfg.model       = whether or not to save the average model (true)
%
% Returns:
% stat.prob         = computed using the specified metric
% stat.significance = w.r.t. the null-hypothesis that we do no better than a
%                     majority classifier
% stat.cv           = the trained crossvalidator
%
% See also: CROSSVALIDATE, CLFPROC, CLASSIFIER
%
% Requires: classification toolbox
%
% Copyright (c) 2007, Marcel van Gerven
% F.C. Donders Centre for Cognitive Neuroimaging, Nijmegen, NL

fieldtripdefs

% specify classification procedure
if isfield(cfg,'cv')
  cv = cfg.cv;
else

  if ~isfield(cfg,'clfproc')
     cfg.clfproc = clfproc({ ...
        standardizer('verbose',true) ...
        svmmethod('verbose',true) ...
        });
   else
     if ~isa(cfg.clfproc,'clfproc')
       cfg.clfproc = clfproc(cfg.clfproc,'verbose',true);
     end
   end
   
   if ~isfield(cfg,'nfolds'), cfg.nfolds = 10; end
   if ~isfield(cfg,'compact'), cfg.compact = true; end
   if ~isfield(cfg,'model'), cfg.model = true; end

   cv = crossvalidator('procedure',cfg.clfproc,'nfolds',cfg.nfolds,'compact',cfg.compact,'model',cfg.model,'verbose',true);

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
    dat{c}    = dataset(tmpdat(cfg.dataset == c,:));
    design{c} = dataset(tmpdesign(cfg.dataset == c,:));
  end
  
else
  
  dat = dataset(dat');
  design = dataset(design');

end

% perform everything ;o)
cv = cv.validate(dat,design);

% the statistic of interest
res = cv.evaluate('metric',cfg.metric);
stat.prob = res;

% is the statistic significant?
stat.significance = cv.significance();

% get the models
if ~cfg.compact  

  [m,desc] = cv.getmodel();

  for c=1:size(m,1) % iterate over parameter types
    
    if size(m,2) > 1 % transfer learning
      
      for j=1:size(m,2)        
        stat.(sprintf('model%d_%d',c,j)) = reshape(m{c,j},cfg.dim);
      end
            
    else      
      stat.(sprintf('model%d',c)) = reshape(m{c},cfg.dim);
    end
  end  
  
  for c=1:size(desc,1) % iterate over parameter types
    stat.(sprintf('desc%d',c)) = desc{c};
  end
  
end

% save crossvalidator object
stat.cv = cv;
