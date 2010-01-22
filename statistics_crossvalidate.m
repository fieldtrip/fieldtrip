function stat = statistics_crossvalidate(cfg, dat, design)
% STATISTICS_CROSSVALIDATE performs cross-validation using a prespecified
% classifier as given by cfg.classifier.
%
% Options:
% cfg.metric        = the metric to report (default = 'accuracy')
% cfg.cv            = crossvalidator object
%  overloads the following
%   cfg.classifier  = a classification procedure (default = myproc; below)
%   cfg.cvfolds     = number of folds (default = 10)
%   cfg.randomize   = permute the examples (default = true)
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
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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
   if ~isfield(cfg,'compact'), cfg.compact = false; end

   cv = crossvalidator('procedure',cfg.clfproc,'nfolds',cfg.nfolds,'compact',cfg.compact,'verbose',true);

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
[res,all] = cv.evaluate('metric',cfg.metric);
stat.prob = mean(cell2mat(all),1);

% is the statistic significant?
stat.significance = cv.significance();

% get the model wrt to each of the labels
if ~cfg.compact  

  m = cv.getmodel();

  if iscell(m) % transfer learning
    
    nlabels = size(m{1},1);
    
    if nlabels > 1
      for j=1:nlabels       
        stat.(sprintf('model%d',j)) = cellfun(@(x)(reshape(x(j,:),cfg.dim)),m,'UniformOutput',false);
      end
    else
      stat.model = cellfun(@(x)(reshape(x,cfg.dim)),m,'UniformOutput',false);
    end
    
  else
    
    nlabels = size(m,1);
    
    if nlabels > 1
      for j=1:nlabels
        stat.(sprintf('model%d',j)) = reshape(m(j,:),cfg.dim);
      end
    else
      stat.model = reshape(m,cfg.dim);
    end
  end
  
else % no model availabe 
  stat.model = [];  
end

% save crossvalidator object
stat.cv = cv;
