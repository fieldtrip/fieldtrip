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
   
   if ~isfield(cfg,'cvfolds'), cfg.cvfolds = 10; end
   if ~isfield(cfg,'randomize'), cfg.randomize = true; end
   if ~isfield(cfg,'compact'), cfg.compact = false; end

   cv = crossvalidator('procedure',cfg.clfproc,'cvfolds',cfg.cvfolds,'randomize',cfg.randomize,'compact',cfg.compact,'verbose',true);

end

% specify metric returned in stat.prob
if ~isfield(cfg,'metric'), cfg.metric = 'accuracy'; end

% check for transfer learning; this is implemented by a design matrix where
% the second vector labels the subjects

if isa(cfg.clfproc.clfmethods{end},'transfer_learner') && size(design,1) == 2

  % split up the datasets
  tmpdat = dat';
  tmpdesign = design';
  n = max(tmpdesign(:,2));
  dat = cell(1,n);
  design = cell(1,n);
  for c=1:n
    dat{c}    = tmpdat(tmpdesign(:,2) == c,:);
    design{c} = tmpdesign(tmpdesign(:,2) == c,1);
  end
  
  % perform everything ;o)
  cv = cv.validate(dat,design);

else

  % perform everything ;o)
  cv = cv.validate(dat',design');
end

% the statistic of interest
stat.prob = cv.evaluate('metric',cfg.metric);

% is the statistic significant?
stat.significance = cv.significance();

% get the model wrt to each of the class labels
if ~cfg.compact
  
  if ~isnan(cv.nclasses())
    for j=1:cv.nclasses()
      stat.(sprintf('model%d',j)) = cv.getmodel(j,cfg.dim);
    end
  else
    stat.model = cv.getmodel(nan,cfg.dim);
  end
  
else % no model availabe 
  stat.model = [];  
end

% save crossvalidator object
stat.cv = cv;
