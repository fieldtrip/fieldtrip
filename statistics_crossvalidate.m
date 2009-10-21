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
% $Log: statistics_crossvalidate.m,v $
% Revision 1.28  2009/07/04 10:32:05  marvger
% changed transfer_classifier to transfer_learner
%
% Revision 1.27  2009/05/19 10:48:56  marvger
% added cfg.compact facility
%
% Revision 1.26  2009/02/18 15:13:38  marvger
% changed default to svmmethod since kernelmethod crashes on some mentats
%
% Revision 1.25  2009/02/17 09:28:47  marvger
% changed default from svmmethod to kernelmethod
%
% Revision 1.24  2009/02/16 18:51:02  marvger
% changed default classifier to kernelmethod
%
% Revision 1.23  2009/02/16 18:39:38  marvger
% regular update
%
% Revision 1.22  2009/01/29 15:43:47  marvger
% added stat.modelx to represent classifier parameters
%
% Revision 1.21  2009/01/29 13:22:08  marvger
% general update; functionality is now stable and useable
%
% Revision 1.20  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.19  2008/09/18 14:23:39  marvger
% integrated classification_module representation
%
% Revision 1.18  2007/12/19 08:42:57  marvger
% changed naming
%
% Revision 1.17  2007/11/23 13:20:12  marvger
% swap dimensions from features x repetitions to repetitions x features
%
% Revision 1.16  2007/11/13 16:14:01  marvger
% crossvalidate is now part of the classification toolbox
%
% Revision 1.15  2007/11/12 18:46:19  marvger
% *** empty log message ***
%
% Revision 1.14  2007/10/10 10:22:07  marvger
% *** empty log message ***

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
