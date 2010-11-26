function [acc,sig,cv] = ft_mv_test(varargin)
% FT_MV_TEST tests a multivariate analysis procedure. It uses a default
% dataset about digit decoding using V1 BOLD response but this can be
% overridden. It handles different types of decoding problems.
%
%   [metric,sigtest,cv] = ft_test_mva(varargin)
%
% options:
% 'mva'       multivariate analysis
% 'X'         input data
% 'Y'         output data
% 'type'      decoding problem {'classification','regression','reconstruction'}
% 'nfolds'    number of crossvalidation folds
% 'parallel'  parallel cross-validation
% 'cv'        cross-validation object
% 'metric'    decoding metric
% 'sigtest'   significance test
% 
%
% EXAMPLE:
%   [a,b,c,d] = ft_mv_test('mva',{ft_mv_standardizer ft_mv_naive},'parallel',true,'type','regression');
%
% Copyright (c) 2009, Marcel van Gerven

  S = [];
  for i=1:2:length(varargin)
    S.(varargin{i}) = varargin{i+1};  
  end
        
  if ~isfield(S,'mva') && ~isfield(S,'cv')
    error('please specify mva or cv');
  end
  if ~isfield(S,'type')
    S.type = 'classification';
  end
  if ~isfield(S,'X') || ~isfield(S,'Y')
    
    fprintf('using default dataset 69digits\n');
    load 69digits;
    
    switch S.type
      case 'classification'
        S.X = X; % BOLD response
        S.Y = Y; % digit class
      case 'regression'
        S.X = I; % digit images
        S.Y = X(:,500); % BOLD response for one voxel
      case 'reconstruction'
        S.X = X; % BOLD response
        S.Y = I; % images

    end
  end
  
  if strcmp(S.type,'reconstruction') && size(S.Y,2) > 256
    % downsample images
    Z = S.Y;
    res = 16; % image resolution
    ores = sqrt(size(Z,2));
    I = zeros(size(Z,1),res^2);
    for j=1:size(Z,1)
      I(j,:) = reshape((imresize(double(reshape(Z(j,:),[ores ores])),[res res],'nearest')),[1 res^2]);
    end
    %S.Y = zscore(I);
  end
  
  if ~isfield(S,'nfolds')
    S.nfolds = 10;
  end
  if ~isfield(S,'parallel')
    S.parallel = false;
  end
  if ~isfield(S,'cv')
     S.cv = ft_mv_crossvalidator('parallel',S.parallel,'balance',false,'mva',S.mva,'nfolds',S.nfolds,'verbose',true,'compact',true,'init',1);
  end
  if ~isfield(S,'metric')
  if strcmp(S.type,'classification')
    S.metric = 'accuracy';
  elseif strcmp(S.type,'regression')
    S.metric = 'invresvar';
  else
    S.metric = 'correlation';
  end
  S.cv.metric = S.metric;
  if ~isfield(S,'sigtest')
    S.sigtest = [];
  end
  S.cv.sigtest = S.sigtest;
  
  cv      = S.cv.train(S.X,S.Y);
  acc     = cv.performance; 
  sig     = cv.significance;
  
end