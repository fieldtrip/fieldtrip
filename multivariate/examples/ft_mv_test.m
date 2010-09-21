function [acc,sig,tim,cv] = ft_mv_test(varargin)
% FT_MV_TEST tests multivariate analysis procedure
%
%   [acc,sig,tim,cv] = ft_test_mva(varargin)
%
%   EXAMPLE:
%    [a,b,c,d] = ft_mv_test('mva',{ft_mv_standardizer ft_mv_naive},'parallel',true);
%
%   Copyright (c) 2009, Marcel van Gerven

  S = [];
  for i=1:2:length(varargin)
    S.(varargin{i}) = varargin{i+1};  
  end
        
  if ~isfield(S,'mva') && ~isfield(S,'cv')
    error('please specify mva or cv');
  end
  if ~isfield(S,'X') || ~isfield(S,'Y')
    fprintf('using default dataset 69digits\n');
    load 69digits;
    S.X = X;
    S.Y = Y;    
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
    S.metric = 'accuracy';
  end
  S.cv.metric = S.metric;
  if ~isfield(S,'sigtest')
    S.sigtest = 'mcnemar';
  end
  S.cv.sigtest = S.sigtest;
  
  tstart  = tic;
  cv      = S.cv.train(S.X,S.Y);
  tim     = toc(tstart);
  acc     = cv.performance; 
  sig     = cv.significance;
  
end