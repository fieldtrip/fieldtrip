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
% 'validator' cross-validation object
% 'metric'    decoding metric
% 'sigtest'   significance test
% 'compact'   only keep model (true)
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
        
  if ~isfield(S,'mva') && ~isfield(S,'validator')
    error('please specify mva or validator');
  end
  if ~isfield(S,'type')
    S.type = 'classification';
  end
  if ~isfield(S,'compact')
    S.compact = true;
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
      %I(j,:) = reshape(im2bw(imresize(double(reshape(Z(j,:),[ores ores])),[res res],'nearest')),[1 res^2]);
      I(j,:) = reshape((imresize(double(reshape(Z(j,:),[ores ores])),[res res],'nearest')),[1 res^2]);
    end
    %S.Y = I+1;
    S.Y = I;
  end
  
  if ~isfield(S,'nfolds'), S.nfolds = 10; end
  
  if ~isfield(S,'parallel'), S.parallel = false; end
  
  if ~isfield(S,'validator')
     S.validator = ft_mv_crossvalidator('parallel',S.parallel,'balance',false,'mva',S.mva,'nfolds',S.nfolds,'verbose',true,'compact',S.compact,'init',1);
  end
  
  if ~isfield(S,'metric')
    if strcmp(S.type,'classification')
      S.metric = 'accuracy';
    elseif strcmp(S.type,'regression')
      S.metric = 'invresvar';
    else
      S.metric = 'correlation';
    end
  end
  
  S.validator.metric = S.metric;
  if ~isfield(S,'sigtest')
    S.sigtest = [];
  end
  S.validator.sigtest = S.sigtest;
  
  cv      = S.validator.train(S.X,S.Y);

  try
    acc     = cv.performance; 
    sig     = cv.significance;
  catch
    fprintf('incorrect metric/sigtest\n');
    acc     = nan;
    sig     = nan;
  end
  
  if strcmp(S.type,'reconstruction')

    design = cell2mat(cv.design');
    post = cell2mat(cv.post');
    
    for j=1:size(design,1)
      subplot(1,2,1);
      imagesc(reshape(design(j,:),[res res]));
      title('stimulus');
      colormap(gray);
      axis off; axis square
      subplot(1,2,2);
      imagesc(reshape(post(j,:),[res res]));
      title('prediction');
      colormap(gray);
      axis off; axis square
      drawnow
      pause
    end
    
  end
    
end