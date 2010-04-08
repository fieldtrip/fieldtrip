classdef basis_recon < reconstructor
% BASIS_RECON basis vector reconstructor class
%
% A basis may be prespecified or it may be learned from the design.
% 
% Prediction is realized by obj.procedure which tries to predict the design
% matrix in its new basis. Afterwards, the prediction is inverted using
% obj.A
%
% the design matrix is standardized prior to transformation to the new basis 
% using and obj.W and unstandardized after prediction; this should be
% disabled when using a classifier as the predicting procedure
%
% Basis transformations may also be applied to the input data. This is
% simply realized by adding the transform to the mva pipeline.
%
% By default, basis_recon uses the identity basis; i.e., it learns to
% predict individual elements of the design matrix using a ridge regressor
%
% EXAMPLES:
%
% load dataset; X = dataset(response); 
% Y = zeros(size(stimuli,1),256); for j=1:size(stimuli,1), Y(j,:) = ...
% reshape(imresize(reshape(stimuli(j,:),[28 28]),[16 16]),[1 256]); end
% Y = dataset(Y);
% 
% % ridge regression on individual elements of the design matrix (i.e., pixel level reconstructions)
% % note however that the ridge regressor already can handle multiple
% % outputs
% p = mva({standardizer basis_recon('procedure',{linreg('L2',1e-2)},'verbose',true)});
% p = p.train(X(1:100,:),Y(1:100,:));
% r = p.test(X(101:111,:));
% 
% % class label probabilities using l2 regularized logistic regression
% Y = dataset((Y.X > 128) + 1); % threshold data
% p = mva({standardizer basis_recon('procedure',{logreg('L2',1e-2)},'verbose',true,'standardize',false)});
% p = p.train(X(1:100,:),Y(1:100,:));
% r = p.test(X(101:111,:));
% 
% % prelearned ICA basis functions
% load ica; 
% p = mva({standardizer basis_recon('A',A,'W',W,'procedure',{linreg('L2',1e-2)},'verbose',true)});
% p = p.train(X(1:100),Y(1:100));
% r = p.test(X(101:111));
% 
% % visualization
% images(r,1:10,[2 5]);
% figure
% images(Y,101:110,[2 5]);
%
% Copyright (c) 2010, Marcel van Gerven


  properties
  
    procedure={linreg('L2',1e-2)}; % the procedure used for prediction of each source
    
    W; % mapping from input to sources; nsources * nfeatures; whitening included in W
    A; % mapping from sources to outputs; nfeatures * nsources; unwhitening included in A
    
    % standardization design before making the basis transformation
    standardize = true;
       
  end

  methods
    
    function obj = basis_recon(varargin)
      
      obj = obj@reconstructor(varargin{:});
      
      assert(~isempty(obj.procedure));
      
      % cast to mva
      obj.procedure = mva(obj.procedure);
      
    end
    
    function p = estimate(obj,X,Y)
      
      % standardize design
      if obj.standardize
        % standardize the design
        p.prepout = mva({standardizer});
        p.prepout = p.prepout.train(Y);
        Y = p.prepout.test(Y);
      end
      
      % compute basis if not prespecified
      if isempty(obj.W)
        [p.W,p.A] = obj.learn_basis(Y);
      else
        p.W = obj.W;
        p.A = obj.A;
      end
       
      % move from reconstructions to sources
      S = Y * p.W';
      
      % predict the sources
      p.proc = cell(1,size(p.W,1));
      for c=1:length(p.proc)
        
        if obj.verbose
          fprintf('training source %d of %d\n',c,length(p.proc));
        end
        
        p.proc{c} = obj.procedure;
        p.proc{c} = p.proc{c}.train(X,S(:,c));
        
      end      
      
    end
    
    function Y = map(obj,X)
      % also returns the estimated sources
      
      % changes source estimate back to native space
      Y = obj.predict_sources(X) * obj.params.A';
      
      Y = reshape(Y,[size(X,1) obj.outdims]);
      
      % invert the standardization
      if obj.standardize
        Y = obj.params.prepout.untest(Y);
      end
      
    end
    
    function S = predict_sources(obj,data)
      % predict the sources based on clf procedure
      
      proc = obj.params.proc;
      
      nsources = size(obj.params.W,1);
      
      S = zeros(size(data,1),nsources);
      
      for c=1:nsources
        if obj.verbose
          fprintf('testing source %d of %d\n',c,length(proc));
        end
        pp = proc{c}.test(data);
        S(:,c) = pp(:,1); % return the first column only (in case of classification)
      end
      
    end
    
  end
  
  methods(Static=true)
  
    function [W,A] = learn_basis(design)
      % called when basis is unspecified; here, we simply
      % return the identity (giving independent pixel reconstruction)
      % child classes may override this behaviour
      
      W = eye(size(design,2));
      A = eye(size(design,2));
      
    end
    
  end
  
end
