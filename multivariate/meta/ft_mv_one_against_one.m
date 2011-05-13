classdef ft_mv_one_against_one < ft_mv_meta
%FT_MV_ONE_AGAINST_ONE one-against-one binary classification. 
%
%   This class evaluates a binary classifier on all possible pairs of class
%   labels. Borda count is used as default combfun. Can be overridden.
%
%   EXAMPLE:
%   ft_mv_one_against_one('mva',ft_mv_svm);
%
%   Copyright (c) 2008, Marcel van Gerven

  properties
    
    nclasses % number of classes
    
    combfun = [];     % by default sums probabilities; can be overridden
    parallel = false; % run train function in parallel mode?

    % balance the data such that both classes are evenly represented during
    % training (note that there may still be an imbalance during testing!)
    balance = false;
    
  end
  
  methods
    
    function obj = ft_mv_one_against_one(varargin)
      
      obj = obj@ft_mv_meta(varargin{:});
            
      if isempty(obj.mva), error('mva not specified'); end
      
      if ~isa(obj.mva,'ft_mv_analysis'), obj.mva = ft_mv_analysis(obj.mva); end
      
    end
    
    function obj = train(obj,X,Y)
      
      if iscell(X)
        error('transfer learning with FT_MV_ONE_AGAINST_ONE not yet supported');
      end
      
      % transform the data such that we have a cell element for each
      % class label pair
      obj.nclasses = max(Y);
      
      % replicate the classifier
      ncomp = obj.nclasses*(obj.nclasses-1)/2;
      
      if ~iscell(obj.mva)
        procedure = obj.mva;
        obj.mva = cell(1,ncomp);
        for j=1:ncomp
          obj.mva{j} = procedure;
        end
      else
        obj.mva = obj.mva;
      end
      
      idx = 1;
      for i=1:obj.nclasses
        for j=(i+1):obj.nclasses
          
          if obj.verbose
            fprintf('training class %d against class %d (pair %d of %d)\n',i,j,idx,floor(0.5*obj.nclasses*(obj.nclasses-1)));
          end
          
          if obj.balance % balance the data between classes

            didx1 = find(Y==i);
            didx2 = find(Y==j);
            
            ntrials = min([numel(didx1) numel(didx2)]);
            
            prm = randperm(numel(didx1));
            didx1 = didx1(prm(1:ntrials));
            
            prm = randperm(numel(didx2));
            didx2 = didx2(prm(1:ntrials));
            
            didx = [didx1; didx2];
            
            if obj.verbose
              fprintf('balancing data using %d trials per class\n',ntrials);
            end
            
          else
            
            didx = (Y == i | Y == j);          
          end
          
          design = Y(didx,:);
          design(design == i) = 1;
          design(design == j) = 2;
          
          obj.mva{idx} = obj.mva{idx}.train(X(didx,:),design);
          
          idx=idx+1;
          
        end
      end
      
    end
    
    function Y = test(obj,X)
      
      Y = cell(1,length(obj.mva));
      
      % get posteriors for all pairs
      idx = 1;
      for i=1:obj.nclasses
        for j=(i+1):obj.nclasses
          
          Y{idx} = zeros(size(X,1),obj.nclasses);
          
          if obj.verbose
            fprintf('testing class %d against class %d\n',i,j);
          end
          
          p = obj.mva{idx}.test(X);
          
          Y{idx}(:,[i j]) = p;
          
          idx = idx+1;
          
        end
      end
      
      % combine the result
      if ~isempty(obj.combfun)
      
        Y = obj.combfun(Y);
        
        % a better combination rule might be voting (to avoid influences of
        % spuriously high probabilities
      
      else
        % sum and normalize probabilities

        Z = Y{1};
        for j=2:length(Y)
          Z = Z + Y{j};
        end
        Y = bsxfun(@rdivide,Z,sum(Z,2));
        
      end
      
    end
    
    function [m,desc] = model(obj)
      
      m = cell((obj.nclasses*(obj.nclasses-1))/2,1);
      desc = cell((obj.nclasses*(obj.nclasses-1))/2,1);
      idx=1;
      for i=1:obj.nclasses
        for j=(i+1):obj.nclasses
          mtd = obj.mva{idx}.mvmethods{end};
          [m{idx},desc{idx}] = mtd.model;
          idx = idx + 1;
        end
      end
      
    end
    
   
  end
end