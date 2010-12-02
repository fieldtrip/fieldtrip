classdef ft_mv_one_against_rest < ft_mv_predictor
%FT_MV_ONE_AGAINST_REST one-against-rest binary classification. 
%
%   This class evaluates a binary classifier on all possible pairs of class
%   labels. Borda count is used as default combfun. Can be overridden.
%
%   EXAMPLE:
%   ft_mv_one_against_rest('mva',ft_mv_svm);
%
%   Copyright (c) 2008, Marcel van Gerven

  properties
    
    nclasses % number of classes
    
    mva               % parallel mva mva to run
    combfun = [];     % by default sums probabilities; can be overridden
    parallel = false; % run train function in parallel mode?


  end
  
  methods
    
    function obj = ft_mv_one_against_rest(varargin)
      
      obj = obj@ft_mv_predictor(varargin{:});
      
      if ~isa(obj.mva,'ft_mv_analysis')
        obj.mva = ft_mv_analysis(obj.mva);
      end
      
    end
    
    function obj = train(obj,X,Y)
      
      if iscell(X)
        error('transfer learning with FT_MV_ONE_AGAINST_REST not yet supported');
      end
      
      % transform the data such that we have a cell element for each
      % class label pair
      obj.nclasses = max(Y);
      
      % transform the data such that we have a cell element for each
      % class label pair
      
      data = cell(1,obj.nclasses);
      design = cell(1,obj.nclasses);
      
      % create new data representation
      
      for i=1:obj.nclasses
        
        didx = (Y == i);
        
        data{i} = X;
        design{i} = Y;
        design{i}(~didx,:) = 1;
        design{i}(didx,:)  = 2;
        
      end
      
      % replicate the classifier
      if ~iscell(obj.mva)
        procedure = obj.mva;
        obj.mva = cell(1,length(data));
        for j=1:length(data)
          obj.mva{j} = procedure;
        end
      else
        obj.mva = obj.mva;
      end
      
      for j=1:length(data)
        
        if obj.verbose
          fprintf('training class %d against rest\n',j);
        end
        
        obj.mva{j} = obj.mva{j}.train(data{j},design{j});
      end
      
    end
    
    function Y = test(obj,X)
      
      % get posteriors for all pairs
      idx = 1;
      for i=1:obj.nclasses
        for j=(i+1):obj.nclasses
          
          
          if obj.verbose
            fprintf('testing class %d against class %d\n',i,j);
          end
          
          Y{idx} = obj.mva{idx}.test(X);
                    
          idx = idx+1;
          
        end
      end
      
      % combine the result
      if ~isempty(obj.combfun)
      
        Y = obj.combfun(Y);
      
      else
    
        Z = zeros(size(Y{1},1),length(Y));
        for j=1:length(Y)
          Z(:,j) = Y{j}(:,2);
        end
        % normalize
        Y = bsxfun(@rdivide,Z,sum(Z,2));
      
      end
      
    end
    
    function [m,desc] = model(obj)
      
      m=[];
      desc=[];
      for c=1:length(obj.mva)
        
        mtd = obj.mva{c}.mvmethods{end};
        [mm,dd] = mtd.model();
        
        if isempty(m)
          m = cell(size(mm,1)*length(obj.mva),1);
          desc = cell(size(mm,1)*length(obj.mva),1);
        end
        
        m((c-1)*size(mm,1)+(1:size(mm,1))) = mm;
        
        for k=1:size(mm,1)
          desc{(c-1)*size(mm,1)+k} = sprintf('class %d against rest; %s\n',c,dd{k});
        end
      end
      
    end
    
    
  end
end