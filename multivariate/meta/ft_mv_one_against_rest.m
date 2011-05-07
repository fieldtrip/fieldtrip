classdef ft_mv_one_against_rest < ft_mv_meta
%FT_MV_ONE_AGAINST_REST one-against-rest binary classification. 
%
%   This class evaluates a binary classifier on all possible pairs of class
%   labels. Borda count is used as default combfun. Can be overridden.
%
%   NOTE:
%   data is balanced such that if we have N trials for a target class then
%   we have at most N randomly selected trials for the other classes.
%
%   EXAMPLE:
%   ft_mv_one_against_rest('mva',ft_mv_svm);
%
%   Copyright (c) 2008, Marcel van Gerven

  properties
    
    nclasses % number of classes
    
    combfun = [];     % by default sums probabilities; can be overridden
    parallel = false; % run train function in parallel mode?


  end
  
  methods
    
    function obj = ft_mv_one_against_rest(varargin)
      
      obj = obj@ft_mv_meta(varargin{:});
          
      if isempty(obj.mva), error('mva not specified'); end
      
      if ~isa(obj.mva,'ft_mv_analysis'), obj.mva = ft_mv_analysis(obj.mva); end
      
    end
    
    function obj = train(obj,X,Y)
      
      if iscell(X)
        error('transfer learning with FT_MV_ONE_AGAINST_REST not supported');
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
        design{i} = nan(size(Y,1),1);
        
        ndidx = find(~didx);
        prm = randperm(numel(ndidx));
        design{i}(ndidx(prm(1:min(numel(prm),sum(didx)))),:) = 1;
        design{i}(didx,:)  = 2;
        
        data{i} = data{i}(~isnan(design{i}),:);
        design{i} = design{i}(~isnan(design{i}),:);
        
      end
      
      % alternative would be to resample the target class such that it is
      % balanced against the rest...
      
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
      for i=1:obj.nclasses
          
          if obj.verbose
            fprintf('testing class %d against rest\n',i);
          end
          
          Y{i} = obj.mva{i}.test(X);

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
      
       m = cell(length(obj.mva),1);
       desc = cell(length(obj.mva),1);
       for idx=1:length(obj.mva)
         mtd = obj.mva{idx}.mvmethods{end};
         [m{idx},desc{idx}] = mtd.model;
       end
      
    end
    
    
  end
end