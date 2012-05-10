classdef one_against_rest < dml.method
% ONE_AGAINST_REST one-against-rest binary classification. 
%
%   DESCRIPTION
%   This class evaluates a binary classifier on all possible pairs of class
%   labels. Borda count is used as default combfun. Can be overridden.
%
%   EXAMPLE
%   X = rand(15,20); Y = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3]';
%   m = dml.one_against_rest('mva',dml.svm);
%   m = m.train(X,Y);
%   Z = m.test(X);
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)

  properties
    
    mva % used multivariate analysis

    nclasses % number of classes
    
    combfun = 'borda' % 'borda' borda count
    
  end
  
  methods
    
    function obj = one_against_rest(varargin)
      
      obj = obj@dml.method(varargin{:});
          
      if isempty(obj.mva), error('mva not specified'); end
      
      if ~isa(obj.mva,'ft_mv_analysis'), obj.mva = dml.analysis(obj.mva); end
      
    end
    
    function obj = train(obj,X,Y)
      
      if iscell(X)
        error('transfer learning with ONE_AGAINST_REST not supported');
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
        
        design{i} = nan(size(Y,1),1);
        
        didx = (Y == i);
        ndidx = find(~didx);

        % downsample
        rand('seed',1); randn('seed',1);
        prm = randperm(numel(ndidx));
        design{i}(ndidx(prm(1:min(numel(prm),sum(didx)))),:) = 1;
        design{i}(didx,:)  = 2;
        data{i} = X;
        data{i} = data{i}(~isnan(design{i}),:);
        design{i} = design{i}(~isnan(design{i}),:);
         
        % upsample
%         didx = find(didx);
%         didx = repmat(didx,[ceil(numel(ndidx)/numel(didx)) 1]);
%         didx = didx(1:numel(ndidx));
%         design{i}(1:numel(didx)) = 1;
%         design{i}(numel(didx)+(1:numel(ndidx))) = 2;
%         data{i} = X(cat(1,didx,ndidx),:);
        
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
      for i=1:obj.nclasses
          
          if obj.verbose
            fprintf('testing class %d against rest\n',i);
          end
          
          Y{i} = obj.mva{i}.test(X);

      end
      
      % combine the result
      if strcmp(obj.combfun,'borda')
        
        % borda count
        
        Z = zeros(size(Y{1},1),length(Y));
        for j=1:length(Y)
          Z(:,j) = Y{j}(:,2);
        end
        % normalize
        Y = bsxfun(@rdivide,Z,sum(Z,2));
       
      else 
      
       error('unknown combination function');
                
       
      end
      
    end
    
    function m = model(obj)
      % MODEL returns the following parameters:
      %
      % m.model{i} the model for the i-th binary classification problem
      
      m.model = cell(length(obj.mva),1);
      for idx=1:length(obj.mva)
        mtd = obj.mva{idx}.method{end};
        m.model{idx} = mtd.model;
      end
      
    end
    
    
  end
end