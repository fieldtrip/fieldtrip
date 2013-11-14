classdef one_against_one < dml.method
% ONE_AGAINST_ONE one-against-one binary classification. 
%
%   DESCRIPTION
%   This class evaluates a binary classifier on all possible pairs of class
%   labels. 
%
%   EXAMPLE
%   X = rand(15,20); Y = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3]';
%   m = dml.one_against_one('mva',dml.svm);
%   m = m.train(X,Y);
%   Z = m.test(X);
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)

  properties
    
    mva % used multivariate analysis
    
    nclasses % number of classes
    
    combfun = 'borda'; % 'majority' takes majority vote; 'borda' count sums outputs
    
    % balance the data such that both classes are evenly represented during
    % training (note that there may still be an imbalance during testing!)
    balance = false;
    
  end
  
  methods
    
    function obj = one_against_one(varargin)
      
      obj = obj@dml.method(varargin{:});
            
      if isempty(obj.mva), error('mva not specified'); end
      
      if ~isa(obj.mva,'dml.analysis'), obj.mva = dml.analysis(obj.mva); end
      
    end
    
    function obj = train(obj,X,Y)
      
      if iscell(X)
        error('transfer learning with ONE_AGAINST_ONE not yet supported');
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
      if strcmp(obj.combfun,'majority')
        
        Z = (Y{1} == repmat(max(Y{1},[],2),[1 size(Y{1},2)]));
        for j=2:length(Y)
          Z = Z + (Y{j} == repmat(max(Y{j},[],2),[1 size(Y{j},2)]));
        end
        Y = bsxfun(@rdivide,Z,sum(Z,2));
        
        
      elseif strcmp(obj.combfun,'borda')
        
        % borda count; rank candidates by summing and normalizing

        Z = Y{1};
        for j=2:length(Y)
          Z = Z + Y{j};
        end
        Y = bsxfun(@rdivide,Z,sum(Z,2));
        
      else 
      
       error('unknown combination function');
                
      end
      
    end
    
    function m = model(obj)
      % MODEL returns the following parameters:
      %
      % m.model{i} the model for the i-th binary classification problem
      % m.pair{i} the i-th pair of class labels belonging to the model
      
      m.model = cell((obj.nclasses*(obj.nclasses-1))/2,1);
      m.pair = cell((obj.nclasses*(obj.nclasses-1))/2,1);
      idx=1;
      for i=1:obj.nclasses
        for j=(i+1):obj.nclasses
          mtd = obj.mva{idx}.method{end};
          m.model{idx} = mtd.model;
          m.pair{idx} = [i j];
          idx = idx + 1;
        end
      end
      
    end
    
   
  end
end