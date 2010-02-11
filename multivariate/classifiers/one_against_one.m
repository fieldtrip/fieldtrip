classdef one_against_one < classifier
%ONE_AGAINST_ONE one-against-one binary classification
%
%   This class evaluates a binary classifier on all possible pairs of class
%   labels
%
%   EXAMPLE:
%   
%   myproc = mva({ ...
%        preprocessor('prefun',@(x)(log10(x))) ...
%        standardizer() ... 
%        one_against_one('procedure',mva({gp()})) ...
%        });
%
%   SEE ALSO:
%   ensemble.m
%
%   Copyright (c) 2008, Marcel van Gerven


    properties        
    
      procedure = []; % mva({da()}); % the used classification procedures
      combination = 'product'; % how to combine classifier output (not how to combine data)
   
    end
    
    methods
        
      function obj = one_against_one(varargin)

        obj = obj@classifier(varargin{:});
        
        assert(~isempty(obj.procedure));
        
        % cast to clf procedure
        if ~isa(obj.procedure,'mva')
          obj.procedure = mva(obj.procedure);
        end
            
      end
      
      function p = estimate(obj,X,Y)
          
        if iscell(X)
          error('transfer learning with one_against_one not yet supported');
        end
        
        % transform the data such that we have a cell element for each
        % class label pair
        p.nclasses = obj.nunique(Y);
        
        % replicate the classifier
        
        ncomp = p.nclasses*(p.nclasses-1)/2;
        
        if ~iscell(obj.procedure)
          procedure = obj.procedure;
          p.procedure = cell(1,ncomp);
          for j=1:ncomp
            p.procedure{j} = procedure;
          end
        else
          p.procedure = obj.procedure;
        end
        
        idx = 1;
        for i=1:p.nclasses
          for j=(i+1):p.nclasses
           
            if obj.verbose
              fprintf('training class %d against class %d\n',i,j);
            end
            
            didx = (Y == i | Y == j);
            
            design = Y(didx,:);
            design(design == i) = 1;
            design(design == j) = 2;
            
            p.procedure{idx} = p.procedure{idx}.train(X(didx,:),design);
        
            idx=idx+1;
            
          end
        end
        
      end
      
      function Y = map(obj,X)
        
        Y = cell(1,length(obj.params.procedure));
        
        % get posteriors for all pairs
        idx = 1;
        for i=1:obj.params.nclasses
          for j=(i+1):obj.params.nclasses
            
            if strcmp(obj.combination,'product')
              % use ones to allow for products of probabilities            
              Y{idx} = ones(size(X,1),obj.params.nclasses);
            else
              Y{idx} = zeros(size(X,1),obj.params.nclasses);
            end
            
            if obj.verbose
              fprintf('testing class %d against class %d\n',i,j);
            end
            
            p = obj.params.procedure{idx}.test(X);
            
            Y{idx}(:,[i j]) = p;
            
            idx = idx+1;
            
          end
        end
        
        % combine the result
        Y = combiner.combine(Y,obj.combination);
        
      end
      
      function [m,desc] = getmodel(obj)
        
        m=[];
        desc=[];
        c=1;
        for i=1:obj.params.nclasses
          for j=(i+1):obj.params.nclasses
            
            mtd = obj.params.procedure{c}.mvmethods{end};
            [mm,dd] = mtd.getmodel();
            
            if isempty(m)
              m = cell(size(mm,1)*length(obj.params.procedure),1);
              desc = cell(size(mm,1)*length(obj.params.procedure),1);
            end
            
            m((c-1)*size(mm,1)+(1:size(mm,1))) = mm;
            
            for k=1:size(mm,1)
              desc{(c-1)*size(mm,1)+k} = sprintf('class %d against %d; %s\n',i,j,dd{k});
            end
            
            c = c+1;
          end
          
        end
        
        
      end
      
      function b = istransfer(obj)
        % return whether or not this method is a transfer learner
        % must be overloaded by e.g., one_against_one
        
        b =  obj.procedure.mvmethods{end}.istransfer();
      end
      
    end
end