classdef one_against_rest < classifier
%ONE_AGAINST_REST one-against-rest binary classification
%
%   This class evaluates a binary classifier on all a class label and the
%   pooled data for the other class labels
%
%   EXAMPLE:
%   
%   myproc = mva({ ...
%        preprocessor('prefun',@(x)(log10(x))) ...
%        standardizer() ... 
%        one_against_rest('procedures',mva({gp()})) ...
%        });
%
%   SEE ALSO:
%   one_against_one
%
%   Copyright (c) 2008, Marcel van Gerven

    properties        
      
      procedure = mva({da()}); % the used classification procedure
  
    end
    
    methods
      function obj = one_against_rest(varargin)
        
        obj = obj@classifier(varargin{:});
        
        if ~isa(obj.procedure,'mva')
          obj.procedure = mva(obj.procedure);
        end
        
      end
      
      function p = estimate(obj,X,Y)
        
        if iscell(X)
          error('transfer learning with one_against_rest not yet supported');
        end
        
         p.nclasses = obj.nunique(Y);
                
        % transform the data such that we have a cell element for each
        % class label pair
        
        data = cell(1,p.nclasses);
        design = cell(1,p.nclasses);
        
        % create new data representation
        
        for i=1:p.nclasses
          
          didx = (Y == i);
          
          data{i} = X;
          design{i} = Y;
          design{i}(~didx,:) = 1;
          design{i}(didx,:)  = 2;
          
        end
        
        % replicate the classifier
        if ~iscell(obj.procedure)
          procedure = obj.procedure;
          p.procedure = cell(1,length(data));
          for j=1:length(data)
            p.procedure{j} = procedure;
          end
        else
          p.procedure = obj.procedure;
        end
        
        for j=1:length(data)
          
          if obj.verbose
            fprintf('training class %d against rest\n',j);
          end
          
          p.procedure{j} = p.procedure{j}.train(data{j},design{j});
        end
        
      end
      
      function Y = map(obj,X)
             
        % get posteriors for all pairs
        Y = zeros(size(X,1),obj.params.nclasses);
        for i=1:obj.params.nclasses
          
          if obj.verbose
            fprintf('testing class %d against rest\n',i);
          end
          
          tpost = obj.params.procedure{i}.test(X);
          Y(:,i) = tpost(:,2);
        end
        
        % normalize
        Y = bsxfun(@rdivide,Y,sum(Y,2));
        
      end
      
       function [m,desc] = getmodel(obj)
       
        m=[];
        desc=[];
        for c=1:length(obj.params.procedure)
          
          mtd = obj.params.procedure{c}.mvmethods{end};
          [mm,dd] = mtd.getmodel();
          
          if isempty(m)
            m = cell(size(mm,1)*length(obj.params.procedure),1);
            desc = cell(size(mm,1)*length(obj.params.procedure),1);
          end
          
          m((c-1)*size(mm,1)+(1:size(mm,1))) = mm;
          
          for k=1:size(mm,1)
            desc{(c-1)*size(mm,1)+k} = sprintf('class %d against rest; %s\n',c,dd{k});
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