classdef one_against_one < classifier
%ONE_AGAINST_ONE one-against-one binary classification
%
%   This class evaluates a binary classifier on all possible pairs of class
%   labels
%
%   EXAMPLE:
%   
%   myproc = clfproc({ ...
%        preprocessor('prefun',@(x)(log10(x))) ...
%        standardizer() ... 
%        one_against_one('procedure',clfproc({gp()})) ...
%        });
%
%   SEE ALSO:
%   ensemble.m
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: one_against_one.m,v $
%

    properties        
    
      procedure = []; % clfproc({da()}); % the used classification procedures
      combination = 'product'; % how to combine classifier output (not how to combine data)
    
      nclasses;
      
    end
    
    methods
        
      function obj = one_against_one(varargin)

        obj = obj@classifier(varargin{:});
        
        assert(~isempty(obj.procedure));
        
        % cast to clf procedure
        if ~isa(obj.procedure,'clfproc')
          obj.procedure = clfproc(obj.procedure);
        end
            
      end
      
      function obj = train(obj,tdata,tdesign)
          
        % transform the data such that we have a cell element for each
        % class label pair
        obj.nclasses = tdesign.nunique;
        
        % replicate the classifier
        
        ncomp = obj.nclasses*(obj.nclasses-1)/2;
        
        if ~iscell(obj.procedure)
          procedure = obj.procedure;
          obj.procedure = cell(1,ncomp);
          for j=1:ncomp
            obj.procedure{j} = procedure;
          end
        end
        
       
        idx = 1;
        for i=1:obj.nclasses
          for j=(i+1):obj.nclasses
           
            if obj.verbose
              fprintf('training class %d against class %d\n',i,j);
            end
            
            didx = (tdesign.X == i | tdesign.X == j);
            
            design = tdesign.X(didx,:);
            design(design == i) = 1;
            design(design == j) = 2;
            
            obj.procedure{idx} = obj.procedure{idx}.train(tdata.subsample(didx),dataset(design));
        
            idx=idx+1;
            
          end
        end
        
      end
      
      function post = test(obj,data)
        
        cpost = cell(1,length(obj.procedure));
        
        % get posteriors for all pairs
        idx = 1;
        for i=1:obj.nclasses
          for j=(i+1):obj.nclasses
            
            if strcmp(obj.combination,'product')
              % use ones to allow for products of probabilities            
              cpost{idx} = ones(data.nsamples,obj.nclasses);
            else
              cpost{idx} = zeros(data.nsamples,obj.nclasses);
            end
            
            if obj.verbose
              fprintf('testing class %d against class %d\n',i,j);
            end
            
            p = obj.procedure{idx}.test(data);
            
            cpost{idx}(:,[i j]) = p.X;
            
            cpost{idx} = dataset(cpost{idx});
            
            idx = idx+1;
            
          end
        end
        
        % combine the result
        post = combiner.combine(cpost,obj.combination);
        
      end
      
    end
end