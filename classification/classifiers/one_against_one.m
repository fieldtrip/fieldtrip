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
        
        if ~isa(obj.procedure,'clfproc')
          obj.procedure = clfproc(obj.procedure);
        end
            
      end
      
      function obj = train(obj,tdata,tdesign)
          
        obj.nclasses = tdesign.nunique;
        
        tdata = tdata.collapse();
        tdesign = tdesign.collapse();
        
        % transform the data such that we have a cell element for each
        % class label pair
        nclasses = obj.nclasses;
        
        data = cell(1,nclasses*(nclasses-1)/2);
        design = cell(size(data));
        
        % create new data representation
        
        idx = 1;
        for i=1:nclasses
          for j=(i+1):nclasses
            
            didx = (tdesign == i | tdesign == j);
            
            data{idx} = tdata(didx,:);
            design{idx} = tdesign(didx,:);
            design{idx}(design{idx} == i) = 1;
            design{idx}(design{idx} == j) = 2;
            
            idx = idx+1;
          end
        end
        
        
        % replicate the classifier
        if ~iscell(obj.procedure)
          procedure = obj.procedure;
          obj.procedure = cell(1,length(data));
          for j=1:length(data)
            obj.procedure{j} = procedure;
          end
        end
        
        for j=1:length(data)
          obj.procedure{j} = obj.procedure{j}.train(dataset(data{j}),dataset(design{j}));
        end
        
      end
      
      function post = test(obj,data)
        
        nclasses = obj.nclasses;
        
        cpost = cell(1,nclasses*(nclasses-1)/2);
        
        % get posteriors for all pairs
        idx = 1;
        for i=1:nclasses
          for j=(i+1):nclasses
            
            if strcmp(obj.combination,'product')
              % use ones to allow for products of probabilities            
              cpost{idx} = ones(data.nsamples,obj.nclasses);
            else
              cpost{idx} = zeros(data.nsamples,obj.nclasses);
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