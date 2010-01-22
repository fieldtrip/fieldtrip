classdef one_against_rest < classifier
%ONE_AGAINST_REST one-against-rest binary classification
%
%   This class evaluates a binary classifier on all a class label and the
%   pooled data for the other class labels
%
%   EXAMPLE:
%   
%   myproc = clfproc({ ...
%        preprocessor('prefun',@(x)(log10(x))) ...
%        standardizer() ... 
%        one_against_rest('procedures',clfproc({gp()})) ...
%        });
%
%   SEE ALSO:
%   one_against_one
%   ensemble.m
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: one_against_rest.m,v $
%

    properties        
      
      procedure = clfproc({da()}); % the used classification procedure
      nclasses;
      
    end
    
    methods
      function obj = one_against_rest(varargin)
        
        obj = obj@classifier(varargin{:});
        
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
        
        data = cell(1,obj.nclasses);
        design = cell(size(data));
        
        % create new data representation
        
        for i=1:obj.nclasses
          
          didx = (tdesign == i);
          
          data{i} = tdata;
          design{i} = tdesign;
          design{i}(~didx,:) = 1;
          design{i}(didx,:)  = 2;
          
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
             
        % get posteriors for all pairs
        post = zeros(data.nsamples,obj.nclasses);
        for i=1:obj.nclasses
          
          tpost = obj.procedure{i}.test(data);
          post(:,i) = tpost.X(:,2);
        end
        
        % normalize
        post = dataset(post ./ repmat(sum(post,2),[1 size(post,2)]));
        
      end
      
    end
end