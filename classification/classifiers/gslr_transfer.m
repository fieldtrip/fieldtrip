classdef gslr_transfer < classifier & transfer_learner
%GSLR_TRANSFER group sparsifying logistic regression for transfer learning
%
%   Group sparsifying logistic regression can be used to regularize a
%   logistic regression model with respect to feature groups. As a special
%   case, for multiple datasets, we are performing transfer learning. This
%   is useful for MEG/EEG analysis when we combine the data of multiple
%   subjects/sessions. 

%   SEE ALSO:
%   gslr.m
%   regularization_example.m
%   transfer_learning_example.m
%   slr_learn.m
%   slr_learn_transfer.m
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: gslr_transfer.m,v $
%

    properties
        
        model
        options
        diagnostics
      
    end

    methods
       function obj = gslr_transfer(varargin)
                            
          obj.options = [];
          
          % parse options
          for i=1:2:length(varargin)
                obj.options.(varargin{i}) = varargin{i+1};
          end

        
       end
       function obj = train(obj,data,design)
           % simply stores input data and design
           
           % transfer learning
           cdata = cell(1,length(data));
           for c=1:length(data)
             cdata{c} = [design{c}.collapse() data{c}.collapse()];
           end

           [obj.model,obj.diagnostics] = slr_learn_transfer(obj.options,cdata);
              
       end
       
       function post = test(obj,data)       
                          
         post = cell(1,length(data));
         for j=1:length(data)
           post{j} = dataset(slr_classify([data{j}.collapse() ones(data{j}.nsamples,1)], obj.model{j}));
         end
         
       end
       
       function [m,desc] = getmodel(obj)
         % return the parameters wrt a class label in some shape 
         % for multiple tasks, the rows represent the models and
         % the columns represent datasets
        
           % return model for all classes           
           m = cell(1,length(obj.model));
           for c=1:length(obj.model)
             m{c} = full(obj.model{c}(:,1:(end-1))); % ignore bias term
           end          
           
           desc = {'unknown'};
                  
       end

    end
end 
