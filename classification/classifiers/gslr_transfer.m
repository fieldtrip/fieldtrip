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

    properties
        
      options
        
    end

    methods
       function obj = gslr_transfer(varargin)
                            
          obj.options = [];
          
          % parse options
          for i=1:2:length(varargin)
                obj.options.(varargin{i}) = varargin{i+1};
          end

        
       end
       
       function p = estimate(obj,X,Y)
           % simply stores input data and design
           
           % transfer learning
           cdata = cell(1,length(X));
           for c=1:length(X)
             cdata{c} = [Y{c} X{c}];
           end

           [p.model,p.diagnostics] = slr_learn_transfer(obj.options,cdata);
              
       end
       
       function Y = map(obj,X)       
                          
         Y = cell(1,length(X));
         for j=1:length(X)
           Y{j} = slr_classify([X{j} ones(size(X{j},1),1)], obj.params.model{j});
         end
         
       end
       
       function [m,desc] = getmodel(obj)
         % return the parameters wrt a class label in some shape 
         % for multiple tasks, the rows represent the models and
         % the columns represent datasets
        
           % return model for all classes           
           m = cell(size(obj.params.model{1},1),length(obj.params.model));
           for c=1:size(m,1)
             for j=1:size(m,2)
               m{c,j} = full(obj.params.model{j}(c,1:(end-1))); % ignore bias term
             end
           end          
           
           desc = {'unknown'};
                  
       end

    end
end 
