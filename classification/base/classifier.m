classdef classifier < predictor
%CLASSIFIER abstract classifier method class
%
% A classifier takes a variable number of arguments upon construction. 
% During operation, the classifier takes data and
% produces posterior probabilities of class labels as an N x C matrix for N
% examples and C classes.
% 
% Subclasses should implement the train and test functions and optionally
% the getmodel function which reshapes parameters to something
% interpretable
%
% SEE ALSO:
% doc classifiers
%
% Copyright (c) 2008, Marcel van Gerven
%
% $Log: classifier.m,v $
%    

  methods
        
        function obj = classifier(varargin)
     
          % parse options
          for i=1:2:length(varargin)
            if ismember(varargin{i},fieldnames(obj))
              obj.(varargin{i}) = varargin{i+1};
            end
          end

        end        
        
        function clf = predict(obj,data)
           % convert posterior dataset into classifications
           
           [tmp,clf] = max(obj.test(data).X,[],2);
        end
        
    end
    
end
