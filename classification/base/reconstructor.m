classdef reconstructor < predictor
%RECONSTRUCTOR reconstructor class
%
% a reconstructor takes multidimensional inputs and outputs and
% reconstructs the outputs. This can be done by simply applying classifiers
% to individual outputs but also by much more sophisticated methods
% 
% Subclasses should implement the train and test functions and optionally
% the getmodel function which reshapes parameters to something
% interpretable
%
% Copyright (c) 2010, Marcel van Gerven
  

  methods
        
        function obj = reconstructor(varargin)
     
          % parse options
          for i=1:2:length(varargin)
            if ismember(varargin{i},fieldnames(obj))
              obj.(varargin{i}) = varargin{i+1};
            end
          end

        end        
        
        function rec = predict(obj,data)
          % convert outcome dataset into predictions
          
          rec = obj.test(data).X; 
        end
        
    end
    
end
