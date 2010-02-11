classdef preprocessor < mvmethod
%PREPROCESSOR preprocessor method class
%
% A preprocessor is a handle class that takes a variable number of arguments 
% upon construction. During operation, the preprocessor takes data and transforms
% this data to a new representation. E.g., it can construct novel features or
% performs a principal components analysis, etc. 
% 
% Subclasses should implement the estimate, map
% and (possibly) unmap functions. These learn the parameters, 
% apply the transform and (possibly) the inverse transform, respectively.
%
% SEE ALSO
% doc preprocessors
%
% Copyright (c) 2008, Marcel van Gerven

  methods
   
    function obj = preprocessor(varargin)
      
      % parse options
      for i=1:2:length(varargin)
        if ismember(varargin{i},fieldnames(obj))
          obj.(varargin{i}) = varargin{i+1};
        end
      end
      
    end
    
  end
  
end
