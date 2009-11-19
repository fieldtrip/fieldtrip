classdef preprocessor < clfmethod
%PREPROCESSOR preprocessor method class
%
% A preprocessor is a handle class that takes a variable number of arguments 
% upon construction. During operation, the preprocessor takes data and transforms
% this data to a new representation. E.g., it can construct novel features or
% performs a principal components analysis, etc. 
% 
% Subclasses should implement the train and test functions 
% (and possibly getmodel function) and must be able to act on 
% cell arrays as well as normal data matrices.
%
% PROPERTIES
%   'prefun'  : custom preprocessing function that can be added
%   'verbose' : output comment if true
%
% EXAMPLE
% A custom preprocessor can be constructed on the fly as follows:
% 
%     obj = preprocessor('function',@myfunction);
%
% SEE ALSO
% doc preprocessors
%
% Copyright (c) 2008, Marcel van Gerven
%
% $Log: preprocessor.m,v $
%
    properties         
        prefun;
    end

    methods
        function obj = preprocessor(varargin)
         % fun is a custom function
         
            % parse options 
            for i=1:2:length(varargin)
              if ismember(varargin{i},fieldnames(obj))
                obj.(varargin{i}) = varargin{i+1};
              end
            end
            
        end
        
        function obj = train(obj,data,design)
        end
        
        function data = test(obj,data)
            
            if iscell(data)
                                
               for c=1:length(data)
                   data{c} = obj.test(data{c});
               end
            else
                data = obj.prefun(data);
            end
        end
        
        
        
    end
end
