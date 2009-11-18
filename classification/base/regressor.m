classdef regressor < predictor
%REGRESSOR abstract regression method class
%
% A regressor takes a variable number of arguments upon construction. 
% During operation, the regressor takes data and
% produces predicted responses as an N x 1 matrix for N examples.
% 
% Subclasses should implement the train and test functions.
%
% OPTIONS
%   'nclasses'  : always nan for a regressor
%   'nfeatures' : number of features (determined from data)
%   'nexamples' : number of examples (determined from data)
%
% SEE ALSO
%   circreg.m
%   linreg.m
%   mulreg.m
%   blinreg.m
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: regressor.m,v $
%        
    
    methods
        
        function obj = regressor(varargin)
               
          % parse options
          for i=1:2:length(varargin)
            if ismember(varargin{i},fieldnames(obj))
              obj.(varargin{i}) = varargin{i+1};
            end
          end

        end
        
        function reg = predict(obj,data)
           % convert posteriors into predictions (mean value)
           
           reg = obj.test(data);
           reg = reg(:,1);
           
        end
    end
    
end 
