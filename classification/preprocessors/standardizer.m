classdef standardizer < preprocessor
%STANDARDIZER standardizes data to have mean 0 and standard deviation 1
%
% OPTIONS
%   'bmean' : if false does not subtract means
%   'bstd'  : if false does not divide by standard deviations
%
% SEE ALSO
%   zscore.m
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: standardizer.m,v $
%

    properties
       means; % mean of each feature based on training data
       stds; % sd of each feature based on training data
       
       % take mean
       bmean = true;
       
       % take standard deviation; if zero then do not take standard
       % deviation; otherwise we multiply the scaled data by bstd
       bstd = 1; 
    end

    methods
        
      function obj = standardizer(varargin)
        
        obj = obj@preprocessor(varargin{:});
        
      end
      
      function obj = train(obj,data,design)
        
        if iscell(data)
          
          cmeans = cell(1,length(data));
          cstds = cell(1,length(data));
          
          for c=1:length(data)
            obj = obj.train(data{c},design{c});
            cmeans{c} = obj.means;
            cstds{c}  = obj.stds;
          end
          
          obj.means = cmeans;
          obj.stds  = cstds;
          
        else
          
          obj.means = mynanmean(data.X);
          
          obj.stds = mynanstd(data.X);
          obj.stds(obj.stds==0) = 1; % bug fix
        
        end
      end
      
      function out = test(obj,data)
        
        if iscell(data)
          
          cmeans = obj.means;
          cstds = obj.stds;
          
          for c=1:length(data)
            
            obj.means = cmeans{c};
            obj.stds  = cstds{c};
            data{c} = obj.test(data{c});
          end
          
          obj.means = cmeans;
          obj.stds  = cstds;
          
        else
          
          Y = bsxfun(@minus,data.X,obj.means);
          Y = bsxfun(@rdivide,Y,obj.stds);
          
          out = dataset(Y);

          out.dims = data.dims;
          out.ndims = data.ndims;
          
        end
      end
      
    end
end
