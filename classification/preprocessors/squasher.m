classdef squasher < preprocessor
% SQUASHER squashes data between 0 and 1
%
%   Copyright (c) 2009, Marcel van Gerven
%
%   $Log: squasher.m,v $
%

  properties

    dmin; % minimum of training data
    dmax; % maximum of training data
    
  end
  
  methods
    
    function obj = squasher(varargin)
      
      obj = obj@preprocessor(varargin{:});
      
    end
    function obj = train(obj,data,design)
      
      if iscell(data)
        
        for c=1:length(data)
          
          obj.dmax = cell(1,length(data));
          obj.dmin = cell(1,length(data));
          
          for c=1:length(data)
            
            obj.dmin{c} = min(data{c});
            data{c} = bsxfun(@minus,data{c}, obj.dmin{c});
            obj.dmax{c} = max(data{c});
            obj.dmax{c}(obj.dmax{c} ==0) = 1;
            
          end
          
        end
      else
        
        obj.dmin = min(data);
        data = bsxfun(@minus,data, obj.dmin);
        obj.dmax = max(data);
        obj.dmax(obj.dmax ==0) = 1;
        
      end
    end
    
    function data = test(obj,data)
      
      if iscell(data)
        
        for c=1:length(data)
          
          data{c} = bsxfun(@minus,data{c}, obj.dmin);
          data{c} = bsxfun(@rdivide,data{c}, obj.dmax);
          
        end
        
      else
        
        data = bsxfun(@minus,data, obj.dmin);
        data = bsxfun(@rdivide,data, obj.dmax);
        
      end
    end
    
  end
end
