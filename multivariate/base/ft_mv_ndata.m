classdef ft_mv_ndata < ft_mv_method
%FT_MV_NDATA wrapper class to make ft_mv_methods handle multiple datasets
%i.e., when iscell(X) || iscell(Y) is true.
%
% Copyright (c) 2010, Marcel van Gerven
    
  properties
    
    mvmethod  % the method which should be replicated 
    
  end

  methods

    function obj = ft_mv_ndata(varargin)

     obj@ft_mv_method(varargin{:});
     
     obj.verbose = obj.mvmethod.verbose;
    
    end

    function obj = train(obj,X,Y)
      
      if ~iscell(X) && ~iscell(Y), error('method expects a cell-array'); end
      
      if iscell(X) && ~iscell(Y)
        Y = repmat({Y},size(X));
      elseif ~iscell(X) && iscell(Y)
        X = repmat({X},size(Y));
      end
      
      if length(X) ~= length(Y), error('X and Y should have the same number of cell-array elements'); end
      
      nx = length(X);
      
      obj.mvmethod = repmat({obj.mvmethod},nx);
      
      for c=1:nx
        
        if obj.verbose
          fprintf('training dataset %d of %d\n',c,nx);
        end
        
        obj.mvmethod{c} = obj.mvmethod{c}.train(X{c},Y{c});
      end
      
    end
    
    function Y = test(obj,X)
      
      if ~iscell(X), error('method expects a cell-array'); end
      
      nx = length(obj.mvmethod);
      
      Y = cell(size(X));
      for c=1:nx

        if obj.verbose
          fprintf('testing dataset %d of %d\n',c,nx);
        end

        Y{c} = obj.mvmethod{c}.test(X{c});

      end
      
    end
    
    function Y = predict(obj,X)
      
      if ~iscell(X), error('method expects a cell-array'); end
   
      nx = length(obj.mvmethod);
      
      Y = cell(size(X));
      for c=1:nx

        if obj.verbose
          fprintf('predicting dataset %d of %d\n',c,nx);
        end

        Y{c} = obj.mvmethod{c}.predict(X{c});

      end
      
    end
    
    function [m,d] = model(obj)
      
      nx = length(obj.mvmethod);
      
      m = cell(size(X));
      for c=1:nx
        
        [m{c},d] = obj.mvmethod{c}.model;
        
      end
      
    end
    
  end

end