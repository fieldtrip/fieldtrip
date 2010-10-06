classdef ft_mv_noutput < ft_mv_method
%FT_MV_NOUTPUT wrapper class to make mvmethods handle multiple outputs
%i.e., size(Y,2) > 1
%
% Copyright (c) 2010, Marcel van Gerven
    
  properties
    
    mvmethod  % the method which should be replicated 
    
  end

  methods

    function obj = ft_mv_noutput(varargin)

     obj@ft_mv_method(varargin{:});

     obj.verbose = obj.mvmethod.verbose;
     
    end

    function obj = train(obj,X,Y)
      
      ny = size(Y,2);
      
      obj.mvmethod = repmat({obj.mvmethod},[1 ny]);
      
      for c=1:ny
        
        if obj.verbose
          fprintf('training output %d of %d\n',c,ny);
        end
        
        obj.mvmethod{c} = obj.mvmethod{c}.train(X,Y(:,c));
      end
      
    end
    
    function Y = test(obj,X)
      
      ny = length(obj.mvmethod);
      
      Y = [];
      for c=1:ny
        
        if obj.verbose
          fprintf('testing output %d of %d\n',c,ny);
        end
        
        Y = [Y obj.mvmethod{c}.test(X)];
        
      end
      
    end
    
    function Y = predict(obj,X)
      
      ny = length(obj.mvmethod);
      
      Y = zeros(size(X,1),ny);
      for c=1:ny
        
        if obj.verbose
          fprintf('predicting output %d of %d\n',c,ny);
        end
        
        Y(:,c) = obj.mvmethod{c}.predict(X);
      end
      
    end
    
    function [m,d] = model(obj)
      
      ny = length(obj.mvmethod);
      
      m = cell(1,ny);
      for c=1:ny
        
        [m{c},d] = obj.mvmethod{c}.model;
        
      end
      
    end
    
  end

end