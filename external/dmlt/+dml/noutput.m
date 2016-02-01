classdef noutput < dml.method
%NOUTPUT wrapper class to make methods handle multiple outputs.
%
%   DESCRIPTION
%   
%
%   REFERENCE
%   
%
%   EXAMPLE
%   
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)

    
  properties
    
    method  % the method which should be replicated 
    
  end

  methods

    function obj = noutput(varargin)

     obj@dml.method(varargin{:});

     if ~iscell(obj.method)
       obj.verbose = obj.method.verbose;
     else
       obj.verbose = obj.method{1}.verbose;
     end
     
    end

    function obj = train(obj,X,Y)
      
      ny = size(Y,2);
      
      if ~iscell(obj.method)
        obj.method = repmat({obj.method},[1 ny]);
      end
      
      for c=1:ny
        
        if obj.method{c}.verbose
          fprintf('training output %d of %d\n',c,ny);
        end
        
        obj.method{c} = obj.method{c}.train(X,Y(:,c));
      end
      
    end
    
    function Y = test(obj,X)
      
      ny = length(obj.method);
      
      Y = [];
      for c=1:ny
        
        if obj.method{c}.verbose
          fprintf('testing output %d of %d\n',c,ny);
        end
        
        Y = [Y obj.method{c}.test(X)];
        
      end
      
    end
    
    function m = model(obj)
      
      ny = length(obj.method);
      
      m = cell(1,ny);
      for c=1:ny
        
        m{c} = obj.method{c}.model;
        
      end
      
    end
    
  end

end