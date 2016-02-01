classdef ndata < dml.method
% NDATA wrapper class to make methods handle multiple datasets.
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)

    
  properties
    
    method  % the method which should be replicated 
    
  end

  methods

    function obj = ndata(varargin)

     obj@dml.method(varargin{:});
     
     obj.verbose = obj.method.verbose;
    
    end

    function obj = train(obj,X,Y)
      
      if ~iscell(X) && ~iscell(Y), error('method expects a cell-array'); end
      if nargin==3 && (length(X) ~= length(Y)), error('X and Y should have the same number of cell-array elements'); end
      
      obj.method = repmat({obj.method},[length(X) 1]);
      
      for c=1:length(X)
        
        if obj.verbose, fprintf('training dataset %d of %d\n',c,length(X)); end
        
        if nargin==3
          obj.method{c} = obj.method{c}.train(X{c},Y{c});
        else
          obj.method{c} = obj.method{c}.train(X{c});
        end
        
      end
      
    end
    
    function Y = test(obj,X)
      
      if ~iscell(X), error('method expects a cell-array'); end
      
      Y = cell(length(X),1);
      for c=1:length(X)

        if obj.verbose, fprintf('testing dataset %d of %d\n',c,length(X)); end

        Y{c} = obj.method{c}.test(X{c});

      end
      
    end
    
    function m = model(obj)
      
      m = cell(length(obj.method),1);
      for c=1:length(obj.method)
        m{c} = obj.method{c}.model;
      end
      
    end
    
  end

end