classdef ft_mv_whitener < ft_mv_preprocessor
%WHITENER whitens the data; input data should be standardized
%
% Copyright (c) 2009, Jason Farquhar, Marcel van Gerven
   
% WE STILL NEED TO ENCODE THE DIM
% ASSUMPTION DIM1 = trials, DIM2 = channels, DIM3 = other crap

  properties
    
    W     % whitening matrix
    invW  % inverse whitening matrix
    
  end

  methods
    
    function obj = ft_mv_whitener(varargin)
      
      obj = obj@ft_mv_preprocessor(varargin{:});
      
    end
    
    function obj = train(obj,X,Y)
       
      if nargin<3, Y = []; end
      
      % multiple datasets
      if iscell(X) || iscell(Y)
        obj = ft_mv_ndata('mvmethod',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      % missing data
      if any(isnan(X(:))) || any(isnan(Y(:))), error('method does not handle missing data'); end
     
      if obj.verbose, fprintf('whitening data\n'); end
      
      % N.B. whitening matrix: W = U*diag(D.^order);
      %      and inverse whitening matrix: W^-1 = U*diag(D.^-order);
      [obj.W,D,wX,U] = whiten(X,2,-0.9,0,0,0,[],1e-6,1,-.5);
      obj.invW = (U * diag(D.^(.5)))';
      
    end
    
    function Y = test(obj,X)
      % whiten
         
      Y = tprod(X,[1 -2 3],obj.W,[-2 2]);
      
    end
    
    function X = invert(obj,Y)
      % invert mapping
      
      X = Y*obj.invW;
      
    end
    
  end
  
end