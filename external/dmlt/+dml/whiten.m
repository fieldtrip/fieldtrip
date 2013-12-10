classdef whiten < dml.method
%WHITEN whitens the data.
%
%  DESCRIPTION
%  Whitening of data through a whitening matrix W. Input data will be 
%  standardized automatically.
%
%   EXAMPLE
%   X = rand(10,200);
%   m = dml.whiten('indim',[2 100],'tdim',1);
%   m = m.train(X);
%   Z = m.test(X);
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)
%   Jason Farquhar (j.farquhar@donders.ru.nl)
   
  properties
    
    W     % whitening matrix
    invW  % inverse whitening matrix
    
    tdim   % target dimension

    indim  % input dimensions
    
    std = dml.standardizer % standardizer object
    
  end

  methods
    
    function obj = whiten(varargin)
      
      obj = obj@dml.method(varargin{:});
      
    end
    
    function obj = train(obj,X,Y)
       
      % handle multiple datasets
      if iscell(X)
        obj = dml.ndata('method',obj);
        obj = obj.train(X);
        return;
      end
      
      % standardization
      obj.std = obj.std.train(X);
      X = obj.std.test(X);
      
      if obj.verbose, fprintf('whitening data\n'); end
      
      % N.B. whitening matrix: W = U*diag(D.^order);
      %      and inverse whitening matrix: W^-1 = U*diag(D.^-order);
      [obj.W,D,wX,U] = whiten(reshape(X,[size(X,1) obj.indim]),obj.tdim+1,-0.9,0,0,0,[],1e-6,1,-.5);
      obj.invW = (U * diag(D.^(.5)))';
      
    end
    
    function Y = test(obj,X)
         
      % standardize
      X = obj.std.test(X);
      
      p1 = 1:(numel(obj.indim)+1);
      p1(obj.tdim+1) = -p1(obj.tdim+1);
      p2 = [(obj.tdim+1) -(obj.tdim+1)];

      Y = tprod(reshape(X,[size(X,1) obj.indim]),p1,obj.W,p2);
      Y = Y(:,:);
      
    end
    
    function X = invert(obj,Y)
      % invert mapping
      
      p1 = 1:(numel(obj.indim)+1);
      p1(obj.tdim+1) = -p1(obj.tdim+1);
      p2 = [(obj.tdim+1) -(obj.tdim+1)];

      X = tprod(reshape(Y,[size(Y,1) obj.indim]),p1,obj.invW,p2);
      
      % invert standardize
      X = obj.std.invert(X(:,:));
    
    end
    
  end
  
end