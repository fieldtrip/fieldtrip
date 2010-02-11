classdef rlda < classifier
%rlda regularized linear discriminant analysis
%
%   Copyright (c) 2008, Pawel Herman


    properties

      kernel_type = 'linear';
      kernel_parameter = 1; % kernel parameter
      C = 1; % regularization parameter
      
      % internal parameters
      
      alpha = 0;
      bias = 0;
      SV;
      
      nclasses;
      
    end
    
    methods
      function obj = rlda(varargin)
        
        obj = obj@classifier(varargin{:});
        
      end
      
      function p = estimate(obj,X,Y)
        
        p.nclasses = obj.nunique(Y);
        
        if p.nclasses~=2, error ('only valid for two-class problems'); end
              
        Y(Y == 1) = -1;
        Y(Y == 2) = 1;
        
        lambda = obj.C;
        p.SV = X;
        
        p.alpha = 0;
        p.bias = 0;
        
        if strcmp(obj.kernel_type,'rbf')
          K = rbf_ker(X,obj.kernel_parameter);
        elseif strcmp(obj.kernel_type,'linear')
          K = X * X';
        else
          disp('Unknown kernel type');
          return
        end
        
        ell = size(K,1);
        ellplus = (sum(Y) + ell)/2;
        yplus = 0.5*(Y + 1);
        ellminus = ell - ellplus;
        yminus = yplus - Y;
        rescale = ones(ell,1)+Y*((ellminus-ellplus)/ell);
        plusfactor = 2*ellminus/(ell*ellplus);
        minusfactor = 2*ellplus/(ell*ellminus);
        B = diag(rescale) - (plusfactor * yplus) * yplus' - (minusfactor * yminus) * yminus';
        p.alpha = (B*K + lambda*eye(ell,ell))\Y;
        p.bias = 0.25*(p.alpha'*K*rescale)/(ellplus*ellminus);
        
        
      end
      function Y = map(obj,X)
        
        if strcmp(obj.kernel_type,'rbf')
          Ktest = rbf_prim(obj.params.SV,X,obj.kernel_parameter);
        elseif strcmp(obj.kernel_type,'linear')
          Ktest = obj.params.SV * X';
        else
          disp('Unknown kernel type');
          post = [];
          return
        end
        
        sgns = sign(Ktest'*obj.params.alpha - obj.params.bias);
        
        Y = zeros(size(sgns,1),obj.params.nclasses);
        Y(sgns == -1,1) = 1;
        Y(sgns == 1,2) = 1;
        
      end
      
    end
end
