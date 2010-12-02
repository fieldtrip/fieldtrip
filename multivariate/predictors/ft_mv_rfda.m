classdef ft_mv_rfda < ft_mv_predictor
%FT_MV_RFDA regularized fisher discriminant analysis
%
%   Copyright (c) 2008, Pawel Herman
%

    properties

      kernel_type = 'linear'; % linear or rbf
      kernel_parameter = 1; % kernel parameter
      C = 1; % regularization parameter
      
      SV
      alpha = 0;
      bias  = 0;
      
    end
    
    methods
      
      function obj = ft_mv_rfda(varargin)
        
        obj = obj@ft_mv_predictor(varargin{:});
        
      end
      
      function obj = train(obj,X,Y)
        
        if max(Y)>2
          error('FT_MV_RFDA only suitable for binary classification');
        end
        
        Y(Y == 1) = -1;
        Y(Y == 2) = 1;
        
        lambda = obj.C;
        obj.SV = X;
        
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
        obj.alpha = (B*K + lambda*eye(ell,ell))\Y;
        obj.bias = 0.25*(obj.alpha'*K*rescale)/(ellplus*ellminus);
        
        
      end
      
      function Y = test(obj,X)
       
        if strcmp(obj.kernel_type,'rbf')
          Ktest = rbf_prim(obj.SV,X,obj.kernel_parameter);
        elseif strcmp(obj.kernel_type,'linear')
          Ktest = obj.SV * X';
        else
          disp('Unknown kernel type');
          Y = [];
          return
        end
        
        sgns = sign(Ktest'*obj.alpha - obj.bias);
        
        Y = zeros(size(sgns,1),2);
        Y(sgns == -1,1) = 1;
        Y(sgns == 1,2) = 1;
                
      end
      
    end
end
