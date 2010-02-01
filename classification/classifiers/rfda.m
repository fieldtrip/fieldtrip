classdef rfda < classifier
%RFDA regularized fisher discriminant analysis
%
%   Copyright (c) 2008, Pawel Herman
%
%   $Log: rfda.m,v $
%

    properties

      kernel_type = 'linear';
      kernel_parameter = 1; % kernel parameter
      C = 1; % regularization parameter
      
    end
    
    methods
      
      function obj = rfda(varargin)
        
        obj = obj@classifier(varargin{:});
        
      end
      
      function p = estimate(obj,data,design)
        
        p.nclasses = design.nunique;
        
        if design.nunique~=2, error ('only valid for two-class problems'); end
        
        data = data.X;
        design = design.X;
        
        design(design == 1) = -1;
        design(design == 2) = 1;
        
        lambda = obj.C;
        p.SV = data;
        
        p.alpha = 0;
        p.bias = 0;
        
        if strcmp(obj.kernel_type,'rbf')
          K = rbf_ker(data,obj.kernel_parameter);
        elseif strcmp(obj.kernel_type,'linear')
          K = data * data';
        else
          disp('Unknown kernel type');
          return
        end
        
        ell = size(K,1);
        ellplus = (sum(design) + ell)/2;
        yplus = 0.5*(design + 1);
        ellminus = ell - ellplus;
        yminus = yplus - design;
        rescale = ones(ell,1)+design*((ellminus-ellplus)/ell);
        plusfactor = 2*ellminus/(ell*ellplus);
        minusfactor = 2*ellplus/(ell*ellminus);
        B = diag(rescale) - (plusfactor * yplus) * yplus' - (minusfactor * yminus) * yminus';
        p.alpha = (B*K + lambda*eye(ell,ell))\design;
        p.bias = 0.25*(p.alpha'*K*rescale)/(ellplus*ellminus);
        
        
      end
      
      function post = map(obj,data)
        
        data = data.X;
        
        if strcmp(obj.kernel_type,'rbf')
          Ktest = rbf_prim(obj.params.SV,data,obj.kernel_parameter);
        elseif strcmp(obj.kernel_type,'linear')
          Ktest = obj.params.SV * data';
        else
          disp('Unknown kernel type');
          post = [];
          return
        end
        
        sgns = sign(Ktest'*obj.params.alpha - obj.params.bias);
        
        post = zeros(size(sgns,1),obj.params.nclasses);
        post(sgns == -1,1) = 1;
        post(sgns == 1,2) = 1;
        
        post = dataset(post);
        
      end
      
    end
end
