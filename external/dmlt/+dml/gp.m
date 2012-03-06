classdef gp < dml.method
% GP Gaussian process.
%
%   DESCRIPTION
%
%   REFERENCE
%
%   EXAMPLE
%
%   DEVELOPER
%   Perry Groot


  properties
    
   gproc % gaussian process object
    
  end
  
  methods
    
    function obj = gp(varargin)

      obj = obj@dml.method(varargin{:});

    end
    
    function obj = train(obj,X,Y)
        
      % GP model definition
      gpoptions.kern = 'rbf';
      gpoptions.lik = 'gauss';
      gpoptions.mean = 'zero';
      gpoptions.approx = 'exact';
      gp = gpCreate(X,Y,gpoptions);
      gp = gpSetHyper(gp, 1, 2);         % bit of a hack: set rbf length-scale to 1
      
       % GP Training
      obj.gproc = gpOptimize(gp);
      
      
    end
    
    function Y = test(obj,X)

      % GP prediction
      [mu_1,s2_1] = gpPredict(obj.gproc, X);
      
    end

    function m = model(obj)
      % returns
      %
      
      
      
    end

  end
  
end
