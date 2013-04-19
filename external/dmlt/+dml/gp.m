classdef gp < dml.method
% GP Gaussian process.
%
%   DESCRIPTION
%
%   Wrapper to the GPstuff Gaussian process code
%
%   REFERENCE
%    Jarno Vanhatalo, Jaakko Riihimäki, Jouni Hartikainen, Pasi Jylänki,
%    Ville Tolvanen, Aki Vehtari (2013). GPstuff: A Toolbox for Bayesian
%    Modeling with Gaussian Processes. In Journal of Machine Learning
%    Research, accepted for publication.
%
%   EXAMPLE
%
%   X = randn(20,2);
%   Y = X(:,1);
%   m = dml.gp;
%   m = m.train(X,Y);
%   Z = m.test(X);
%
%   DEVELOPER
%   Jarno Vanhatalo, Jaakko Riihimäki, Jouni Hartikainen, Pasi Jylänki, Ville Tolvanen, Aki Vehtari


  properties
    
   gproc % gaussian process object
    
  end
  
  methods
    
    function obj = gp(varargin)

      obj = obj@dml.method(varargin{:});
      
      lik = lik_gaussian;
      gpcf = gpcf_sexp;
      pn=prior_logunif();
      lik = lik_gaussian(lik, 'sigma2_prior', pn);
      pl = prior_unif();
      pm = prior_sqrtunif();
      gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
      obj.gproc = gp_set('lik', lik, 'cf', gpcf);
      
    end
    
    function obj = train(obj,X,Y)
        
      opt=optimset('TolFun',1e-3,'TolX',1e-3,'Display','iter');
      obj.gproc=gp_optim(obj.gproc,X,Y,'opt',opt);
      
      obj.gproc.Xtrain = X;
      obj.gproc.Ytrain = Y;

    end
    
    function Y = test(obj,X)

       [xt1,xt2]=meshgrid(-1.8:0.1:1.8,-1.8:0.1:1.8);
       xt=[xt1(:) xt2(:)];
       [Y, Varft_map] = gp_pred(obj.gproc, obj.gproc.Xtrain, obj.gproc.Ytrain, X);

    end

  end
  
end
