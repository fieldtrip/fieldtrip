classdef cca_recon < reconstructor
% CCA_RECON canonical correlation analysis reconstructor class
%
% both data and design will be whitened and reduced in dimensions such that
% they are not rank deficient; dimensions can be further reduced through
% indim and outdim parameters
%
% EXAMPLES:
%
% load 69data; X = response; Y = stimuli;
% p = mva({cca_recon('indim',50,'outdim',50,'verbose',true)});
% p = p.train(X(1:100,:),Y(1:100,:));
% r = p.test(X(101:111,:));
% images(r,1:10,[2 5]);
% figure
% images(Y,101:110,[2 5]);
%
% Copyright (c) 2010, Marcel van Gerven


  properties
  
    A; % the matrix from input to latent variables
    B; % the matrix from output to latent variables
    
    dims; % dimensions of the output data
       
    prepin; % save all transformations to original data
    prepout; % save all transformations to original design
    
    indim;  % reduced number of dimensions for the data
    outdim; % reduced number of dimensions for the design
      
  end

  methods
        
    function obj = cca_recon(varargin)
      
      obj = obj@reconstructor(varargin{:});
      
    end
    
    function p = estimate(obj,X,Y)
      
      if isempty(obj.indim) || obj.indim > size(X,1)-1
        p.indim = min(size(X,1)-1,size(X,2));
        if obj.verbose
          fprintf('reducing data dimensions to %d\n',p.indim);
        end
      else 
        p.indim = obj.indim;
      end
      
      if isempty(obj.outdim) || obj.outdim > size(Y,1)-1
        p.outdim = min(size(Y,1)-1,size(Y,2));
        if obj.verbose
          fprintf('reducing design dimensions to %d\n',p.outdim);
        end
      else
        p.outdim = obj.outdim;
      end
      
      % standardize, pca and whiten the data
      if isempty(obj.prepin)
        p.prepin = mva({standardizer pcanalyzer('proportion',p.indim) whitener});
      else
        p.prepin = obj.prepin;
      end
      p.prepin = p.prepin.train(X);
      
      % standardize, pca and whiten the design
      if isempty(obj.prepout)
        p.prepout = mva({standardizer pcanalyzer('proportion',p.outdim) whitener});
      else
        p.prepout = obj.prepout;
      end
      p.prepout = p.prepout.train(Y);
      
      % perform canonical correlation analysis on preprocessed data
      [p.A,p.B] = canoncorr(p.prepin.test(X),p.prepout.test(Y));
      
    end
    
    function Y = map(obj,X)
      % also returns the estimated sources
      
      % standardize and whiten the input data
      % changes source estimate back to native space
      Y = obj.params.prepin.test(X) * obj.params.A * pinv(obj.params.B);
      
      % invert the whitening and standardization of the design matrix
      Y = obj.params.prepout.untest(Y);
      
      % reshape to old dimensions
      Y = reshape(Y,[size(X,1) obj.outdims]);
            
    end
    
  end
  
end
