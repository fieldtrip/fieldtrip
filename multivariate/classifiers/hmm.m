classdef hmm < dynamic_classifier
%HMM hidden Markov model
%
%   this model has as extra parameters:
%   obj.mixture specifies the number of gaussian mixtures to use for each
%       feature variable
%   obj.ar if true specifies an autoregressive HMM (self connections only)   
%
%   SEE ALSO:
%   asynchronous_hmm_example.m
%   synchronous_hmm_example.m    
%
%   Copyright (c) 2008, Marcel van Gerven

    properties       
        mixture = 1; % increase to use mixtures of gaussians as features
        ar = false;  % set to true to allow autoregressive (self) connections
    end

    methods
      
      function obj = hmm(varargin)
        
        obj = obj@dynamic_classifier(varargin{:});
      end
      
      function p = estimate(obj,X,Y)
        
        data = X;
        design = Y;
        
        % data must accommodate hidden variables
        if obj.mixture > 1
          
          numvar = 2*(size(data,2)/obj.horizon);
          ncont = (numvar)/2;
          
          hdata = nan(size(data,1),obj.horizon*numvar);
          
          for j=1:obj.horizon
            hdata(:,((j-1)*numvar+1):((j-1)*numvar+ncont)) = data(:,((j-1)*ncont+1):(j*ncont));
          end
          
          p = obj.estimate@dynamic_classifier(hdata,design);
          
        else
          p = obj.estimate@dynamic_classifier(data,design);
        end
      end
      
      function Y = map(obj,X)
        
        % data must accommodate hidden variables
        if obj.mixture > 1
          
          numvar = 2*(size(X,2)/obj.horizon);
          ncont = (numvar)/2;
          
          hdata = nan(size(X,1),obj.horizon*numvar);
          
          for j=1:obj.horizon
            hdata(:,((j-1)*numvar+1):((j-1)*numvar+ncont)) = X(:,((j-1)*ncont+1):(j*ncont));
          end
          
          Y = obj.map@dynamic_classifier(hdata);
          
        else
          Y = obj.map@dynamic_classifier(X);
        end
        
      end
      
      function factors = construct_factors(obj)
        
        numvar = obj.params.numvar;
        nclasses = obj.params.nclasses;
        
        % definition of HMM factors
        
        factors = cell(1,2*numvar);
        
        % discrete part
        
        factors{1} = multinomial_cpd(1,[],ones([nclasses 1]));
        factors{numvar+1} = multinomial_cpd(numvar+1,1,ones(nclasses));
        
        if obj.mixture > 1
          
          % continuous part
          ncont = ((numvar-1)/2);
          
          % first slice
          
          % continuous variables
          for f=2:(ncont+1)
            factors{f} = gaussian_cpd(f,[],[1 f+ncont],randn([nclasses obj.mixture]),cell([nclasses obj.mixture]),ones([nclasses obj.mixture]));
          end
          
          % hidden discrete variables
          for f=(ncont+2):numvar
            factors{f} = multinomial_cpd(f,[], 1/obj.mixture - 1e-2 + 2e-2*rand([obj.mixture 1]));
          end
          
          % second slice
          
          if obj.ar
            
            % add autoregressive links
            
            % continuous variables
            for f=2:(ncont+1)
              
              betas = cell([nclasses obj.mixture]);
              betas(:) = { randn };
              
              factors{numvar+f} = gaussian_cpd(numvar + f,f,[numvar+1 numvar+f+ncont], ...
                randn([nclasses obj.mixture]),betas,ones([nclasses obj.mixture]));
            end
            
            % hidden discrete variables
            for f=(ncont+2):numvar
              factors{numvar+f} = multinomial_cpd(numvar+f,[],1/obj.mixture - 1e-2 + 2e-2*rand([obj.mixture 1]));
            end
            
          else
            
            % naive bayes like structure
            
            % continuous variables
            for f=2:(ncont+1)
              factors{numvar+f} = gaussian_cpd(numvar + f,[],[numvar+1 numvar+f+ncont], ...
                randn([nclasses obj.mixture]),cell([nclasses obj.mixture]),ones([nclasses obj.mixture]));
            end
            
            % hidden discrete variables
            for f=(ncont+2):numvar
              factors{numvar+f} = multinomial_cpd(numvar+f,[],1/obj.mixture - 1e-2 + 2e-2*rand([obj.mixture 1]));
            end
            
          end
          
        else % no mixture components
          
          % continuous part
          
          % first slice
          for f=2:numvar
            factors{f} = gaussian_cpd(f,[],1,randn([nclasses 1]),cell([nclasses 1]),rand([nclasses 1]));
          end
          
          % second slice
          if obj.ar
            
            % add autoregressive links
            for f=(numvar+2):(2*numvar)
              
              betas = cell([nclasses 1]);
              betas(:) = { randn };
              
              factors{f} = gaussian_cpd(f,f-numvar,numvar+1,randn([nclasses 1]),betas,rand([nclasses 1]));
            end
          else
            
            % naive bayes like structure
            for f=(numvar+2):(2*numvar)
              factors{f} = gaussian_cpd(f,[],numvar+1,randn([nclasses 1]),cell([nclasses 1]),rand([nclasses 1]));
            end
          end
        end
      end
      
    end
end
