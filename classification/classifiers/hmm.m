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
%
%   $Log: hmm.m,v $
%

    properties       
        mixture = 1; % increase to use mixtures of gaussians as features
        ar = false;  % set to true to allow autoregressive (self) connections
    end

    methods
       function obj = hmm(varargin)
                  
           obj = obj@classifier(varargin{:});
       end
       function obj = train(obj,data,design)
            
           % data must accommodate hidden variables
           if obj.mixture > 1
             
               numvar = 2*(size(data,2)/obj.horizon);
               ncont = (numvar)/2;

               hdata = nan(size(data,1),obj.horizon*numvar);
               
               for j=1:obj.horizon
                   hdata(:,((j-1)*numvar+1):((j-1)*numvar+ncont)) = data(:,((j-1)*ncont+1):(j*ncont));
               end
            
               obj = obj.train@dynamic_classifier(hdata,design);
               
           else
               obj = obj.train@dynamic_classifier(data,design);
           end
       end
       function post = test(obj,data)
          
           % data must accommodate hidden variables
           if obj.mixture > 1
             
               numvar = 2*(size(data,2)/obj.horizon);
               ncont = (numvar)/2;

               hdata = nan(size(data,1),obj.horizon*numvar);
               
               for j=1:obj.horizon
                   hdata(:,((j-1)*numvar+1):((j-1)*numvar+ncont)) = data(:,((j-1)*ncont+1):(j*ncont));
               end
               
               post = obj.test@dynamic_classifier(hdata);
               
           else
               post = obj.test@dynamic_classifier(data);
           end
           
       end
       function factors = construct_factors(obj)
                 
           % definition of HMM factors

           factors = cell(1,2*obj.numvar);

           % discrete part

           factors{1} = multinomial_cpd(1,[],ones([obj.nclasses 1]));
           factors{obj.numvar+1} = multinomial_cpd(obj.numvar+1,1,ones(obj.nclasses));

           if obj.mixture > 1

               % continuous part
               ncont = ((obj.numvar-1)/2);
               
               % first slice
               
               % continuous variables
               for f=2:(ncont+1)
                   factors{f} = gaussian_cpd(f,[],[1 f+ncont],randn([obj.nclasses obj.mixture]),cell([obj.nclasses obj.mixture]),ones([obj.nclasses obj.mixture]));
               end

               % hidden discrete variables
               for f=(ncont+2):obj.numvar
                   factors{f} = multinomial_cpd(f,[], 1/obj.mixture - 1e-2 + 2e-2*rand([obj.mixture 1]));
               end
               
               % second slice
               
               if obj.ar
                   
                   % add autoregressive links

                   % continuous variables
                   for f=2:(ncont+1)
                       
                       betas = cell([obj.nclasses obj.mixture]);
                       betas(:) = { randn };
                       
                       factors{obj.numvar+f} = gaussian_cpd(obj.numvar + f,f,[obj.numvar+1 obj.numvar+f+ncont], ...
                           randn([obj.nclasses obj.mixture]),betas,ones([obj.nclasses obj.mixture]));
                   end

                   % hidden discrete variables
                   for f=(ncont+2):obj.numvar
                       factors{obj.numvar+f} = multinomial_cpd(obj.numvar+f,[],1/obj.mixture - 1e-2 + 2e-2*rand([obj.mixture 1]));
                   end
               
               else

                   % naive bayes like structure
                   
                   % continuous variables
                   for f=2:(ncont+1)
                       factors{obj.numvar+f} = gaussian_cpd(obj.numvar + f,[],[obj.numvar+1 obj.numvar+f+ncont], ...
                           randn([obj.nclasses obj.mixture]),cell([obj.nclasses obj.mixture]),ones([obj.nclasses obj.mixture]));
                   end

                   % hidden discrete variables
                   for f=(ncont+2):obj.numvar
                       factors{obj.numvar+f} = multinomial_cpd(obj.numvar+f,[],1/obj.mixture - 1e-2 + 2e-2*rand([obj.mixture 1]));
                   end
                   
               end
           
           else % no mixture components
               
               % continuous part

               % first slice
               for f=2:obj.numvar
                   factors{f} = gaussian_cpd(f,[],1,randn([obj.nclasses 1]),cell([obj.nclasses 1]),rand([obj.nclasses 1]));
               end

               % second slice
               if obj.ar

                   % add autoregressive links
                   for f=(obj.numvar+2):(2*obj.numvar)
                       
                       betas = cell([obj.nclasses 1]);
                       betas(:) = { randn };
                       
                       factors{f} = gaussian_cpd(f,f-obj.numvar,obj.numvar+1,randn([obj.nclasses 1]),betas,rand([obj.nclasses 1]));
                   end
               else

                   % naive bayes like structure
                   for f=(obj.numvar+2):(2*obj.numvar)
                       factors{f} = gaussian_cpd(f,[],obj.numvar+1,randn([obj.nclasses 1]),cell([obj.nclasses 1]),rand([obj.nclasses 1]));
                   end
               end
           end
       end                       

    end
end 
