classdef gnb < classifier
%GNB generalized naive Bayes classifier
%
%   can use various conditional distributions
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: gnb.m,v $
%

    properties

        priors
        conditional = 'normal'; % existing pdf or function handle
        mle; % can be used to specify custom mle function
        params
        
    end

    methods
       function obj = gnb(varargin)
                  
           obj = obj@classifier(varargin{:});                      
       end
       
       function obj = train(obj,data,design)
            
         [data,design] = obj.check_input(data,design);
         
         if isnan(obj.nclasses), obj.nclasses = max(design(:,1)); end
         
         nfeatures = size(data,2);
         
         % estimate class priors
         obj.priors = zeros(obj.nclasses,1);
         for j=1:obj.nclasses
           obj.priors(j) = sum(design(:,1)==j)/size(design,1);
         end
         
         % parameters per class/feature pair
         obj.params = cell(obj.nclasses,nfeatures);
         
         if ~isa(obj.conditional,'function_handle')
           
           % estimate class-conditional means
           for j=1:nfeatures
             for k=1:obj.nclasses
               obj.params{k,j} = mle(data(design(:,1) == k,j),'distribution',lower(obj.conditional));
             end
           end
           
         else
           if isempty(obj.mle)
             
             for j=1:nfeatures
               for k=1:obj.nclasses
                 
                 obj.params{k,j} = mle(data(design(:,1) == k,j),'pdf',obj.conditional);
               end
             end
           else
             for j=1:nfeatures
               for k=1:obj.nclasses
                 
                 obj.params{k,j} = obj.mle(data(design(:,1) == k,j));
               end
             end
           end
         end
         
       end
       function post = test(obj,data)

         data = obj.check_input(data);
         
         post = nan(size(data,1),obj.nclasses);
         
         nparams = length(obj.params{1,1});
         
         for m=1:size(post,1) % iterate over examples
           
           for c=1:obj.nclasses
             
             % compute conditional
             conditional = zeros(1,size(data,2));
             
             if ~isa(obj.conditional,'function_handle')
               
               for j=1:size(data,2)
                 if nparams == 1
                   conditional(j) = pdf(obj.conditional,data(m,j),obj.params{c,j}(1));
                 elseif nparams == 2
                   conditional(j) = pdf(obj.conditional,data(m,j),obj.params{c,j}(1),obj.params{c,j}(2));
                 else
                   conditional(j) = pdf(obj.conditional,data(m,j),obj.params{c,j}(1),obj.params{c,j}(2),obj.params{c,j}(3));
                 end
               end
             else
               for j=1:size(data,2)
                 if nparams == 1
                   conditional(j) = obj.conditional(data(m,j),obj.params{c,j}(1));
                 elseif nparams == 2
                   conditional(j) = obj.conditional(data(m,j),obj.params{c,j}(1),obj.params{c,j}(2));
                 else
                   conditional(j) = obj.conditional(data(m,j),obj.params{c,j}(1),obj.params{c,j}(2),obj.params{c,j}(3));
                 end
               end
             end
             
             % degenerate cases
             if ~obj.priors(c) || any(isinf(conditional)) || ~all(conditional)
               post(m,c) = 0;
               break
             end
             
             % compute probability
             post(m,c) = log(obj.priors(c)) + nansum(log(conditional));
             
           end
           
           % compute normalizing term using log-sum-exp trick
           
           mx = max(post(m,:));
           
           nt = 0;
           for c=1:obj.nclasses
             nt = nt + exp(post(m,c) - mx);
           end
           nt = log(nt) + mx;
           
           % normalize
           post(m,:) = exp(post(m,:) - nt);
           
         end
         
         
       end
       
    end
end
