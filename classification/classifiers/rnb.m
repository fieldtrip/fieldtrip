classdef rnb < nb
%RNB regularized naive Bayes classifier
%
%   EXAMPLE:
%
%   regularization with a prespecified regularization path, tolerance level
%   and step size:
%  
%   rnb('lambas',[100 10 1 0],'tolerance',0,'epsilon',1e-2)
%
%   SEE ALSO:
%   regularize_nb.m
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: rnb.m,v $
%

    properties

        options % cfg parameter format
        diagnostics
        path % models in the regularization path
        accuracies % accuracies of each model in the path
        
    end

    methods
       function obj = rnb(varargin)
                  
           obj.options = struct(varargin{:});
           
       end
       function obj = train(obj,data,design)
 
         obj.nclasses = design.nunique;
         
         data = data.collapse();
         design = design.collapse();
         
         % create hold-out set
         if ~isfield(obj.options,'cvfolds'), obj.options.cvfolds = 0.9; end
         assert(obj.options.cvfolds > 0 && obj.options.cvfolds < 1);
         
         prm = randperm(size(data,1));
         trainidx = prm(1:ceil(obj.options.cvfolds*size(data,1)));
         testidx = prm((1+ceil(obj.options.cvfolds*size(data,1))):end);
         
         % estimate class priors
         obj.priors = zeros(obj.nclasses,1);
         for j=1:obj.nclasses
           obj.priors(j) = sum(design(:,1)==j)/size(design,1);
         end
         
         [obj.path,obj.diagnostics] = regularize_nb(obj.options,[design(trainidx,1) data(trainidx,:)]);
         
         % evaluate performance of all models
         
         nfeatures = size(data,2);
         obj.means = zeros(obj.nclasses,nfeatures);
         obj.stds = zeros(obj.nclasses,nfeatures);
         obj.accuracies = zeros(1,length(obj.diagnostics.activeset));
         for i=1:length(obj.diagnostics.activeset)
           
           for j=1:nfeatures
             for k=1:obj.nclasses
               obj.means(k,j) = obj.path{i}(k,j);
             end
           end
           
           for j=1:nfeatures
             for k=1:obj.nclasses
               obj.stds(k,j) = obj.path{i}(k+obj.nclasses,j);
             end
           end
           
           post = obj.test(dataset(data(testidx,:)));
           
           obj.accuracies(i) = validator.eval(post.X,design(testidx,1),'metric','accuracy');
           
         end
         
         % get best model
         
         [a,midx] = max(obj.accuracies);
         
         for j=1:nfeatures
           for k=1:obj.nclasses
             obj.means(k,j) = obj.path{midx}(k,j);
           end
         end
         
         for j=1:nfeatures
           for k=1:obj.nclasses
             obj.stds(k,j) = obj.path{midx}(k+obj.nclasses,j);
           end
         end
         
       end
       
       function post = test(obj,data)
         
         data = data.collapse();

         post = zeros(size(data,1),obj.nclasses);
         
         for m=1:size(post,1) % iterate over examples
           
           for c=1:obj.nclasses
             
             % compute probability
             post(m,c) = log(obj.priors(c)) + mynansum(log(1./(sqrt(2*pi)*obj.stds(c,:)) .* exp(- (data(m,:) - obj.means(c,:)).^2./(2*obj.stds(c,:).^2))));
             
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
         
         post = dataset(post);
         
       end
       
    end
end
