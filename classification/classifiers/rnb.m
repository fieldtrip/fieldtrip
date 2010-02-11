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

    properties
      
      options % cfg parameter format
      diagnostics
      path % models in the regularization path
      accuracies % accuracies of each model in the path
      
      nclasses
      
    end

    methods
      
      function obj = rnb(varargin)
        
        obj.options = struct(varargin{:});
        
      end
      
       function p = estimate(obj,X,Y)
 
         p.nclasses = obj.nunique(Y);
       
         % create hold-out set
         if ~isfield(obj.options,'cvfolds'), obj.options.cvfolds = 0.9; end
         assert(obj.options.cvfolds > 0 && obj.options.cvfolds < 1);
         
         prm = randperm(size(X,1));
         trainidx = prm(1:ceil(obj.options.cvfolds*size(X,1)));
         testidx = prm((1+ceil(obj.options.cvfolds*size(X,1))):end);
         
         % estimate class priors
         p.priors = zeros(p.nclasses,1);
         for j=1:p.nclasses
           p.priors(j) = sum(Y(:,1)==j)/size(Y,1);
         end
         
         [p.path,p.diagnostics] = regularize_nb(obj.options,[Y(trainidx,1) X(trainidx,:)]);
         
         % evaluate performance of all models
         
         nfeatures = size(X,2);
         p.means = zeros(obj.nclasses,nfeatures);
         p.stds = zeros(obj.nclasses,nfeatures);
         p.accuracies = zeros(1,length(p.diagnostics.activeset));
         for i=1:length(p.diagnostics.activeset)
           
           for j=1:nfeatures
             for k=1:p.nclasses
               p.means(k,j) = p.path{i}(k,j);
             end
           end
           
           for j=1:nfeatures
             for k=1:p.nclasses
               p.stds(k,j) = p.path{i}(k+p.nclasses,j);
             end
           end
           
           obj.params = p;
           
           Z = obj.map(X(testidx,:));
           
           p.accuracies(i) = obj.evaluate(Z,Y(testidx,1),'metric','accuracy');
           
         end
         
         % get best model
         
         [a,midx] = max(p.accuracies);
         
         for j=1:nfeatures
           for k=1:p.nclasses
             p.means(k,j) = p.path{midx}(k,j);
           end
         end
         
         for j=1:nfeatures
           for k=1:p.nclasses
             p.stds(k,j) = p.path{midx}(k+p.nclasses,j);
           end
         end
         
       end
       
       function Y = map(obj,X)
         
         Y = zeros(size(X,1),obj.params.nclasses);
         
         for m=1:size(Y,1) % iterate over examples
           
           for c=1:obj.params.nclasses
             
             % compute probability
             Y(m,c) = log(obj.params.priors(c)) + mynansum(log(1./(sqrt(2*pi)*obj.params.stds(c,:)) .* exp(- (X(m,:) - obj.params.means(c,:)).^2./(2*obj.params.stds(c,:).^2))));
             
           end
           
           % compute normalizing term using log-sum-exp trick
           
           mx = max(Y(m,:));
           
           nt = 0;
           for c=1:obj.params.nclasses
             nt = nt + exp(Y(m,c) - mx);
           end
           nt = log(nt) + mx;
           
           % normalize
           Y(m,:) = exp(Y(m,:) - nt);
           
         end
         
       end
       
    end
end
