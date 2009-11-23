classdef lr < classifier
%LR unregularized logistic regression
%
% Tries to perform logistic regression using any of the below methods.
%
% SEE ALSO
%   minFunc
%   gpml-matlab
%   mnrfit
%   regularize_lr
%
%   Copyright (c) 2008, Marcel van Gerven, Ali Bahramisharif
%
%   $Log: lr.m,v $
%

    properties

        model; % the weight vector
        
        options; % only used by minFunc       
    end

    methods
       function obj = lr(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       function obj = train(obj,data,design)
            % simply stores input data and design
           
            if iscell(data), error('LR does not take multiple datasets as input'); end

            % remove unlabeled data
            data = data(~isnan(design),:);
            design = design(~isnan(design));
            
            if isnan(obj.nclasses), obj.nclasses = max(design(:,1)); end

            if exist('minFunc','dir') % external code: minFunc

                nexamples = size(data,1);
                nfeatures = size(data,2)+1;

                classidxs = (1:nexamples)' + (design(:,1) - 1) .* nexamples;
                ptargets = zeros(nexamples,obj.nclasses);
                ptargets(classidxs) = 1;
                targets = (1:nexamples)' + (design(:,1) - 1) * nexamples;

                w = zeros(obj.nclasses,nfeatures);
                
                if isempty(obj.options)
                  obj.options.Method='lbfgs';
                  obj.options.display = 0;
                  
                  % parameters for exact minimization (overfits)
%                   obj.options.TolFun = 0;
%                   obj.options.TolX = 0;
%                   obj.options.MaxIter = 1e4;
%                   obj.options.MaxFunEvals = 1e4;
                end
                
                w = minFunc(@(w)logreg(w(:),[data ones(size(data,1),1)],targets,ptargets,obj.nclasses),w(:),obj.options);

                obj.model = reshape(w,obj.nclasses,nfeatures);

            elseif exist('gpml-matlab','dir') % external code: gpml-matlab

                nexamples = size(data,1);
                nfeatures = size(data,2)+1;

                classidxs = (1:nexamples)' + (design(:,1) - 1) .* nexamples;
                ptargets = zeros(nexamples,obj.nclasses);
                ptargets(classidxs) = 1;
                targets = (1:nexamples)' + (design(:,1) - 1) * nexamples;

                w = zeros(obj.nclasses,nfeatures);
                w = minimize(w(:),'logreg',50,[data ones(size(data,1),1)],targets,ptargets,obj.nclasses);

                obj.model = reshape(w,obj.nclasses,nfeatures);

            elseif license('test','statistics_toolbox') && size(data,2) <= 300 && obj.nclasses == 2

              % try matlab native code
              B = mnrfit(data,design(:,1));

              % construct model
              obj.model = [ transpose([B(2:end,:); B(1,:)]); zeros(1,size(B,1))];

            else 
              error('no suitable toolboxes found for LR');
            end
                     
       end

       function post = test(obj,data)
           
           if iscell(data)
               post = cell(1,length(data));
               for j=1:length(data)
                   post{j} = slr_classify([data{j} ones(size(data{j},1),1)], obj.model{j});
               end
           else       
               post = slr_classify([data ones(size(data,1),1)], obj.model);
           end
       end              
       
       function m = getmodel(obj,label,dims)
         % return the parameters wrt a class label in some shape 
         
         if nargin < 2 || isempty(label) || isnan(label)        

           % return model for all classes           
           m = obj.model(:,1:(end-1)); % ignore bias term

         else

           m = obj.model(label,1:(end-1)); % ignore bias term
         
         end
         
         if nargin == 3 && numel(m) == prod(dims)
           m = reshape(m,dims);
         end
         
       end
       
    end
end 
