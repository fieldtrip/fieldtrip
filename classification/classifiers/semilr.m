classdef semilr < classifier
% Semi-supervised logistic regression, 
% forces the posteriors for predictions to be the same 
% by the mean of 'nu' coefficient. The other parameter 'lambda', is for 
% L2 regularization over logistic regressions.
%
%
%   Copyright (c) 2009, Ali Bahramisharif, Marcel van Gerven
%
%   $Log: da.m,v $

    properties
        nu=0;       % coefficient for regularizer over posteriors of unlabeled data
        lambda=0;   % coefficient for L2 regularizer over weights of Logistic Regression
        partition;  % cell-array which partitions the features
        model;      % the weight vector
        w;          % raw representation of weight vector
    end

    methods
       function obj = semilr(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       function obj = train(obj,data,design)
         
         assert(any(exist('minFunc','dir'))) % external code: minFunc

         obj.nclasses = max(design(:,1));
         
         nparts = length(obj.partition);
         
         dat = data;
         data = cell(1,nparts);
         for c=1:nparts
           data{c} = dat(:,obj.partition{c});
         end
         
         unlabel=cell(1,nparts);         
         for i=1:nparts
           unlabel{i}=data{i}(isnan(design),:);
           unlabel{i}=[unlabel{i} ones(size(unlabel{i},1),1)];
           data{i}=data{i}(~isnan(design),:);
           data{i}=[data{i} ones(size(data{i},1),1)];
         end
         
         nexamples = size(data{1},1);
          
         % this second is necessary for design.
          
         design=design(~isnan(design),:);
          
         classidxs = (1:nexamples)' + (design(:,1) - 1) .* nexamples;
         ptargets = zeros(nexamples,obj.nclasses);
         ptargets(classidxs) = 1;
         targets = (1:nexamples)' + (design(:,1) - 1) * nexamples;
         
         % numer of weight elements and their indices per partition
         psz = cellfun(@length,obj.partition)+1;
         pidx = cell(1,nparts);
         for i=1:nparts
           pidx{i}=(1+(obj.nclasses*sum(psz(1:(i-1))))):obj.nclasses*sum(psz(1:i));
         end
         
         if isempty(obj.w)
           obj.w = zeros(obj.nclasses,sum(psz));
         end
         
         options.MaxIter=5000;
         options.MaxFunEvals=10000;
         options.Display='off';
         options.Method='lbfgs';
         %options.TolFun=1e-10;
         %options.TolX=1e-15;
         
         obj.w = minFunc(@(x)semilogreg(x(:),data,targets,ptargets,pidx,psz,obj.nclasses,unlabel,obj.nu,obj.lambda),obj.w(:),options);

         weight=cell(1,nparts);
         for i=1:nparts
           weight{i}=obj.w(pidx{i});
           weight{i} = reshape(weight{i},[obj.nclasses psz(i)]);
         end
         
         obj.model = weight;

       end
       
       function post = test(obj,data)

         post=0;
         for i=1:length(obj.partition)
           post = post+slr_classify([data(:,obj.partition{i}) ones(size(data,1),1)], obj.model{i});
         end
         post = bsxfun(@rdivide,post,sum(post,2));
         
       end
       
    end
    
end