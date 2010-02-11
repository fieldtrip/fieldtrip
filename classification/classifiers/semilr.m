classdef semilr < classifier
% Semi-supervised logistic regression, 
% forces the posteriors for predictions to be the same 
% by the mean of 'nu' coefficient. The other parameter 'lambda', is for 
% L2 regularization over logistic regressions.
%
%   Copyright (c) 2009, Ali Bahramisharif, Marcel van Gerven


    properties
       
      partition;  % cell-array which partitions the features
      nu=0;       % coefficient for regularizer over posteriors of unlabeled data
      lambda=0;   % coefficient for L2 regularizer over weights of Logistic Regression
      model;      % the weight vector
      w;          % raw representation of weight vector
      
      nclasses;
      
    end
    
    methods
      
      function obj = semilr(varargin)
        
        obj = obj@classifier(varargin{:});
        
        assert(~isempty(obj.partition));
        
      end
      
      function p = estimate(obj,X,Y)
        
        assert(any(exist('minFunc','dir'))) % external code: minFunc
        
        p.nclasses = obj.nunique(Y);
        
        nparts = length(obj.partition);
        
        dat = X;
        X = cell(1,nparts);
        for c=1:nparts
          X{c} = dat(:,obj.partition{c});
        end
        
        unlabel=cell(1,nparts);
        for i=1:nparts
          unlabel{i}=X{i}(isnan(Y),:);
          unlabel{i}=[unlabel{i} ones(size(unlabel{i},1),1)];
          X{i}=X{i}(~isnan(Y),:);
          X{i}=[X{i} ones(size(X{i},1),1)];
        end
        
        nexamples = size(X{1},1);
        
        % this second is necessary for Y.
        
        Y=Y(~isnan(Y),:);
        
        classidxs = (1:nexamples)' + (Y(:,1) - 1) .* nexamples;
        ptargets = zeros(nexamples,p.nclasses);
        ptargets(classidxs) = 1;
        targets = (1:nexamples)' + (Y(:,1) - 1) * nexamples;
        
        % numer of weight elements and their indices per partition
        psz = cellfun(@length,obj.partition)+1;
        pidx = cell(1,nparts);
        for i=1:nparts
          pidx{i}=(1+(p.nclasses*sum(psz(1:(i-1))))):p.nclasses*sum(psz(1:i));
        end
        
        if isempty(obj.w)
          p.w = zeros(p.nclasses,sum(psz));
        else
          p.w = obj.w;
        end
        
        options.Display='off';
        options.Method='lbfgs';
        %options.MaxIter=5000;
        %options.MaxFunEvals=10000;
        %options.TolFun=1e-10;
        %options.TolX=1e-15;
        
        p.w = minFunc(@(x)semilogreg(x(:),X,targets,ptargets,pidx,psz,p.nclasses,unlabel,obj.nu,obj.lambda),p.w(:),options);
        
        weight=cell(1,nparts);
        for i=1:nparts
          weight{i}=p.w(pidx{i});
          weight{i} = reshape(weight{i},[p.nclasses psz(i)]);
        end
        
        p.model = weight;
        
      end
      
      function Y = map(obj,X)
        
        Y=0;
        for i=1:length(obj.partition)
          Y = Y+slr_classify([X(:,obj.partition{i}) ones(size(X,1),1)], obj.params.model{i});
        end
        Y = bsxfun(@rdivide,Y,sum(Y,2));
        
      end
      
    end
    
end