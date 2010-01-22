classdef regdif < classifier
% REGDIF regularization of differences
%
%   Copyright (c) 2008, Ali Bahramisharif, Marcel van Gerven
%
%   $Log: regdif.m,v $
%

  properties
  
    model; % the weight vector
    method=''; % '' (standard logistic regression), 'l2', 'rov' (regularization of variation)
    lambda=0; % regularization parameter
    divnum;   % the number of divisions required for ROV and ROVN
    model_p=NaN; % start model default is using Welch's weights
    fx=1;
    maxiter=500;
    nclasses;
    
  end
  
  methods
    function obj = regdif(varargin)
      
      obj = obj@classifier(varargin{:});
      
    end
    function obj = train(obj,data,design)
      % simply stores input data and design
      
      obj.nclasses = design.nunique;
      
      data = data.collapse();
      design = design.collapse();
      
      % proprietary code: minFunc
      if exist('minFunc','dir')
        
        nexamples = size(data,1);
        nfeatures = size(data,2)+1;
        
        classidxs = (1:nexamples)' + (design(:,1) - 1) .* nexamples;
        ptargets = zeros(nexamples,obj.nclasses);
        ptargets(classidxs) = 1;
        targets = (1:nexamples)' + (design(:,1) - 1) * nexamples;
        
        options.Method='lbfgs';
        options.TolFun=1e-16;
        options.TolX=1e-16;
        options.Display=obj.verbose;
        options.MaxIter=obj.maxiter;
        options.MaxFunEvals=15000;
        
        w =zeros(obj.nclasses,nfeatures);
        
        if strcmp(obj.method,'l2')
          if ~isnan(obj.model_p)
            w=obj.model_p;
          end
          [w,obj.fx] = minFunc(@(w)logregr(w(:),[data ones(size(data,1),1)],targets,ptargets,obj.nclasses,obj.lambda),w(:),options);
        elseif strcmp(obj.method,'rov')
          %calculating the bartlett's weight as the start point
          if isnan(obj.model_p)
            w=bartletweight(data,design,w,obj.divnum,obj.nclasses);
          else
            w=obj.model_p;
          end
          [w,obj.fx] = minFunc(@(w)regularizedlogreg(w(:),[data ones(size(data,1),1)],targets,ptargets,obj.nclasses,obj.lambda,obj.divnum),w(:),options);
        else
          if ~isnan(obj.model_p)
            w=obj.model_p;
          end
          [w,obj.fx] = minFunc(@(w)logreg(w(:),[data ones(size(data,1),1)],targets,ptargets,obj.nclasses),w(:),options);
        end
        obj.model = reshape(w,obj.nclasses,nfeatures);
        
      else
        error('Please install the minFunc toolbox')
      end
      
    end
    
    function post = test(obj,data)
      
      data = data.collapse();
      
      if iscell(data)
        post = cell(1,length(data));
        for j=1:length(data)
          post{j} = slr_classify([data{j} ones(size(data{j},1),1)], obj.model{j});
        end
      else
        post = slr_classify([data ones(size(data,1),1)], obj.model);
      end
      
      post = dataset(post);
      
    end
    
  end
end