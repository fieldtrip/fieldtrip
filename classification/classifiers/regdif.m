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
    
    function p = estimate(obj,data,design)
      % simply stores input data and design
      
      p.nclasses = design.nunique;
      
      data = data.X;
      design = design.X;
      
      % proprietary code: minFunc
      if exist('minFunc','dir')
        
        nexamples = size(data,1);
        nfeatures = size(data,2)+1;
        
        classidxs = (1:nexamples)' + (design(:,1) - 1) .* nexamples;
        ptargets = zeros(nexamples,p.nclasses);
        ptargets(classidxs) = 1;
        targets = (1:nexamples)' + (design(:,1) - 1) * nexamples;
        
        options.Method='lbfgs';
        options.TolFun=1e-16;
        options.TolX=1e-16;
        options.Display=obj.verbose;
        options.MaxIter=obj.maxiter;
        options.MaxFunEvals=15000;
        
        w =zeros(p.nclasses,nfeatures);
        
        if strcmp(obj.method,'l2')
          if ~isnan(obj.model_p)
            w=obj.model_p;
          end
          [w,obj.fx] = minFunc(@(w)logregr(w(:),[data ones(size(data,1),1)],targets,ptargets,p.nclasses,obj.lambda),w(:),options);
        elseif strcmp(obj.method,'rov')
          %calculating the bartlett's weight as the start point
          if isnan(obj.model_p)
            w=bartletweight(data,design,w,obj.divnum,obj.nclasses);
          else
            w=obj.model_p;
          end
          [w,obj.fx] = minFunc(@(w)regularizedlogreg(w(:),[data ones(size(data,1),1)],targets,ptargets,p.nclasses,obj.lambda,obj.divnum),w(:),options);
        else
          if ~isnan(obj.model_p)
            w=obj.model_p;
          end
          [w,obj.fx] = minFunc(@(w)logreg(w(:),[data ones(size(data,1),1)],targets,ptargets,p.nclasses),w(:),options);
        end
        p.model = reshape(w,p.nclasses,nfeatures);
        
      else
        error('Please install the minFunc toolbox')
      end
      
    end
    
    function post = map(obj,data)
      
      data = data.X;
      
      if iscell(data)
        post = cell(1,length(data));
        for j=1:length(data)
          post{j} = slr_classify([data{j} ones(size(data{j},1),1)], obj.params.model{j});
        end
      else
        post = slr_classify([data ones(size(data,1),1)], obj.params.model);
      end
      
      post = dataset(post);
      
    end
    
  end
end