classdef kernelmethod < classifier
%KERNELMETHOD kernel method
%
% Kernel method using L2 regularization. The method can be any of
% @l2svm_cg, @rkls, or @klr_cg. Options for this method
% are the type of kernel and kernel parameters, as well as the
% regularization constant C. 
%
% Options:
% 'method' : kernel method used
% 'kernel' : type of kernel
% 'C' : regularization parameter; use optimizer class to find an optimal
%       value
%
% REQUIRES:
% Farquhar toolbox
%
% EXAMPLE:
%
% % definition of the classification procedure 
%
% myproc = clfproc({ kernelmethod('method',@l2svm_cg)});
%
% SEE ALSO:
% l2svm_cg.m
% rkls.m
% klr_cg.m
%
% Copyright (c) 2008, Marcel van Gerven, Jason Farquhar
%
% $Log: kernelmethod.m,v $
%

properties
  method = @l2svm_cg;
  kernel = 'linear';
  p1 = nan;
  p2 = nan;
  C = nan;
  
  weights;
  f;
  J;
  traindata;
  primal;
  
end

methods
  
  function obj = kernelmethod(varargin)
    
    obj = obj@classifier(varargin{:});
    
  end
  
  function obj = train(obj,data,design)
    % simply stores input data and design
  
    if design.nunique ~= 2, error('l2svm only makes binary classifications'); end
    
    data = data.X;
    design = design.X;
    
    % transform elements of the design matrix to class labels
    targets = design(:,1);
    targets(design == 1) = -1;
    targets(design == 2) = 1;
    
    % set input for our classifier; compute kernel matrix
    
    switch lower(obj.kernel)
      
      case {'poly','npoly'};     % polynomial
        
        if isnan(obj.p1), obj.p1 = 2/3; end
        if isnan(obj.p2), obj.p2 = 1; end
        
        K = compKernel(data,data,obj.kernel,obj.p1,obj.p2);
        
      case {'rbf','nrbf'};       % Radial basis function, a.k.a. gaussian
        
        if isnan(obj.p1), obj.p1 = .1*(mean(diag(data))-mean(data(:))); end
        
        K = compKernel(data,data,obj.kernel,obj.p1);
        
      otherwise
        K = compKernel(data,data,obj.kernel);
    end
    
    % regularization parameter
    if isnan(obj.C), obj.C = .1*(mean(diag(K))-mean(K(:))); end
    
    if obj.verbose
      fprintf('regularization parameter was set to %.2g\n',obj.C);
    end
    
    [obj.weights,obj.f,obj.J] = obj.method(K,targets,obj.C,'verb',-1);
    obj.traindata = data;
    
    % weights in primal form
    obj.primal = 0;
    for j=1:size(data,1), obj.primal = obj.primal + obj.weights(j)*data(j,:); end
    
  end
  
  function post = test(obj,data)
    
    data = data.X;
    
    % deal with empty data
    if size(data,2) == 0
      
      % random assignment
      post = rand([size(data,1) obj.nclasses]);
      post = double(post == repmat(max(post,[],2),[1 obj.nclasses]));
      
      return
    end
    
    switch lower(obj.kernel)
      
      case {'poly','npoly'};     % polynomial
        K = compKernel(data,obj.traindata,obj.kernel,obj.p1,obj.p2);
        
      case {'rbf','nrbf'};       % Radial basis function, a.k.a. gaussian
        K = compKernel(data,obj.traindata,obj.kernel,obj.p1);
        
      otherwise
        K = compKernel(data,obj.traindata,obj.kernel);
    end
    
    probs = K * obj.weights(1:end-1) + obj.weights(end);
    
    % post is just the sign and does not have a probabilistic interpretation
    post = zeros(size(probs,1),2);
    post(:,1) = (probs < 0);
    post(:,2) = (probs > 0);
    
    post = dataset(post);
    
  end
  
  function [m,desc] = getmodel(obj)
    % return the parameters wrt a class label in some shape
    
    m = {obj.primal}; % only one vector for kernelmethod
    desc = {'unknown'};
    
  end
  
end
end
