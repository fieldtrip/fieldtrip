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
% PARAMETERS:
% weights
% f
% J
% traindata
% primal  
%
%
% EXAMPLE:
%
% % definition of the classification procedure 
%
% myproc = mva({ kernelmethod('method',@l2svm_cg)});
%
% SEE ALSO:
% l2svm_cg.m
% rkls.m
% klr_cg.m
%
% Copyright (c) 2008, Marcel van Gerven, Jason Farquhar


  properties
    method = @l2svm_cg;
    kernel = 'linear';
    p1 = nan;
    p2 = nan;
    C = nan;
    
  end
  
  methods
    
    function obj = kernelmethod(varargin)
      
      obj = obj@classifier(varargin{:});
      
    end
    
    function p = estimate(obj,X,Y)
      % simply stores input data and design
      
      if obj.nunique(Y) ~= 2, error('l2svm only makes binary classifications'); end
      
      % transform elements of the design matrix to class labels
      targets = Y(:,1);
      targets(Y == 1) = -1;
      targets(Y == 2) = 1;
      
      % set input for our classifier; compute kernel matrix
      
      switch lower(obj.kernel)
        
        case {'poly','npoly'};     % polynomial
          
          if isnan(obj.p1), obj.p1 = 2/3; end
          if isnan(obj.p2), obj.p2 = 1; end
          
          K = compKernel(X,X,obj.kernel,obj.p1,obj.p2);
          
        case {'rbf','nrbf'};       % Radial basis function, a.k.a. gaussian
          
          if isnan(obj.p1), obj.p1 = .1*(mean(diag(X))-mean(X(:))); end
          
          K = compKernel(X,X,obj.kernel,obj.p1);
          
        otherwise
          K = compKernel(X,X,obj.kernel);
      end
      
      % regularization parameter
      if isnan(obj.C), obj.C = .1*(mean(diag(K))-mean(K(:))); end
      
      if obj.verbose
        fprintf('regularization parameter was set to %.2g\n',obj.C);
      end
      
      [p.weights,p.f,p.J] = obj.method(K,targets,obj.C,'verb',-1);
      p.traindata = X;
      
      % weights in primal form
      p.primal = 0;
      for j=1:size(X,1), p.primal = p.primal + p.weights(j)*X(j,:); end
      
    end
    
    function Y = map(obj,X)
      
      % deal with empty data
      if size(X,2) == 0
        
        % random assignment
        Y = rand([size(X,1) obj.nclasses]);
        Y = double(Y == repmat(max(Y,[],2),[1 obj.nclasses]));
        
        return
      end
      
      switch lower(obj.kernel)
        
        case {'poly','npoly'};     % polynomial
          K = compKernel(X,obj.params.traindata,obj.kernel,obj.p1,obj.p2);
          
        case {'rbf','nrbf'};       % Radial basis function, a.k.a. gaussian
          K = compKernel(X,obj.params.traindata,obj.kernel,obj.p1);
          
        otherwise
          K = compKernel(X,obj.params.traindata,obj.kernel);
      end
      
      probs = K * obj.params.weights(1:end-1) + obj.params.weights(end);
      
      % post is just the sign and does not have a probabilistic interpretation
      Y = zeros(size(probs,1),2);
      Y(:,1) = (probs < 0);
      Y(:,2) = (probs > 0);
      
    end
    
    function [m,desc] = getmodel(obj)
      % return the parameters wrt a class label in some shape
      
      m = {obj.params.primal}; % only one vector for kernelmethod
      desc = {'primal form parameters; positive values indicate condition 2'};
      
    end
    
  end
end
