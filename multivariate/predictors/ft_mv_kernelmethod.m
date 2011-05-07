classdef ft_mv_kernelmethod < ft_mv_predictor
%FT_MV_KERNELMETHOD abstract kernel method class
%
% Note: this is not the most efficient implementation since the kernel is
% recomputed each time. 
%
% Copyright (c) 2008, Marcel van Gerven, Jason Farquhar


  properties
    
    kernel = 'linear'; % kernel; linear, polynomial or rbf
    p1 = nan;
    p2 = nan;
    C = nan;
    
    weights
    f
    J
    traindata
    primal
     
  end
  
  methods
    
    function obj = ft_mv_kernelmethod(varargin)
      
      obj = obj@ft_mv_predictor(varargin{:});
      
    end
    
    function obj = train(obj,X,Y)
      % simply stores input data and design
      
      % multiple datasets
      if iscell(X) || iscell(Y)
        obj = ft_mv_ndata('mvmethod',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      % convert data to double
      if ~isa(X(1),'double'), X = double(X); end
      if ~isa(Y(1),'double'), Y = double(Y); end
     
      % missing data
      if any(isnan(X(:))) || any(isnan(Y(:))), error('method does not handle missing data'); end
     
      % multiple outputs
      if size(Y,2) > 1
        obj = ft_mv_noutput('mvmethod',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      if max(Y(:)) > 2, error('kernel method only makes binary classifications'); end
      
      % transform elements of the design matrix to class labels
      idx1 = (Y==1);
      Y(idx1) = -1;
      Y(~idx1) = 1;

      K = obj.getkernel(X,X);
   
      % C parameter
      if isnan(obj.C), obj.C = .1*(mean(diag(K))-mean(K(:))); end      
      if obj.verbose, fprintf('C parameter was set to %.2g\n',obj.C); end

      % store training data
      obj.traindata = X;
   
      % kernelmethod specific call
      [obj.weights,obj.f,obj.J] = obj.estimate(K,Y);
            
      % weights in primal form
      obj.primal = 0;
      for j=1:size(X,1), obj.primal = obj.primal + obj.weights(j)*X(j,:); end
      
    end
    
    function Y = test(obj,X)
      
      probs = obj.getkernel(X,obj.traindata) * obj.weights(1:end-1) + obj.weights(end);
      
      % post is just the sign and does not have a probabilistic interpretation
      Y = zeros(size(probs,1),2);
      Y(:,1) = (probs < 0);
      Y(:,2) = (probs > 0);
      
    end
    
    function [m,desc] = model(obj)
      % return the parameters wrt a class label in some shape
      
      m = {obj.primal}; % only one vector for kernelmethod
      desc = {'primal form parameters; positive values indicate condition 2'};
      
    end
    
    function K = getkernel(obj,X,Z)
      
      switch lower(obj.kernel)
        
        case {'poly','npoly'};     % polynomial
          
          if isnan(obj.p1), obj.p1 = 2/3; end
          if isnan(obj.p2), obj.p2 = 1; end
          
          K = compKernel(X,Z,obj.kernel,obj.p1,obj.p2);
          
        case {'rbf','nrbf'};       % Radial basis function, a.k.a. gaussian
          
          if isnan(obj.p1), obj.p1 = .1*(mean(diag(X))-mean(X(:))); end
          
          K = compKernel(X,Z,obj.kernel,obj.p1);
          
        otherwise % linear kernel
          K = compKernel(X,Z,obj.kernel);
      end
      
    end
    
  end
  
  methods(Abstract)
    
    estimate(obj,K,Y);
    
  end
  
end
