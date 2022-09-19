classdef svm < dml.method
% SVM support vector machine.
%
% DESCRIPTION
% Linear support vector machine classifier
%
%   EXAMPLE
%   X = rand(10,20); Y = [1 1 1 1 1 2 2 2 2 2]';
%   m = dml.svm
%   m = m.train(X,Y);
%   Z = m.test(X);
%
%   DEVELOPER
%   Jason Farquhar (j.farquhar@donders.ru.nl)

  properties
    
    X % training data
    
    primal % weights in primal form

    dual % weights in dual form
    
    C % regularization parameter
    
    kernel = 'linear'; % type of kernel
    
    Ktrain % precomputed kernel for training data
    Ktest  % precomputed kernel for test data
    
    distance = false; % return distance from decision boundary instead of class label
    native   = false; % uses native Bioinformatics toolbox SVM implementation if true
    issqrtmK = false;
  end
  
  methods
    
    function obj = svm(varargin)

      obj = obj@dml.method(varargin{:});

    end
    
    function obj = train(obj,X,Y)
   
      opts = {'verb' -1};
        
      % handle multiple datasets
      if iscell(X)
        obj = dml.ndata('method',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      if obj.native
        obj.Ktrain = fitcsvm(X,Y);
        return
      end
      
      if obj.restart || isempty(obj.dual)
        
        obj.X      = X;
        obj.Ktrain = compKernel(X,X,obj.kernel,obj.Ktrain);
        obj.Ktest  = [];
        
        if isempty(obj.C)
          if obj.issqrtmK
            % the actual kernel will be K*K'
            diagK = sum(obj.Ktrain.^2,2);
            meanK = mean(reshape(obj.Ktrain*obj.Ktrain',[],1));
          else
            diagK = diag(obj.Ktrain);
            meanK = mean(obj.Ktrain(:));
          end
          obj.C = .1*(mean(diagK)-meanK);
          if obj.verbose, fprintf('using default C=%.2f\n',obj.C); end
        end
        
        if obj.issqrtmK
          opts = cat(2, opts, {'issqrtmK' true});
        end
        obj.dual = l2svm_cg(obj.Ktrain,2*(Y-1)-1,obj.C, opts);
        
      else
        
        % facilitate warm restarts
        obj.dual = l2svm_cg(obj.Ktrain,2*(Y-1)-1,obj.C, opts); %,'alphab',obj.dual);
        
      end
      
      obj.primal = 0;
      for j = 1:size(X,1) 
        obj.primal = obj.primal + obj.dual(j)*X(j,:);
      end
      
    end
    
    function Y = test(obj,X)
      
      if obj.native
        C = svmclassify(obj.Ktrain,X);
        Y = zeros(size(C,1),2);
        Y(C + (0:2:(2*(size(C,1)-1)))')=1;
        return
      end
      
      if isempty(obj.Ktest)
        obj.Ktest = compKernel(X,obj.X,'linear');
      end
      
      probs = obj.Ktest * obj.dual(1:end-1) + obj.dual(end);
      
      Y = zeros(numel(probs),2);
      if obj.distance
      
        Y(probs<0,1) = abs(probs(probs<0));
        Y(probs>0,2) = abs(probs(probs>0));
      
      else
        % post is just the sign and does not have a probabilistic interpretation
        Y(:,1) = (probs < 0);
        Y(:,2) = (probs > 0);
      end
      
    end

    function m = model(obj)
    % returns
    %
    % m.primal primal form SVM parameters
    
      m.primal = obj.primal;
      
    end
  
  end
 
  
end
