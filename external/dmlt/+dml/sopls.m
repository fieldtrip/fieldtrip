classdef sopls < dml.method
% SOPLS sparse orthonormalized partial least squares.
%
%   DESCRIPTION
%
%   SOPLS estimates a model of the form Z = PX + C, Y = QZ.
%
%   lambda_min specifies the fraction of the maximal lambda value (computed
%   by glmnet) for which the output will be returned. 
%
%   alpha specifies the balance between ridge regression (alpha=0) and lasso
%   regression (alpha=1).
%
%   The number of steps taken to traverse the regularization path is given by nlambda.
% 
%   if nlambda=0 then standard partial least squares will be used.
%
%   pmax specifies the maximum number of variables included in the solution.
% 
%   nhidden specifies the number of used components. If it is decreased
%   during testing then the restricted number of components will be used
%   (only in case of sparse partial least squares!)
%
%   EXAMPLE
%   
%   load data;
%   p = dml.sopls('nhidden',3,'verbose',true);
%   p = p.train(X(1:40,1:10),I(1:40,:));
%   r = p.test(X(41:50,1:10));
%   for j=1:10, subplot(2,5,j); imagesc(reshape(r(j,:),[16 16])); colormap(gray); axis off, end
%   figure
%   for j=1:10, subplot(2,5,j); imagesc(reshape(I(40+j,:),[16 16])); colormap(gray); axis off, end
%
%   v = dml.enet.lambdapath(X,Y,'gaussian',50,1e-3);
%   m = dml.gridsearch('cv',dml.crossvalidator('type','split','stat','correlation','mva',dml.enet('family','gaussian','restart',false)),'vars','L1','vals',v,'verbose',true);
%   p = dml.sopls('nhidden',3,'verbose',true,'regularizer',m);
%   p = p.train(X(1:40,1:10),I(1:40,:));
%   r = p.test(X(41:50,1:10));
%   for j=1:10, subplot(2,5,j); imagesc(reshape(r(j,:),[16 16])); colormap(gray); axis off, end
%   figure
%   for j=1:10, subplot(2,5,j); imagesc(reshape(I(40+j,:),[16 16])); colormap(gray); axis off, end
%
%   REFERENCE
%  
%   When using this method please refer to the following:
%
%   van Gerven MAJ, Heskes T. Sparse Orthonormalized Partial Least Squares. In: BNAIC. 2010. 
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl), Tom Heskes (tomh@cs.ru.nl)


  properties
  
    nhidden = 5;      % number of hidden units; estimated using leave-one-out if an array
    
    Q             % output matrix
    P             % input matrix
    C             % hidden bias vector
    
    score         % leave-one-out score
    
    G             % trained ft_mv_glmnet objects

    method = 'sopls' % 'pls' (non-sparse), 'opls' (non-sparse orthonormal), 'sopls' (sparse opls)
    
    regularizer =  dml.enet('family','gaussian','verbose',true); % default regularizer
         
  end

  methods
        
    function obj = sopls(varargin)
      
      obj = obj@dml.method(varargin{:});
      
    end
    
    function obj = train(obj,X,Y)
        
      if ~isempty(obj.P) && obj.nhidden < size(obj.P,2)
        % used when performing a grid search
        % NOTE: grid search needs to count down from maximum nhidden
        obj.P = obj.P(:,1:obj.nhidden);
        obj.Q = obj.Q(:,1:obj.nhidden);
        if ~isempty(obj.G)
          obj.G = obj.G(1:obj.nhidden);
        end
        obj.C = obj.P(end,:) * obj.Q';
        
      else
      
        if strcmp(obj.method,'sopls') % sparse orthonormalized partial least squares
          
          [obj.Q,obj.P,bias,obj.G] = sopls(X',Y',obj.nhidden,obj.regularizer,obj.verbose);
          
          obj.C = bias * obj.Q';
          
        elseif strcmp(obj.method,'opls') % orthonormalized partial least squares
          
          [obj.Q,obj.P] = opls([X ones(size(X,1),1)]',Y',obj.nhidden);
          
          obj.C = obj.P(end,:) * obj.Q';
          obj.P = obj.P(1:(end-1),:);
          
        else % Matlab's partial least squares (SIMPLS algorithm)
          
          [XL,YL,XS,YS,beta,pctvar,mse,stats] = plsregress(X,Y,obj.nhidden);
          
          obj.P = stats.W;
          obj.Q = YL;
          obj.C = mean(Y,1) - mean(X,1)*obj.P*obj.Q';
          
        end
      
      end
      
    end
    
    function Y = test(obj,X)

      Y = bsxfun(@plus,X * obj.P(:,1:obj.nhidden) * obj.Q(:,1:obj.nhidden)',obj.C);
       
    end
    
    function m = model(obj)
      % returns
      %
      % m.weights regression weights (P*Q')
      % m.bias bias term
      
      m.weights = obj.P(:,1:obj.nhidden) * obj.Q(:,1:obj.nhidden)';
      m.bias = obj.C;
      
    end
    
  end
  
end
