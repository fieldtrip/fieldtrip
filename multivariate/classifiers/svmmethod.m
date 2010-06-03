classdef svmmethod < classifier
%SVMMETHOD different svm methods
%
% It is a wrapper for different types of SVM: L2-SVM, L1-SVM, SMO-SVM and RFD. They
% are referred to as @svm_l2_slow, @svm_l1_slow (alternatively, @svm_km_l2 and @svm_km_l1), 
% and @rfda, respectively. Options include the type of kernel, 
% kernel parameters and regularization constant C. If the latter is a vector then 
% all values are tried using an optimizer t§  hat uses a validator and evaluation criterion.
% 
% REMARK:
%     1) l2 and l1 norms refer to the penalisation component
%     2) It is possible to evaluate probabilistic outputs using Platt's approach
%       - either in a standard processing pipeline by specifying ratio4estplatt > 0 and using proportion of
%         training dataset for estimation of the sigmoid (during training)
%       - or, by calling estplatt method with some external validation set and then testing
%
%
% Options:
% 'method'  : SVM method used
% 'kernel'  : type of kernel (see calc_kernel.m) - default 'linear'
% 'kerparam': kernel parameters (see calc_kernel.m) - default 1
% 'C' : regularization parameter; if a vector then the optimal value will be computed using some validation procedure
% 'ratio4estplatt': proportion of data for estimating Platt's probabilistic outputs (default 0, see below in object's properties)
%
% REQUIRES:
% 1) MATLAB Optimization Toolbox (QP)
% 2) svm_km_l1 and svm_km_l2 require 2 external functions for QP: monqp.m and monqpCinfty.m 
%
%
% EXAMPLE:
%
% definition of the classification procedure while optimizing C
%
% myproc = mva({standardizer svmmethod('method',@svm_km_l2,'C',90:10:100,'validator',crossvalidator('cvfolds',0.8,'verbose',true,'procedure',{standardizer svmmethod('method',@svm_km_l2,'kernel','rbf','kerparam',1)}),'kernel','rbf','kerparam',1)});
%
% SEE ALSO:
%   svm_l2_slow.m
%   svm_l1_slow.m
%   svm_km_l2.m
%   svm_km_l1.m
%   rfda_svm.m
% 
% Copyright (c) 2008, Pawel Herman

    properties
      
      method = @svm_km_l2;
      kernel = 'linear'; % 'poly','rbf'('gaussian'),'htrbf' - see calc_kernel.m
      kerparam = 1;      % each kernel is associated with native parameters - see calc_kernel.m
      C = nan;
      validator;% = crossvalidator('procedure',mva({}));
      criterion = 'accuracy';
      ratio4estplatt = 0;     %>=0 && <=1 (0-default)
 
    end
    
    methods
      
      function obj = svmmethod(varargin)
        
        obj = obj@classifier(varargin{:});
        
      end
      
      function p = estimate(obj,X,Y)
        % simply stores input data and design
            
        idx = obj.labeled(X);
        
        nclasses = obj.nunique(Y);
        
        % remove any missing data
        data = X(idx,:);
        design = Y(idx,:);
        
        p.platt_sigmoid = [];
        p.platt_sigmoid.A = [];
        p.platt_sigmoid.B = [];
      
        if nclasses ~= 2, error('svm only makes binary classifications'); end
        
        % transform elements of the design matrix to class labels
        labels = design(:,1);
        labels(design(:,1) == 1) = -1;
        labels(design(:,1) == 2) = 1;
        
        % calculate kernel matrix
        K = calc_kernel(data,data,obj.kernel,obj.kerparam);
        
        % regularization parameter
        if isnan(obj.C), obj.C = .1*(mean(diag(K))-mean(K(:))); end
        
        if obj.verbose
          fprintf('regularisation parameter was set to %.2f\n',obj.C);
        end
        
        [p.alpha,p.bias,p.margin] = obj.method(data,labels,K,obj.C);
        p.traindata = data;
        p.alpha = p.alpha.*labels;  %explicitly assign sign to alphas
        
        % weight vector in primal form - only for linear svm
        p.wv = 0;
        if strcmp(obj.kernel,'linear') && ~strcmp(obj.method,'@rfda')
          for j=1:size(data,1)
            p.wv = p.wv + p.alpha(j)*data(j,:);
          end
        end
        
        p.sv_model.sv = p.traindata(p.alpha~=0,:);
        p.sv_model.weights = p.alpha(p.alpha~=0);
        p.sv_model.b = p.bias;
        
        % evaluation of sigmoid for platt's probabilistic outputs on a subset of training set
        if obj.ratio4estplatt > 0
          
          [data_vld,labels_vld] = stratified_division(data,labels,obj.ratio4estplatt);
          
          svm_out = svm_eval(data,p.sv_model.sv,p.sv_model.weights,p.sv_model.b,obj.kernel,obj.kerparam);
          
          [p.platt_sigmoid.A,p.platt_sigmoid.B] = platt_sigmoidest(svm_out,labels);
          
        end
        
      end
      
      function post = map(obj,X)
    
        % remove missing data
        X = X(obj.labeled(X),:);
                
        probs = svm_eval(X,obj.params.sv_model.sv,obj.params.sv_model.weights,obj.params.sv_model.b,obj.kernel,obj.kerparam);
        
        if ~isempty(obj.params.platt_sigmoid.A) && ~isempty(obj.params.platt_sigmoid.B)
          
          % probabilistic outputs if platts sigmoid has been estimated
          probs = platt_svmproboutput(probs,obj.params.platt_sigmoid.A,obj.params.platt_sigmoid.B);
          post = [ 1-probs probs];
          
        else
          
          % probs is just the sign and does not have a probabilistic interpretation
          post = zeros(size(probs,1),2);
          post(:,1) = (probs < 0);
          post(:,2) = (probs > 0);
        
        end
        
      end
     
      
      function [m,desc] = getmodel(obj)
        % return the parameters
        
        m = {obj.params.wv}; % only one vector for svmmethod
        desc = {'primal form parameters; positive values indicate condition two'};
        
      end
      
    end
end
