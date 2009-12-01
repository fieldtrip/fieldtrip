classdef svmmethod < classifier
%SVMMETHOD different svm methods
%
% It is a wrapper for different types of SVM: L2-SVM, L1-SVM, SMO-SVM and RFD. They
% are referred to as @svm_l2_slow, @svm_l1_slow (alternatively, @svm_km_l2 and @svm_km_l1), 
% and @rfda, respectively. Options include the type of kernel, 
% kernel parameters and regularization constant C. If the latter is a vector then 
% all values are tried using an optimizer that uses a validator and evaluation criterion.
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
% myproc = clfproc({standardizer svmmethod('method',@svm_km_l2,'C',90:10:100,'validator',crossvalidator('cvfolds',0.8,'verbose',true,'procedure',{standardizer svmmethod('method',@svm_km_l2,'kernel','rbf','kerparam',1)}),'kernel','rbf','kerparam',1)});
%
% SEE ALSO:
%   svm_l2_slow.m
%   svm_l1_slow.m
%   svm_km_l2.m
%   svm_km_l1.m
%   rfda_svm.m
% 
% Copyright (c) 2008, Pawel Herman
%
% $Log: svmmethod.m,v $
%

    properties
      method = @svm_km_l2;
      kernel = 'linear'; % 'poly','rbf'('gaussian'),'htrbf' - see calc_kernel.m
      kerparam = 1;      % each kernel is associated with native parameters - see calc_kernel.m
      C = nan;
      validator;% = crossvalidator('procedure',clfproc({}));
      criterion = 'accuracy';
      alpha;
      bias;
      wv;  %weight vector in primal form
      margin;
      traindata;
      sv_model;   %struct
      % .sv
      % .weights
      % .b
      ratio4estplatt = 0;     %>=0 && <=1 (0-default)
      % if >0 then estimation is done on
      % randomly selected subset of the data
      platt_sigmoid = [];  %with fields .A and .B
      
    end
    
    methods
      function obj = svmmethod(varargin)
        
        obj = obj@classifier(varargin{:});
        
      end
      function obj = train(obj,data,design)
        % simply stores input data and design
        
        % remove missing data
        data = data(~isnan(design),:);
        design = design(~isnan(design));
        
        % check for consistency
        [data,design] = obj.check_input(data,design);
        
        if ~exist('platt_sigmoid','var')
          obj.platt_sigmoid = [];
          obj.platt_sigmoid.A = [];
          obj.platt_sigmoid.B = [];
        end
        
        if ~isnan(any(obj.C)) && ~isscalar(obj.C)
          obj = optimize(obj,data,design,'variables','C','values',obj.C,'validator', ...
            obj.validator,'criterion',obj.criterion);
        end
        
        if iscell(data), error('svmmethod does not take multiple datasets as input'); end
        
        if isnan(obj.nclasses), obj.nclasses = max(design(:,1)); end
        
        if obj.nclasses ~= 2, error('svm only makes binary classifications'); end
        
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
        
        [obj.alpha,obj.bias,obj.margin] = obj.method(data,labels,K,obj.C);
        obj.traindata = data;
        obj.alpha = obj.alpha.*labels;  %explicitly assign sign to alphas
        
        % weight vector in primal form - only for linear svm
        obj.wv = 0;
        if strcmp(obj.kernel,'linear') && ~strcmp(obj.method,'@rfda')
          for j=1:size(data,1)
            obj.wv = obj.wv + obj.alpha(j)*data(j,:);
          end
        end
        
        obj.sv_model.sv = obj.traindata(obj.alpha~=0,:);
        obj.sv_model.weights = obj.alpha(obj.alpha~=0);
        obj.sv_model.b = obj.bias;
        
        % evaluation of sigmoid for platt's probabilistic outputs on a subset of training set
        if obj.ratio4estplatt > 0
          [data_vld,labels_vld] = stratified_division(data,labels,obj.ratio4estplatt);
          
          obj = estplatt(obj,data_vld,labels_vld);
        end
        
        
      end
      function post = test(obj,data)
        
        data = obj.check_input(data);
        
        % deal with empty data
        if size(data,2) == 0
          
          % random assignment
          post = rand([size(data,1) obj.nclasses]);
          post = double(post == repmat(max(post,[],2),[1 obj.nclasses]));
          
          return
        end
        
        probs = svm_eval(data,obj.sv_model.sv,obj.sv_model.weights,obj.sv_model.b,obj.kernel,obj.kerparam);
        
        if ~isempty(obj.platt_sigmoid.A) && ~isempty(obj.platt_sigmoid.B)
          % probabilistic outputs if platts sigmoid has been estimated
          probs = platt_svmproboutput(probs,obj.platt_sigmoid.A,obj.platt_sigmoid.B);
          post = [ 1-probs probs];
        else
          % probs is just the sign and does not have a probabilistic interpretation
          post = zeros(size(probs,1),2);
          post(:,1) = (probs < 0);
          post(:,2) = (probs > 0);
        end
        
      end
      
      function obj = estplatt(obj,data,design)
        if iscell(data), error('classifier does not take multiple datasets as input'); end
        
        % transform elements of the design matrix to class labels
        if unique(design(:,1))==[-1; 1]
          labels = design(:,1);
        elseif unique(design(:,1))==[1; 2]
          labels = design(:,1);
          labels(design(:,1) == 1) = -1;
          labels(design(:,1) == 2) = 1;
        else
          error('Design matrix should contain labels -1,1 or 1,2');
        end
        
        svm_out = svm_eval(data,obj.sv_model.sv,obj.sv_model.weights,obj.sv_model.b,obj.kernel,obj.kerparam);
        
        [obj.platt_sigmoid.A,obj.platt_sigmoid.B] = platt_sigmoidest(svm_out,labels);
        
      end
      
      function m = getmodel(obj,label,dims)
        % return the parameters wrt a class label in some shape
        
        m = obj.wv; % only one vector for svmmethod
        
        if nargin > 2 && numel(m) == prod(dims)
          m = reshape(m,dims);
        end
        
      end
      
    end
end
