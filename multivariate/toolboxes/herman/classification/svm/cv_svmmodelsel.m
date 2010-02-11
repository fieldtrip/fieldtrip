
function [cv_acc,opt_params,cv_array] = cv_svmmodelsel(data,labels,kernel_cfg,svm_cfg,nfolds,mtimes)


% CV_SVMMODELSEL implements cross-validation based selection of SVM's
% parameters (hyperparameter C and kernel parameters)
%
% INPUT
%     data    - input data features
%     labels  - their class labels
%      - a set of grid parameters to verify using cross-validation estimate of generalization error (or accuracy: 1 - error_rate) as a criterion
%        is passed in two structures: kernel_cfg and svm_cfg
%      - nfolds and mtimes specify the details of a cross-validation algorithm (mtimes x nfolds-fold CV)
%      
% OUTPUT 
%      cv_acc     - cross-validation estimates of classification accuracy for given sets of parameters 
%      opt_params - optimal set of parameters [C, kernel_params], i.e. for which the maximum classification accuracy 
%                   was obtained in cross-validation
%      cv_array   - full array with cross-validation accuracies and corresponding parameter sets
%       
% USE
%      - it is an auxiliary function used within SVMMETHOD object (see svmmethod.m)

% Pawel Herman, 2009

if nargin < 3
    kernel_cfg.kernel = 'linear';
    kernel_cfg.kerparam = 1;
end
if nargin < 4
    svm_cfg.C = [1 10 20 50 100 200];
    svm_cfg.method = 'svm_km_l2';
end
if nargin < 5
    nfolds = 5;
end
if nargin < 6
    mtimes = 10;
end

param_cfg = myvariation(svm_cfg.C,kernel_cfg.kerparam);
res = zeros(1,size(param_cfg,1));

for j=1:size(param_cfg,1)

    C_j = param_cfg(j,1);
    kerparam_j = param_cfg(j,2:end);

    local_proc = clfproc({ svmmethod('method',str2func(svm_cfg.method),'C',C_j,'kernel',kernel_cfg.kernel,'kerparam',kerparam_j) });
    cv = crossvalidator('procedure',local_proc,'cvfolds',nfolds,'randomize',true,'verbose',true);

    acc = zeros(mtimes,1);
    for i=1:mtimes
        cv = cv.validate(data,labels);
        acc(i) = cv.evaluate('metric','accuracy');
    end
    res(j) = mean(acc);
end

cv_array = [res' param_cfg];

[cv_acc optind] = max(res);
Csvm =  param_cfg(optind,1);  
kerparam = param_cfg(optind,2:end); 
opt_params = [{Csvm} {kerparam}];