
function [alpha,bias,margin] = svm_km_l1(data,labels,kernel,C,epsilon)

% SVM_KM_L1 implements an SVM classifier with linear norm of penalization terms
%
%  Use as
%    1) [alpha,bias,margin] = svm_km_l1(data,labels,kernel,C,epsilon)
%    2) As part of SVMMETHOD object (see svmmethod.m) - RECOMMENDED
%
%  INPUT
%           data    - input data features
%           labels  - their class labels
%           kernel  - kernel matrix
%           C       - SVM's hyperparameter (default 1)
%           eps     - threshold for zeroing alphas (default 1e-12)
%
%  OUTPUT
%           alphas  - SVM coefficients (in a dual form)
%           bias    - bias term
%           margin  - margin size
%
%
%  REQUIRES
%           this function uses monqp.m and monqpCinfty.m (look there for references)
%
%  SEE ALSO
%           svm_km_l2.m
%           svmmethod.m

% Pawel Herman, 2009

n_data = size(data,1);

alpha = zeros(n_data,1);
bias = 0;

if nargin < 4
    C = 1;
end
if nargin < 5
    epsilon = 1e-12;
end

if all(sort(unique(labels))~=[-1 ; 1])
    error('Labels should be transformed to values -1 and 1');
end

H = kernel.*(labels*labels');

Aeq = labels;
beq = 0;
f = ones(n_data,1);      

verbose = 0;
lambda = 1e-10;
[weights , bias , index_sv] = monqp(H,f,Aeq,beq,C,lambda,verbose,data,kernel);       

alpha = zeros(n_data,1);
alpha(index_sv) = weights;

% boundary SVs (f(x)=+/-1):  0+e < alpha < C-e
index_bound = find( (alpha > epsilon) & (alpha < (C - epsilon)));
n_svbound = length(index_bound);

bias = sum( labels(index_bound) - H(index_bound,index_sv) * alpha(index_sv) .* labels(index_bound) ) / n_svbound;

% margin
margin = 2 / sqrt( alpha(index_sv)' * H(index_sv,index_sv) * alpha(index_sv) );
