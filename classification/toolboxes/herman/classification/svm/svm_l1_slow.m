
function [alpha,bias,margin] = svm_l1_slow(data,labels,kernel,C,epsilon)

% SVM_L1_SLOW implements an SVM classifier with quadratic norm of penalization terms 
% It is a slower realization than svm_km_l1.m (uses Matlab QP solver).
%
%  Use as
%    1) [alpha,bias,margin] = svm_l1_slow(data,labels,kernel,C,epsilon)
%    2) As part of SVMMETHOD object (see svmmethod.m) - RECOMMENDED
%
%  INPUT
%           data    - input data features
%           labels  - their class labels
%           kernel  - kernel matrix
%           C       - SVM's hyperparameter (default 1)
%           epsilon - constant for boundary support vectors (default 1e-12)
%
%  OUTPUT
%           alpha   - SVM coefficients (in a dual form)
%           bias    - bias term
%           margin  - margin size
%       


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

% add small numbers to diagonal 
lambda = 1e-12;
H = kernel.*(labels*labels') + lambda*eye(size(kernel));

Aeq = labels';
beq = 0;
f = -ones(n_data,1);       % alpha

LB = zeros(n_data,1);      % 0 <= alpha
x0 = zeros(n_data,1);      % starting point

if length(C) == 1,
    UB = C*ones(n_data,1);
elseif length(C) == n_data,
    UB=C;
else
    error('regularisation parameter C should be a scalar or a vector with the number of elements corresponding to the number of data vectors');
end

% Matlab Optimization Toolbox - QP
dispoff = optimset('Display','off');
[alpha fval exflag]= quadprog(H, f, [],[],Aeq, beq, LB, UB, x0, dispoff);

% all legitimate SVs
alpha_all = alpha;
alpha = zeros(length(alpha_all),1);
index_sv = find( alpha_all > epsilon);
alpha(index_sv) = alpha_all(index_sv);

% boundary SVs (f(x)=+/-1):  0+e < alpha < C-e
index_bound = find( (alpha_all > epsilon) & (alpha_all < (C - epsilon)));
n_svbound = length(index_bound);

% bias
bias = sum( labels(index_bound) - H(index_bound,index_sv) * alpha(index_sv) .* labels(index_bound) ) / n_svbound;

% margin
margin = 2 / sqrt( alpha(index_sv)' * H(index_sv,index_sv) * alpha(index_sv) );
