
function [alpha,bias,margin] = rfda_svm(data,labels,kernel,C)

% RFDA_SVM performs regularised Fischer's discriminant analysis (or linear discriminant analysis) 
% for classification
%
%  Use as
%    1) [alpha,bias,margin] = rfda_svm(data,labels,kernel,C)
%    2) As part of SVMMETHOD object (see svmmethod.m) - RECOMMENDED
%
%  INPUT
%           data    - input data features
%           labels  - their class labels
%           kernel  - kernel matrix
%           C       - SVM's hyperparameter (default 1)
%
%  OUTPUT
%           alphas  - SVM coefficients (in a dual form)
%           bias    - bias term
%           margin  - margin size

% Pawel Herman, 2009


[n_data]= size(data,1);

alpha = zeros(n_data,1);
bias = 0;

if nargin < 4
    C = 1;
end

if all(sort(unique(labels))~=[-1 ; 1])
    error('Labels should be transformed to values -1 and 1');
end

ell = n_data;
ellplus = (sum(labels) + ell)/2;
yplus = 0.5*(labels + 1);
ellminus = ell - ellplus;
yminus = yplus - labels;
rescale = ones(ell,1)+labels*((ellminus-ellplus)/ell);
plusfactor = 2*ellminus/(ell*ellplus);
minusfactor = 2*ellplus/(ell*ellminus);
B = diag(rescale) - (plusfactor * yplus) * yplus' - (minusfactor * yminus) * yminus';
A = B * kernel + C * eye(ell,ell);

% A * alpha = labels
alpha = A \ labels;  % or simply inv(A) * labels
bias = -0.25 * (alpha' * kernel * rescale) / (ellplus * ellminus);

margin = 0;

% in order to fit the framework of svmmethod, which performs additional
% operation of multiplying alphas by labels (assigning signs) required for
% SVMs, here this multiplication is performed to cancel out with the
% corresponding operation in svmmethod
%--------------------------------------------------------------------------
alpha = alpha .* labels;
