
function proboutput = platt_svmproboutput(svm_out,A,B)

% Auxiliary function for estimating Platt's probabilistic outputs for SVM 
%
%  USE
%       to be used within SVMMETHOD object (see also platt_sigmoidtest.m)

% Pawel Herman, 2009


% with assumption that svm's targets are -1 and 1
proboutput = 1./(1 + exp(A*svm_out+B));