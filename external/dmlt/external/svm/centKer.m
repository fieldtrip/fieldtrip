function K=centKer(K)

%function K=centKer(K)
%
% Centeres the matrix K
%
%INPUTS
% K = the kernel matrix to be centered
%
%OUTPUTS
% Kc = the centered kernel matrix
%
%
%For more info, see www.kernel-methods.net


% original kernel matrix stored in variable K
% output uses the same variable K
% K is of dimension ell x ell
% D is a row vector storing the column averages of K
% E is the average of all the entries of K
ell = size(K,1);
D = sum(K) / ell;
E = sum(D) / ell;
J = ones(ell,1) * D;
K = K - J - J' + E;
