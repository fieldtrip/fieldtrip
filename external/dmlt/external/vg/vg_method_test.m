function Y = vg_method_test(X, vg)

% Variational Garrote: calculates error for a given model
%   See file test_vg.m for an example of use
%
% required parameters (n input dimension, p samples)
%   X     : n x p (test set, input)
%   vg    : learned model created using vg_method_train

pt = size(X,2);

% normalize test set    
xt=X-mean(X,2)*ones(1,pt);
dx=sqrt(1/pt*sum(xt.^2,2));

Y = vg.v_mf*X;
