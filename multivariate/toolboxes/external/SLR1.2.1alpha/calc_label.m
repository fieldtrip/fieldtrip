function [label_est, P] = calc_label(X, w);
% Estimate label for feature matrix X and linear weight w
%
% -- Input
% X : [Nsample, Nfeature]
% w : [Nfeature, Nclass]
%
% 2009/06/01 OY add header comments
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

Nclass = size(w,2);

if Nclass == 1,   % binary classification with parsimonious parametrization
    [label_est] = X * w > 0 ; 
    label_est = label_est + 1; % {0,1} -> {1,2} 
    P = 1 ./(1+exp(-X* w));
else
    [tmp, label_est] = max(X * w,[],2);
    eY = exp(X*w); % Nsamp*Nclass
    P = eY ./ repmat(sum(eY,2), [1, Nclass]); % Nsamp*Nclass
end