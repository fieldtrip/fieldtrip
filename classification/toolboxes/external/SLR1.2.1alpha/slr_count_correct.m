function [Ncorrect, label_est, prob ,y] = slr_count_correct(label, X, w)
% Count the number of correct label 
%
% function [Ncorrect, label_est, prob ,y] = slr_count_correct(label, X, w)
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

y = X*w;
prob = 1 ./ (1+exp(-y));
label_est = prob > 0.5;

% ~xor  0 0 -> 1, 
%       0 1 -> 0
%       1 0 -> 0 
%       1 1 -> 1

[v] = ~xor(label, label_est);
Ncorrect = sum(v);
