function p = smlog2prob(q),
% Calculate probability of softmax function from log probability 
% 
% -- Input
% q = log(p) before normalization (i.e. sum(exp(p)) need not to be 1.)
% [Ndata Ncls]
%
% 2009/06/05 OY
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

[Ndata,Ncls] = size(q);

maxq = max(q,[],2);
qadj = q - repmat(maxq,[1,Ncls]);

p = exp(qadj) ./ repmat(sum(exp(qadj),2),[1,Ncls]);

