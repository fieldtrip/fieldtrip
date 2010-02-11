function [N, bin] = slr_hist_feature(IX_EFF, Nfeature) 
% Histogram of effective (survived) features
%
% Nfeature : # of features
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

ix_eff = []
for i = 1 : length(IX_EFF)
    ix_eff = [ix_eff; IX_EFF{i}];
end

[N, bin] = hist(ix_eff, Nfeature);
hist(ix_eff, Nfeature);
