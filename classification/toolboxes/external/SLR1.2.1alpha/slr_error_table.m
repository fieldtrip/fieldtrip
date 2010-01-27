function Table = slr_error_table(label_true, label_est)
% Make an error table
% The i-th row of 'Table' ecorresponds to histogram of estimated label when true label is
% i.
%
% -- Usage
% Table = error_table(label_true, label_est);
%
% -- Input
% label_true : true label     [Nsamp*1]
% label_est  : estimated (predicted) label [Nsamp*1] 
%
% 2006/09/07  OY
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

label_pair = [label_true, label_est];

label_names = unique(label_pair(:));
Nclass = length(label_names);

for ii = 1 : Nclass
    for jj = 1 : Nclass
          ix = find(label_pair(:,1) == label_names(ii) & label_pair(:,2) == label_names(jj));        
          Table(ii, jj) = length(ix);
    end
end