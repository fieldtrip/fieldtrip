function [Pcorrect] = calc_percor(errTable);
% Calculate percent correct from an error table (a confusion matrix)
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.
  
Nsamp = sum(errTable(:));
Ncor  = sum(diag(errTable));

Pcorrect = Ncor/Nsamp * 100;  % percent 

