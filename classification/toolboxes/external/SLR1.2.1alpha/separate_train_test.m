function [ix_train, ix_test] = separate_train_test(label, Rtrain)
% Separate whole data into train and test data set so that number of
% samples in each class is balanced.
%
% [ix_train, ix_test] = separate_train_test(label, Rtrain)
%
% -- Input
% label : Label of each sample.
% Rtrain : Ratio of trainig data to whole data 
%
% 2006/09/29 OY
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

label = label(:);

% label_names (must be real value)
label_names = unique(label);
Nclass = length(label_names);

Nsamp_tot = length(label);
Nsamp_tr = floor(Nsamp_tot*Rtrain);
Nsamp_tr_cl = floor(Nsamp_tr/Nclass);
fprintf('Number of samples in training data set %d ...\n', Nsamp_tr_cl * Nclass); 
fprintf('Number of samples in test data set %d ...\n', Nsamp_tot - Nsamp_tr_cl * Nclass); 

ix_train = [];

for cc = 1 : Nclass
    
    ix_samp_cl = find(label == label_names(cc)); % samples index belonging to class 'cc'
    ix_cl_rnd = randperm(length(ix_samp_cl));
    
    ix_train = [ix_train; ix_samp_cl(ix_cl_rnd(1:Nsamp_tr_cl))];    % column vector
    
end

ix_train = sort([ix_train]);
ix_test = setdiff([1:Nsamp_tot]', ix_train);  






