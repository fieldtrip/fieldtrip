function [t,x] = gen_sin_data(N, mu, S);
% Generate artificial binary-labeled data from sin shaped boundary
% This code is only valid for two dimension feature space. 
%
% N  : Number of samples to be generated
% mu : coordinate fo center of gravity for two feature vectors [2x(#feature)]
%      each row vector represents the coordinate. 
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

[Nclass,Nd] = size(mu);

t = zeros(N,1);

Nsamp = floor(N/2);    % sample per each class

xtmp0 = randmn(mu(:,1),S,Nsamp);
xtmp1 = randmn(mu(:,2),S,Nsamp);

xtmp0 = xtmp0';
xtmp1 = xtmp1';

ix01 = find(xtmp0(:,2) > sin(xtmp0(:,1)));
t(ix01) = 1;
ix00 = find(xtmp0(:,2) < sin(xtmp0(:,1)));
t(ix00) = 0;

ix11 = find(xtmp1(:,2) > sin(xtmp1(:,1)));
t(ix11+Nsamp) = 1;
ix10 = find(xtmp1(:,2) < sin(xtmp1(:,1)));
t(ix10+Nsamp) = 0;



x = [xtmp0; xtmp1];