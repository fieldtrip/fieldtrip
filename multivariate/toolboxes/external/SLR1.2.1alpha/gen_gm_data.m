function [t,x] = gen_gm_data(N, mu);
% Generate artificial binary-labeled data from Gaussian mixture model
% 
% N  : Number of samples to be generated
% mu : coordinate fo center of gravity for two feature vectors [2x(#feature)]
%      each row vector represents the coordinate. 
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

[Nclass,Nd] = size(mu);


r = rand(N,1); 
t = r > 0.5; % label

xtmp0 = randn(N,Nd) + repmat(mu(1,:),[N,1]);
xtmp1 = randn(N,Nd) + repmat(mu(2,:),[N,1]);

ix0 = find(t == 0);
ix1 = find(t == 1);


x(ix0,:) = xtmp0(ix0,:);
x(ix1,:) = xtmp1(ix1,:); %% data