function MVNR = randmn(mu,S,N)
% RANDMN : generate a random number from multi-variate normal distribution N(mu,S)
%
% -- Usage
% MVNR = randmn(mu,S,N)
%
% -- Input
% mu : average vector (d*1)
%  S : covariance matrix (d*d)(should a be symmetric positive definite matrix)
%  N : the number of sample generated
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

d=size(mu,1); % dimension

if cond(S) > 1e+10  
    MVNR=zeros(d,N);
else
    L=chol(S);
    NR=randn(d,N);
    MVNR= L'* NR + repmat(mu,[1,N]);
end
