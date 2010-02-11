function s = invSMW(X, a, b)
% Nsamp << Nfeat 
%    (A+X'BX)^(-1) = iA - iA*X'*(iB+X*iA*X')^(-1)*X*iA
%

%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

[Nsamp, Nfeat] = size(X);

ia = 1./a;
ib = 1./b;
iAX = ia(:,ones(1,Nsamp)) .* X';
XiAX = X*iAX;

C = diag(ib) + XiAX;
iC = inv(C);

for jj = 1 : Nfeat
    ic(jj) = iAX(jj,:)*iC*iAX(jj,:)';
end
s = ia - ic';
