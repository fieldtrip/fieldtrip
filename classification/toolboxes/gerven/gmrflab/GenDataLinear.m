function [X,Y,beta_c,noise] = GenDataLinear(m,beta,v)


X         = randn(m, length(beta));
X         = X./repmat(sqrt(sum(X.^2,1)),size(X,1),1);
noise     = randn(m,1)/sqrt(v);
beta_c    = beta(:);
Y         = X*beta_c+noise;