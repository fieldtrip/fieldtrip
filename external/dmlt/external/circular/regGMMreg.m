function [fv,grad] = regGMMreg(param,theta,X,lambda)
mu=param(1);
beta=param(2:end);
u = X * beta;
lx = 2*atan(u); % link function
% g=X'*(theta-mu-lx);
% dg=-2*X'*(X./repmat((1+u.^2),1,size(beta,1)));
% fv=g'*g+ lambda * sum(beta.*beta);
% grad=2*dg'*g+ 2*lambda * beta;
% grad=[-2*sum(X*g);grad];
% return
% 
g=X'*sin(theta-mu-lx);
dg=-2*X'*(X./repmat((1+u.^2)./cos(theta-mu-lx),1,size(beta,1)));
fv=g'*g+ lambda * sum(beta.*beta);
grad=2*dg'*g+ 2*lambda * beta;
grad=[-2*sum((X.*repmat(cos(theta-mu-lx),1,size(beta,1)))*g);grad];