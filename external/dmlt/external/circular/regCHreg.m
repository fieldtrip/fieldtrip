function [fv,grad] = regCHreg(param,theta,X,lambda)

mu=param(1);
beta=param(2:end);

u = X * beta;
lx = 2*atan(u); % link function
dlx = 2./(1+u.^2); % derivative of link function
DG = X' * (dlx .* sin(theta - mu - lx));
DG2 = 2*X' * (dlx .* sin(2*(theta - mu - lx)));

fv = -sum(cos(theta-mu-lx)+cos(2*(theta-mu-lx))) + lambda * sum(beta.*beta);
grad = 2 * lambda * beta - DG-DG2;
grad=[-sum(sin(theta-mu-lx)+2*sin(2*(theta-mu-lx)));grad];