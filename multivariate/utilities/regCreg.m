function [fv,grad] = regCreg(beta,theta,X,lambda)

u = X * beta;
lx = 2*atan(u); % link function
dlx = 2./(1+u.^2); % derivative of link function
S = mean(sin(theta - lx));
C = mean(cos(theta - lx));

mu = angle(1i*S+C);
DG = X' * (dlx .* sin(theta - mu - lx));

fv = -sum(cos(theta-mu-lx)) + lambda * sum(beta.*beta);
grad = 2 * lambda * beta - DG;