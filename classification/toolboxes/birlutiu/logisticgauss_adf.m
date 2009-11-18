function prior = logisticgauss_adf(prior, feat, d)

npoints = length(d);
for k = 1:npoints
    [mean_post, cov_post] = redDim_approx(prior.mean, prior.cov, d(k).*feat(k,:)');
    prior = struct('mean', mean_post ,'cov', cov_post);
end
end


function [m_post, cov_post, norm] = redDim_approx(m, C, a)

Ca = C*a;
aCa = Ca'*a;
am = a'*m;
sqa = sqrt(aCa);

%--- Gaussian quadrature
[norm, mom1_norm, mom2_norm] = mom_norm_gaussQuad((@(t)f_norm(t, am, sqa)), 5);

scalar1 = mom1_norm/sqa;
scalar2 = (mom2_norm-1)/aCa;
m_post = m + Ca*scalar1;
cov_post = sum_matrix_c(C, Ca, scalar1^2, scalar2);
end


function f = f_norm(w, am, sqa)
f = g(w*sqa + am);
end


function y = g(z)
y = 1 / ( 1 + exp(-z));
end


function [Z, mhat, vhat] = mom_norm_gaussQuad(fname, nQuad)

%w = [0.0000 0.0008 0.0191 0.1355 0.3446 0.3446 0.1355 0.0191 0.0008 0.0000];
%x = [-4.8595 -3.5818 -2.4843 -1.4660 -0.4849 0.4849 1.4660 2.4843 3.5818 4.8595]';

w = [0.0113 0.2221 0.5333 0.2221 0.0113];   % quadrature weights
x = [-2.8570 -1.3556 -0.0000 1.3556 2.8750]'; % quadrature points

fval = zeros(nQuad,1);
for i=1:nQuad
    fval(i) = feval(fname, x(i));
end
Z = w*fval;
mhat = diag(x*w)'*fval/Z;
vhat = diag(diag(x*w)*x')'*fval/Z;
end
