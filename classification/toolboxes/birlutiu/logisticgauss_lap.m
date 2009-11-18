function x = logisticgauss_lap(prior, feat, d)

% Laplace approximation

dim = size(feat,2);
options = optimset('GradObj', 'on', 'Hessian', 'on','Display', 'off');%, 'DerivativeCheck', 'on');
x = fminunc(@(w)(likeprior_all(w, prior, feat, d)), zeros(dim,1), options);

end


function [f, grad, hess] = likeprior_all(w, prior, feat, d)

nterms = size(feat,1); 
dim = length(w);
f=0; grad=zeros(dim,1); hess=zeros(dim);
for k = 1:nterms
    [f1, grad1, hess1] = likeprior_term(w, d(k) * feat(k,:)');
    
    f = f + f1;
    grad = grad + grad1;
    hess = hess + hess1;
end

if ~isempty(prior)
    [f1, g1, h1] = loggauss(w, prior.mean, prior.cov);
    f = f - f1;
    grad = grad - g1;
    hess = hess - h1;
end
end


function [f, grad, hess] = likeprior_term(w, expo)

wm = w'*expo;
f = log(1 + exp(wm));

% Gradient.
gr = 1 + exp(-wm);
grad = expo./gr;

% Hessian matrix
dimW = length(w);
hess = zeros(dimW);
for i =1:dimW
    for j =1:dimW
        if (i <= j) % the Hessian is symmetric
            hess(i,j) = (expo(i)*expo(j)*(gr-1))/(gr^2);
            hess(j,i) = hess(i,j);
        end
    end
end
end


function [p, grad, hess] = loggauss(x, Mu, A)
n = length(Mu);
% log(sqrt(2*pi)) = 0.9189
p = -0.5 * (x-Mu)' * A * (x-Mu) - n * 0.9189 - (0.5 * log(det(A)));
grad = -A * (x - Mu);
hess = -0.5 * (A + A');
end

