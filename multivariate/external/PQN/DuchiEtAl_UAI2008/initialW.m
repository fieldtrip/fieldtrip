function W = initialW(Sigma,lambda)

alpha = 1;
W     = (alpha - 1) * Sigma + (1-alpha) * diag(diag(Sigma)) + diag(lambda);

