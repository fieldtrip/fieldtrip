function prior =  priorstandard(dim, c)

prior = struct('mean',zeros(dim,1), 'cov', c*eye(dim));

