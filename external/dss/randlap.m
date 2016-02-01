function X = randlap(n,N); 
% RANDLAP  Laplace distributed random numbers.
%   s = RANDLAP(n, M)
%     n  dimension
%     N  sample size

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

X0=log(rand(n,N)).*((rand(n,N)<.5)*2-1);
for t=1:n
  X0(t,:) = X0(t,:) - mean(X0(t,:));
  X(t,:)=X0(t,:)/std(X0(t,:),1);
end

