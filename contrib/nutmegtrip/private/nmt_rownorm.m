function norm = nmt_rownorm(A)
% NORM = NMT_ROWNORM(A)
% Treat each row of A as vector and take norm.

norm = sqrt(sum(A.^2,2));
