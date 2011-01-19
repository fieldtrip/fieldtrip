function w = hh_sam_weight(Rinv,L)
% hh_sam_weight - This routine computes the weight vector of SAM
% 
% Usage: w = hh_sam_weight(Rinv,L)
%
% $AUTHOR: Copyright (C) by Hung Dang$
% $DATE: Thu Sep 30 23:32:12 MDT 2010$
% $LOG$
%  
% Revision 1.2 Thu Sep 30 23:32:30 MDT 2010, hungptit
% Update the comment

%% Compute the optimum orientation
X = Rinv * L;
P = L' * X;
Q = X' * X;
[V,D] = eig(P,Q);
d = diag(D);
max_pos = find(d == max(d));
vopt = V(:,max_pos);
vopt = vopt / norm(vopt);

%% Compute the weight vector
g = L * vopt;
w = Rinv * g / (vopt' * P * vopt);
% $$$ s = w' * y;
% $$$ plot(s);