function w = hh_smvb_weight(Rinv,L)
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
P = L' * Rinv * L;
Q = L' * L;
[V,D] = eig(Q,P);
d = diag(D);
max_pos = find(d == max(d));
vopt = V(:,max_pos);
vopt = vopt / norm(vopt);

%% Compute the weight vector
g = L * vopt;
w = Rinv * g / (vopt' * P * vopt);
% $$$ s = w' * y;
% $$$ plot(s);