function W = hh_lcmv_weight(Rinv,L)
% hh_lcmv_weight - This routine computes the weight matrix of the
% LCMV 
%
% Usage: W = hh_lcmv_weight(Rinv,L)
%
% $AUTHOR: Copyright (C) by Hung Dang$
% $DATE: Thu Sep 30 23:57:41 MDT 2010$
% $LOG$ 
% Revision 1.1 Thu Sep 30 23:58:51 MDT 2010, hungptit
% First create
% 

P = L' * Rinv * L;
Pinv = inv(P);
W = Rinv * L * Pinv;