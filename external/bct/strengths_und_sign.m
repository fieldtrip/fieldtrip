function [Spos,Sneg,vpos,vneg] = strengths_und_sign(W)
%STRENGTHS_UND_SIGN        Strength and weight
%
%   [Spos Sneg] = strengths_und_sign(W);
%   [Spos Sneg vpos vneg] = strengths_und_sign(W);
%
%   Node strength is the sum of weights of links connected to the node.
%
%   Inputs:     W,              undirected connection matrix with positive
%                               and negative weights
%
%   Output:     Spos/Sneg,      nodal strength of positive/negative weights
%               vpos/vneg,      total positive/negative weight
%
%
%   2011, Mika Rubinov, UNSW

%   Modification History:
%   Mar 2011: Original


n = length(W);              %number of nodes
W(1:n+1:end) = 0;           %clear diagonal
Spos = sum( W.*(W>0));      %positive strengths
Sneg = sum(-W.*(W<0));      %negative strengths

if nargout>2
    vpos = sum(Spos);       %positive weight
    vneg = sum(Sneg);       %negative weight
end