function [VIn,MIn] = partition_distance(Cx,Cy)
%PARTITION_DISTANCE     Distance between community partitions
%
%   This function quantifies the distance between pairs of community
%   partitions with information theoretic measures.
%
%   VIn       = partition_distance(Cx,Cy)
%   [VIn MIn] = partition_distance(Cx,Cy)
%
%   Inputs:         Cx,         Community affiliation vector X
%                   Cy,         Community affiliation vector Y
%
%   Outputs:        VIn,        Normalized variation of information
%                   MIn,        Normalized mutual information
%
%   (Definitions:
%       VIn = [H(X) + H(Y) - 2MI(X,Y)]/log(n)
%       MIn = 2MI(X,Y)/[H(X)+H(Y)]
%    where H is entropy, MI is mutual information and n is number of nodes)
%
%   Reference: Meila M (2007) J Multivar Anal 98, 873-895.
%
%   2011, Mika Rubinov, UNSW

%   Modification History:
%   Mar 2011: Original

%#ok<*ASGLU>

n = numel(Cx);                                  %n

[dum,dum,Cx]  = unique(Cx);
[dum,dum,Cy]  = unique(Cy);
[dum,dum,Cxy] = unique([Cx(:) Cy(:)],'rows');

Px  = hist(Cx,1:max(Cx))/n;                     %P(x)
Py  = hist(Cy,1:max(Cy))/n;                     %P(y)
Pxy = hist(Cxy,1:max(Cxy))/n;                   %P(x,y)

Hx  = -sum(Px.*log(Px));                        %H(x)
Hy  = -sum(Py.*log(Py));                        %H(y)
Hxy = -sum(Pxy.*log(Pxy));                      %H(x,y)

VIn = (2*Hxy-Hx-Hy)/log(n);                     %VIn
MIn = 2*(Hx+Hy-Hxy)/(Hx+Hy);                    %MIn