function [Ppos,Pneg]=participation_coef_sign(W,Ci)
%PARTICIPATION_COEF_SIGN     Participation coefficient
%
%   [Ppos Pneg] = participation_coef_sign(W,Ci);
%
%   Participation coefficient is a measure of diversity of intermodular
%   connections of individual nodes.
%
%   Inputs:     W,      undirected connection matrix with positive and
%                       negative weights
%
%               Ci,     community affiliation vector
%
%   Output:     Ppos,   participation coefficient from positive weights
%
%               Pneg,   participation coefficient from negative weights
%
%   Reference: Guimera R, Amaral L. Nature (2005) 433:895-900.
%
%
%   2011, Mika Rubinov, UNSW

%   Modification History:
%   Mar 2011: Original
%   Sep 2012: Fixed treatment of nodes with no negative strength
%             (thanks to Alex Fornito and Martin Monti)


n=length(W);                                %number of vertices

Ppos = pcoef( W.*(W>0));
Pneg = pcoef(-W.*(W<0));

    function P=pcoef(W_)
        S   = sum(W_,2);                    %strength
        Gc  = (W_~=0)*diag(Ci);             %neighbor community affiliation
        Sc2 = zeros(n,1);                   %community-specific neighbors

        for i = 1:max(Ci);
            Sc2 = Sc2 + (sum(W_.*(Gc==i),2).^2);
        end

        P = ones(n,1) - Sc2./(S.^2);
        P(isnan(P)) = 0;
        P(~P) = 0;                            %p_ind=0 if no (out)neighbors
    end
end
