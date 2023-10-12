function [Hpos,Hneg] = diversity_coef_sign(W, Ci)
%DIVERSITY_COEF_SIGN     Shannon-entropy based diversity coefficient
%
%   [Hpos Hneg] = diversity_coef_sign(W,Ci);
%
%   The Shannon-entropy based diversity coefficient measures the diversity
%   of intermodular connections of individual nodes and ranges from 0 to 1.
%
%   Inputs:     W,      undirected connection matrix with positive and
%                       negative weights
%
%               Ci,     community affiliation vector
%
%   Output:     Hpos,   diversity coefficient based on positive connections
%               Hneg,   diversity coefficient based on negative connections
%
%   References: Shannon CE (1948) Bell Syst Tech J 27, 379-423.
%               Rubinov and Sporns (2011) NeuroImage.
%
%
%   2011-2012, Mika Rubinov, U Cambridge

%   Modification History:
%   Mar 2011: Original
%   Sep 2012: Fixed treatment of nodes with no negative strength
%             (thanks to Alex Fornito and Martin Monti)


n = length(W);                                  %number of nodes
m = max(Ci);                                    %number of modules

Hpos = entropy(W.*(W>0));
Hneg = entropy(-W.*(W<0));

    function H = entropy(W_)
        S = sum(W_,2);                          %strength
        Snm = zeros(n,m);                       %node-to-module degree
        for i = 1:m                             %loop over modules
            Snm(:,i) = sum(W_(:,Ci==i),2);
        end
        pnm = Snm ./ S(:,ones(1,m));
        pnm(isnan(pnm)) = 0;
        pnm(~pnm) = 1;
        H = -sum(pnm.*log(pnm),2)/log(m);
    end
end
