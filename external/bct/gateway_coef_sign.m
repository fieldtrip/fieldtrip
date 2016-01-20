function [GWpos,GWneg] = gateway_coef_sign(W,Ci,centtype)
%GATEWAY_COEF_SIGN     Gateway coefficient
%
%   [Gpos,Gneg] = gateway_coef_sign(W,Ci,centtype);
%
%   Gateway coefficient is a variant of participation coefficient. Similar 
%   to participation coefficient, gateway coefficient measures the 
%   diversity of intermodular connections of individual nodes, but this is
%   weighted by how critical these connections are to intermodular
%   connectivity (e.g., if a node is the only connection between it's
%   module and another module, it will have a higher gateway coefficient). 
%
%   Inputs:     W,        undirected connection matrix with positive and
%                         negative weights
%
%               Ci,       community affiliation vector
%
%               centtype, centrality measure to use
%                         1 = Node Strength
%                         2 = Betweenness Centrality
%
%   Output:     Gpos,     gateway coefficient for positive weights
%
%               Gneg,     gateway coefficient for negative weights
%
%   Reference: Vargas ER, Wahl LM. Eur Phys J B (2014) 87:1-10.
%
%   Jeff Spielberg, Boston University

%   Modification History:
%   May 2015: Original (adapted from participation_coef_sign.m)


[~,~,Ci] = unique(Ci);
n=length(W);                                %number of vertices
W(1:(n+1):end) = 0;

GWpos = gcoef( W.*(W>0));
GWneg = gcoef(-W.*(W<0));

    function GW = gcoef(W_)
        S     = sum(W_,2);                    %strength
        Gc    = (W_~=0)*diag(Ci);             %neighbor community affiliation
        Sc2   = zeros(n,1);                   %community-specific neighbors
        ksm   = zeros(n,1);                   % Weighting by extra-modular importance
        centm = zeros(n,1);                   % Weighting by intra-modular importance
        
        switch centtype % Which centrality measure to use?
            case 1      % Node Strength
                cent = S;
            case 2      % Betweenness Centrality
                L = weight_conversion(W_,'lengths');
                cent = betweenness_wei(L);
        end
        
        for i = 1:max(Ci);
            ks = sum(W_.*(Gc==i),2);
            Sc2 = Sc2 + (ks.^2);
            for j = 1:max(Ci);
                ksm(Ci==j) = ksm(Ci==j) + ks(Ci==j)/sum(ks(Ci==j)); % Calculate extra-modular weights
            end
            centm(Ci==i) = sum(cent(Ci==i));                        % Calculate intra-modular weights
        end
        
        centm = centm./max(centm); % Re-map to be between 0 & 1
        gs = (1-(ksm.*centm)).^2;  % Calculate total weights
        
        GW = ones(n,1) - (Sc2./(S.^2).*gs);
        GW(isnan(GW)) = 0;
        GW(~GW) = 0;                            %p_ind=0 if no (out)neighbors
    end
end