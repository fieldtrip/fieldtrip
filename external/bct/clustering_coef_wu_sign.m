function [C_pos,C_neg,Ctot_pos,Ctot_neg] = clustering_coef_wu_sign(W,coef_type)
%CLUSTERING_COEF_WU_SIGN     Multiple generalizations of the clustering coefficient 
%
%   [C_pos,C_neg,Ctot_pos,Ctot_neg] = clustering_coef_wu_sign(W,coef_type);
%
%   The weighted clustering coefficient is the average weight or intensity
%   of all triangles associated with each node.
%
%   Inputs:
%       W,          
%           Weighted undirected connection matrix
%
%       corr_type,
%           Desired type of clustering coefficient.
%           Options:  
%           1,  (default) Onnela et al. formula, used in original
%               clustering_coef_wu.m. Computed separately for positive &
%               negative weights.
%           2,  Zhang & Horvath formula, similar to Onnela formula except
%               denominator of Onnela formula relies on binarizing the
%               network whereas this denominator is based on weight value,
%               which reduces the sensitivity of this measure to the
%               weights directly connected to the node of interest.
%               Computed separately for positive & negative weights.
%           3,  Constantini & Perugini's generalization of the Zhang &
%               Horvath formula. This formula takes both positive &
%               negative weights into account simultaneously, & is
%               particularly sensitive to non-redundnacy in path
%               information based on sign (i.e., when two weights are
%               positive & one negative, or all three are negative, both of
%               which indicate that the weight of the third path is not
%               redundant information). Produces only one value.
%
%
%   Outputs: 
%       C_pos/C_neg,
%           Clustering coefficient vector for positive/negative weights.
%           For the third option, only one vector is outputted (as C_pos). 
%       Ctot_pos/Ctot_neg,
%           Mean clustering coefficient for positive and negative weights.
%
%   References: 
%       Onnela et al. (2005) Phys Rev E 71:065103
%       Zhang & Horvath (2005) Stat Appl Genet Mol Biol 41:1544-6115
%       Costantini & Perugini (2014) PLOS ONE 9:e88669
%
%
%   Contributor: Jeff Spielberg, Boston University, 2014-2015
%                (script based on clustering_coef_wu.m)

%
%   Modification History:
%   May 2014: Added computation of pos & neg weights separately & 
%             computation of mean coefficient (Jeff Spielberg)
%   May 2015: Added computation of Zhang & Horvath and Constantini & 
%             Perugini formulas (Jeff Spielberg)

if ~exist('coef_type','var')
    coef_type = 1;
end

n            = length(W);                   %number of nodes
W(1:n+1:end) = 0;

switch coef_type
    case 1
        W_pos                = W.*(W>0);
        K_pos                = sum(W_pos~=0,2);
        cyc3_pos             = diag((W_pos.^(1/3))^3);
        K_pos(cyc3_pos == 0) = inf;             %if no 3-cycles exist, make C=0 (via K=inf)
        C_pos                = cyc3_pos./(K_pos.*(K_pos-1));         %clustering coefficient
        Ctot_pos             = mean(C_pos);
        
        W_neg                = -W.*(W<0);
        K_neg                = sum(W_neg~=0,2);
        cyc3_neg             = diag((W_neg.^(1/3))^3);
        K_neg(cyc3_neg == 0) = inf;             %if no 3-cycles exist, make C=0 (via K=inf)
        C_neg                = cyc3_neg./(K_neg.*(K_neg-1));         %clustering coefficient
        Ctot_neg             = mean(C_neg);
    case 2
        W_pos    = W.*(W>0);
        cyc3_pos = zeros(n,1);
        cyc2_pos = zeros(n,1);
        for i = 1:n
            for j = 1:n
                for q = 1:n
                    cyc3_pos(i) = cyc3_pos(i)+(W_pos(j,i)*W_pos(i,q)*W_pos(j,q));
                    if j~=q
                        cyc2_pos(i) = cyc2_pos(i)+(W_pos(j,i)*W_pos(i,q));
                    end
                end
            end
        end
        cyc2_pos(cyc3_pos == 0) = inf;             %if no 3-cycles exist, make C=0 (via K=inf)
        C_pos                   = cyc3_pos./cyc2_pos;         %clustering coefficient
        Ctot_pos                = mean(C_pos);
        
        W_neg    = -W.*(W<0);
        cyc3_neg = zeros(n,1);
        cyc2_neg = zeros(n,1);
        for i = 1:n
            for j = 1:n
                for q = 1:n
                    cyc3_neg(i) = cyc3_neg(i)+(W_neg(j,i)*W_neg(i,q)*W_neg(j,q));
                    if j~=q
                        cyc2_neg(i) = cyc2_neg(i)+(W_neg(j,i)*W_neg(i,q));
                    end
                end
            end
        end
        cyc2_neg(cyc3_neg == 0) = inf;             %if no 3-cycles exist, make C=0 (via K=inf)
        C_neg                   = cyc3_neg./cyc2_neg;         %clustering coefficient
        Ctot_neg                = mean(C_neg);
    case 3
        cyc3         = zeros(n,1);
        cyc2         = zeros(n,1);
        
        for i = 1:n
            for j = 1:n
                for q = 1:n
                    cyc3(i) = cyc3(i)+(W(j,i)*W(i,q)*W(j,q));
                    if j~=q
                        cyc2(i) = cyc2(i)+(W(j,i)*W(i,q));
                    end
                end
            end
        end
        
        cyc2(cyc3 == 0) = inf;             %if no 3-cycles exist, make C=0 (via K=inf)
        C_pos           = cyc3./cyc2;         %clustering coefficient
        Ctot_pos        = mean(C_pos);
        C_neg           = nan(size(C_pos));
        Ctot_neg        = nan(size(Ctot_pos));
end
