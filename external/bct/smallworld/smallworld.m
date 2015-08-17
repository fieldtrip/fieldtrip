function [sw sw_node] = smallworld(W,type,n_rand)

% this function computes smallworld nes using the formula
% 
% SW=Cn/Ln, where Ln=L/Lr and Cn=C/Cr;
% 
% L = the characterstic path length or original network
% Lr = the characterstic path length of ranom network with same
%      degree/weight distirbtuion
% C = the mean clustering coefficient
% Cr = the mean clustering coefficient of ranom network with same
%      degree/weight distirbtuion
%
% W = the connectivity matrix
%
% type = how to treat the matrix ('binund', 'bindir', 'weiund',
%       'weidir'). 
%
% n_rand = the numder of random network from which to estimate Lr and Cr.
%         Default=100.
%
% The function can also compute a regional smallworldness value (sw_node),
% using the local path length and the local clustering ocefficients in the
% equivilent way as standard smallworldness
%
% (C) Mark Drakesmith, Cardiff University

if ~exist('n_rand','var')
    n_rand=100;
end

for i=1:n_rand+1
    if i==1 % this is the unrandomised netowrk
        disp('Computing paramaters for original network');
        W_iter=W;
    else  % do randomisation
        fprintf('Computing paramaters for random network %u\n',i-1);
        switch type
            case 'binund'
                W_iter=randomise_graph_bin_und(W);
            case 'bindir'
                W_iter=randomise_graph_bin_dir(W);
            case 'weiund'
                W_iter=randomise_graph_wei_und(W);
            case 'weidir'
                W_iter=randomise_graph_wei_dir(W);
        end
    end
    
    
    
    % compute char path
    switch type
        case {'binund','bindir'}
            dist=distance_bin(W_iter);
        case {'weiund','weidir'}
            dist=distance_wei(W_iter);
    end
    
    
    L_iter(i)=charpath(dist);
    
    % compute mean min path length for each node
    dist(isinf(dist))=NaN;
    L_iter_node(:,i)=nanmean(dist,2);
    
    % compute clustering coefficient (per node)
    switch type
        case 'binund'
            C_iter_node(:,i)=clustering_coef_bu(W_iter);
        case 'bindir'
            C_iter_node(:,i)=clustering_coef_bd(W_iter);
        case 'weiund'
            C_iter_node(:,i)=clustering_coef_wu(W_iter);
        case 'weidir'
            C_iter_node(:,i)=clustering_coef_wd(W_iter);
    end
    
    % compute mean clustering coefficient
    C_iter(i)=mean(C_iter_node(:,i));
end
disp('Finished randomisation. Computing smallworldness.');

% get median of path length / clust coeffs across randomisations

Cr_node=median(C_iter_node(:,2:end),2);
Cr=median(C_iter(2:end));

Lr_node=median(L_iter_node(:,2:end),2);
Lr=median(L_iter(:,2:end));  

% get original of path length / clust coeffs across randomisations
C_node=C_iter_node(:,1);
C=C_iter(:,1);  

L_node=L_iter_node(:,1);
L=L_iter(:,1);  
    
% compute normalised charpath / clust coeff
Ln=L./Lr;
Cn=C./Cr;

Ln_node=L_node./Lr_node;
Cn_node=C_node./Cr_node;

% finally compute smallworldness
sw=Cn./Ln;
sw_node=Cn_node./Ln_node;




       