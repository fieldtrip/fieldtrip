function [C, q]=core_periphery_dir(W,gamm,C)
%CORE_PERIPHERY_DIR     Core/periphery structure and core-ness statistic
%
%   C     = core_periphery_dir(W)
%   [C,q] = core_periphery_dir(W,gamm,C0)
%
%   The optimal core/periphery subdivision is a partition of the network
%   into two non-overlapping groups of nodes, a core group and a periphery
%   group, in a way that maximizes the number/weight of within core-group
%   edges, and minimizes the number/weight of within periphery-group edges.
%
%   The core-ness is a statistic which quantifies the goodness of the
%   optimal core/periphery subdivision.
%
%   Input:      W       directed (weighted or binary) connection matrix.
%               gamma,  core-ness resolution parameter (optional)
%                       gamma>1     detects small core/large periphery
%                       0<=gamma<1  detects large core/small periphery
%                       default is gamma=1
%
%   Outputs:    C,      binary vector of optimal core structure
%                       C = 1 represents nodes in the core
%                       C = 0 represents nodes in the periphery
%               q,      maximized core-ness statistic
%
%   Algorithm: A version of Kernighan-Lin algorithm for graph partitioning
%   used in community detection (Newman, 2006) applied to optimize a
%   core-structure objective described in Borgatti and Everett (2000).
%
%   Reference: Borgatti and Everett (2000) Soc Networks 21:375–395.
%              Newman (2006) Phys Rev E 74:036104, PNAS 23:8577-8582.
%              Rubinov, Ypma et al. (2015) PNAS 112:10032-7
%
%   2015, Mika Rubinov, U Cambridge

n = length(W);                              % number of nodes
W = double(W);                              % convert from logical
W(1:n+1:end) = 0;                           % clear diagonal
if ~exist('gamm','var')
    gamm = 1;
end
if ~exist('C','var')
    C = (rand(1,n)<0.5);
else
    C = logical(reshape(C,1,n));
end

% Methodological note: cf. community detection, the core-detection
% null model is not corrected for degree (to enable detection of hubs).
s = sum(W(:));
p = mean(W(:));
b = W - gamm*p;
B = (b+b.')/(2*s);                          % directed core-ness matrix
q = sum(sum(B(C,C))) - sum(sum(B(~C,~C)));  % core-ness statistic

f=1;                                        % loop flag
while f;
    f=0;
    Idx = 1:n;                              % initial node indices
    Ct = C;
    while any(Idx);
        Qt = zeros(1,n);                    % check swaps of node indices
        q0 = sum(sum(B(Ct,Ct))) - sum(sum(B(~Ct,~Ct)));
        Qt( Ct) = q0 - 2*sum(B( Ct, :),2);
        Qt(~Ct) = q0 + 2*sum(B(~Ct, :),2);
        
        %%% verification that the above update is equivalent to:
        % for u=Idx
        %     Ct(u) = ~Ct(u);
        %     Qt(u) = sum(sum(B(Ct,Ct))) - sum(sum(B(~Ct,~Ct)));
        %     Ct(u) = ~Ct(u);
        % end
        
        max_Qt = max(Qt(Idx));              % make swap with maximal
        u = find(abs(Qt(Idx)-max_Qt)<1e-10);% increase in core-ness
        u = u(randi(numel(u)));
        Ct(Idx(u)) = ~Ct(Idx(u));
        Idx(u)=[];                          % remove index from consideration
        
        if max_Qt-q>1e-10;                  % recompute core-ness statistic
            f = 1;
            C = Ct;
            q = sum(sum(B(C,C))) - sum(sum(B(~C,~C)));
        end
    end
end

q = sum(sum(B(C,C))) - sum(sum(B(~C,~C)));  % return core-ness statistic