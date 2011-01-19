function cost = hh_mcsmb_map(data,compress,decompress,node_sizes,R,Rinv,L0)
% hh_smvb - Find topographic map for multiple-correllated source
% model beamformer.
% 
% Usage: cost = hh_mcsmb_map(data,compress,decompress,node_sizes,R,Rinv,L0)
%
% $AUTHOR: Copyright by Hung Dang$
% $ID: hh_mcsmb_map.m$
% $DATE: Mon Sep 27 15:29:01 MDT 2010$
%
% $LOG$
% Revision 1.0 Mon Sep 27 15:29:12 MDT 2010, hungptit
% First update 

% Welcome message
fprintf('\n------ %s -------\n',['Compute the tomographic map of ' ...
                    'multiple-correlated source model']);

% Calculate covariance matrix and other related parameters
NNODE = length(decompress);
cost = zeros(NNODE,1);

% Calculate cost values for all nodes
t0 = clock;
if (isempty(L0))
    % No constraint is given then we only use SMVB
    for nnode = 1:NNODE
        npos = decompress(nnode);
        L = hh_leadfield_new(data,compress,node_sizes,npos);
        X = Rinv * L;
        P = L' * X;
        Q = L' * L;
        % Cost function of SMVB
        cost(nnode) = max(eig(Q,P));
    end        
else
    % Constraint is given then compute the the output map of MCSMB    
    for nnode = 1:NNODE
        npos = decompress(nnode);
        L = hh_leadfield_new(data,compress,node_sizes,npos);
        G = [L,L0];
        [U,S,V] = svd(G,0);
        X = Rinv * U;
        A = U' * X;
        % Cost function of SMVB
        cost(nnode) = 1 / min(eig(A));
    end        
end
% Display the running time
t = etime(clock,t0);
fprintf('Total running time: %f seconds\n---------------------------\n',t);
