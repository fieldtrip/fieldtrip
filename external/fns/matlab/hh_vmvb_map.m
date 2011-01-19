function cost = hh_vmvb_map(data,compress,decompress,node_sizes, R, Rinv)
% hh_vmvb_map - Find output power map of VMVB
%
% Usage: cost = hh_vmvb_map(data,compress,decompress,node_sizes, R, Sigma_inv)
%
% $AUTHOR: Copyright (C) by Hung Dang$
% $ID: hh_vmvb_map$
% $DATE: Tue Sep 21 12:43:56 MDT 2010$
%
% $LOG: hh_vmvb_map.m$
%
% Revision 1.1 Wed Jul 21 23:24:07 MDT 2010 hungptit
% This is a normalized output power of VMVB with Pout =
% trace(inv(P)) / trace(inv(L'*L)), this metrics has been proposed
% in Dang (2010).
%
    
% Welcome message
fprintf('\n------ %s -------\n',['Compute the normalized output power ' ...
                    'of VMVB');

% Init parameters
NNODE = length(decompress);
cost = zeros(NNODE,1);

% Calculate cost values for all nodes
t0 = clock;
for nnode = 1:NNODE
    npos = decompress(nnode);
    L = hh_leadfield_new(data,compress,node_sizes,npos);
    X = Rinv * L;
    P = L' * X;
    Q = L' * L;
    Pinv = inv(P);
    Qinv = inv(Q);
    % The normalized output power of VMVB
    cost(nnode) = trace(Pinv) / trace(Qinv);
end

% Display the running time
t = etime(clock,t0);
fprintf('Total running time: %f seconds\n---------------------------\n',t);
