function cost = hh_vsmvb_map(data,compress,decompress,node_sizes,R, Rinv)
% hh_smvb - Find topographic map for VSMVB
%
% $AUTHOR: Copyright (C) by Hung Dang$
% $ID: hh_smvb$
% $DATE: Fri Sep 17 22:03:18 MDT 2010$
%
% Revision 1.1 Fri Sep 17 22:03:52 MDT 2010, hungptit
% First update 
%     

% Welcome message
fprintf('\n------ %s -------\n','Compute the tomographic map of VSMVB');

% Init variables
NNODE = length(decompress);
cost = zeros(NNODE,1);

% Compute the compressed tomographic map of V-WNMVB
t0 = clock;
for nnode = 1:NNODE
    npos = decompress(nnode);
    L = hh_leadfield_new(data,compress,node_sizes,npos);
    X = Rinv * L;
    P = L' * X;
    Q = X' * X;
    % Cost function of V-WNMVB
    cost(nnode) = 1 / sum(eig(Q,P));
end

% Display the running time
t = etime(clock,t0);
fprintf('Total running time: %f seconds\n---------------------------\n',t);


