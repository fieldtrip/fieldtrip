function cost = hh_vwnmvb_map(data,compress,decompress,node_sizes,R, Rinv)
% hh_smvb - Find topographic map for SMVB
%
% $author: Copyright by Hung Dang$
% $Id: hh_smvb$
% $Date: June 27, 2009$
%
% Revision  1.2 Fri Sep 17 22:05:52 MDT 2010, hungptit
% Fix the program based on the update formula
%
% Revision 1.1 
% First update 
    
% Welcome message
fprintf('\n------ %s -------\n','Compute the tomographic map of VWNMVB');

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
    cost(nnode) = sum(eig(P,Q));
end

% Display the running time
t = etime(clock,t0);
fprintf('Total running time: %f seconds\n---------------------------\n',t);


