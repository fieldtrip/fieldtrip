function cost = hh_smvb_map(data,compress,decompress,node_sizes,R,Rinv)
% hh_smvb - Find topographic map for SMVB
%
% $AUTHOR: Copyright by Hung Dang$
% $ID: hh_smvb$
% $DATE: June 27, 2009$
%
% $LOG$
% Revision 1.0 Thu Sep 16 23:25:11 MDT 2010, hungptit
% First update 

%% Welcome message
fprintf('\n------ %s -------\n','Compute the tomographic map of SMVB');

%% Calculate covariance matrix and other related parameters
NNODE = length(decompress);
cost = zeros(NNODE,1);

%% Calculate cost values for all nodes
t0 = clock;
for nnode = 1:NNODE
    npos = decompress(nnode);
    L = hh_leadfield_new(data,compress,node_sizes,npos);
    X = Rinv * L;
    P = L' * X;
    Q = L' * L;
    % Cost function of SMVB
    cost(nnode) = max(eig(Q,P));
end        

% Display the running time
t = etime(clock,t0);
fprintf('Total running time: %f seconds\n---------------------------\n',t);

%% EOF
