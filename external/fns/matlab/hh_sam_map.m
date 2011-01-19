function cost = hh_sam_map(data,compress,decompress,node_sizes,R,Rinv)
% hh_smvb - Find topographic map for SMVB
%
% $AUTHOR: Copyright by Hung Dang$
% $ID: hh_smvb$
% $DATE: June 27, 2009$
%
% $LOG$
% 
% Revision 1.2 Thu Sep 16 19:41:37 MDT 2010, hungptit
% Update the function based on the new formula, in which the
% inverse of the covariance matrix is replaced by its regularised 
% inverse in the final equation. 
% 
% Revision 1.1 Wed Aug  4 14:27:26 MDT 2010
% Create the routine to compute the tomographic map of SAM using
% the formulation in Vrba and Robinson (1999) (Need to verify this)
% 

%% Welcome message
fprintf('\n------ %s -------\n','Compute the tomographic map of SAM');

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
    Q = X' * X;
    % Cost function of SMVB
    cost(nnode) = max(eig(P,Q));
end        

% Display the running time
t = etime(clock,t0);
fprintf('Total running time: %f seconds\n---------------------------\n',t);

%% EOF
