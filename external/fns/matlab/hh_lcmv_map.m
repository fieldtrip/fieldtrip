function cost = hh_lcmv_map(data,compress,decompress,node_sizes, R, Rinv)
% hh_lcmv_map - Find topographic map of LCMV
%
% Usage: cost = hh_lcmv_map(data,compress,decompress,node_sizes, R, Sigma_inv)
%
% $Author: Copyright (C) by Hung Dang$
% $Id: hh_lcmv_map$
% $Date: June 27, 2009$
% $Log: hh_lcmv_map.m$
%
% Revision 1.2 Thu Sep 16 16:03:44 MDT 2010, hungptit
% Update routine based on the new equation (only used the
% regularised inverse in the final equation): Pout = trace(Pinv) 
%
% Revision 1.1 Wed Jul 21 23:24:07 MDT 2010 hungptit
% Copy the lcmv routine from the HHSIM code 
%
    
% Welcome message
fprintf('\n------ %s -------\n','Compute the tomographic map of LCMV');

% Calculate covariance matrix and other related parameters
NNODE = length(decompress);
cost = zeros(NNODE,1);

% Calculate cost values for all nodes
t0 = clock;
for nnode = 1:NNODE
    npos = decompress(nnode);
    L = hh_leadfield_new(data,compress,node_sizes,npos);
    X = Rinv * L;
    P = L' * X;
    Q = X' * X;
    Pinv = inv(P);
    Omega = Pinv * Q * Pinv;       
    % Cost function of VMVB
    cost(nnode) = trace(Pinv) / trace(Omega);
end

% Display the running time
t = etime(clock,t0);
fprintf('Total running time: %f seconds\n---------------------------\n',t);
