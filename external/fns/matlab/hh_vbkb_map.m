function cost = hh_vbkb_map(data,compress,decompress,node_sizes,R, Rinv)
% hh_vbkb_map - Find topographic map for V-BKB
%
% $author: Copyright by Hung Dang$
% $Id: hh_vbkb_map.m$
% $Date: June 27, 2009$
% 
% $LOG$
% 
% Revision 1.1 Fri Sep 17 17:55:37 MDT 2010, hungptit
% Revise the old routine based on the new formula for VBKB where
% the regularised inverse is used in the final equation.
% 

% Welcome message
fprintf('\n------ %s -------\n',['Compute the tomographic map of the ' ...
                    'Borgiotti-Kaplan beamformer']);

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
    % Cost function of VBKB
    cost(nnode) = trace(Pinv ./ Omega);
end

% Display the running time
t = etime(clock,t0);
fprintf('Total running time: %f seconds\n---------------------------\n',t);

