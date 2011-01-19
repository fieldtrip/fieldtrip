function cost = hh_vwnmvb_mapsvd(data,compress,decompress,node_sizes,R,Sigma_inv)
% hh_smvb - Find topographic map for SMVB
%
% $author: Copyright by Hung Dang$
% $Id: hh_smvb$
% $Date: June 27, 2009$
% Revision 1.1 
    
% Welcome message
fprintf('\n------ %s -------\n',['Compute the tomographic map of ' ...
                    'VWNMVB using svd approach']);

% Init variables
NNODE = length(decompress);
cost = zeros(NNODE,1);

% Compute the compressed tomographic map of V-WNMVB
t0 = clock;
for nnode = 1:NNODE
    npos = decompress(nnode);
    L = hh_leadfield_new(data,compress,node_sizes,npos);
    X = Sigma_inv * L;
    [U,S,V] = svd(X,0);
    A = U' * R * U;
    % Cost function of VWNMVB
    cost(nnode) = sum(eig(A));
end

% Display the running time
t = etime(clock,t0);
fprintf('Total running time: %f seconds\n---------------------------\n',t);

% EOF
