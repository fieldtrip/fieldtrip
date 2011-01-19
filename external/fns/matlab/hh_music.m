function cost = hh_music(data,compress,decompress,node_sizes,R,NSOURCE,Nk)
% hh_music - Calculate the cost values for MUSIC algorithm. This
% routine is the kernel for RAP-MUSIC algorithm
%
% Usage: cost = hh_music(data,compress,decompress,node_sizes,R,NSOURCE,Nk);
%
% $Id: hh_music.m$
% $Date: December 04, 2008$
% $author: Copyright by Hung Dang$
%
% Revision 1.1 Thu Aug  5 11:02:55 MDT 2010, hungptit
% Created based on the algorithm described in Mosher et. al. (1992)
% This approach uses the SVD to convert the optimization problem
% into symmetric eigen problem.    

% Welcome message
fprintf('\n------ %s -------\n',['Compute the tomographic map ' ...
                    'of MUSIC']);

% Find the signal subspace
[U,S,V] = svd(R);
Es = U(:,1:NSOURCE);
if (isempty(Nk))
    Esk = Es;
else
    Esk = Nk * Es;
end

% Update NNODE and cost
NNODE = length(decompress);
cost = zeros(NNODE,1); 

% Calculate topographic map of MUSIC
t0 = clock;
% Find the source here
if (isempty(Nk))
    for nnode = 1:NNODE
        npos = decompress(nnode);
        L = hh_leadfield_new(data,compress,node_sizes,npos);
        [U,S,V] = svd(L,0);
        X = Esk' * U;
        A = X' * X;
        % Cost value of MUSIC
        cost(nnode) = max(eig(A));
    end
else
    for nnode = 1:NNODE
        npos = decompress(nnode);
        L = hh_leadfield_new(data,compress,node_sizes,npos);
        Y = Nk * L;
        [U,S,V] = svd(Y,0);
        X = Esk' * U;
        A = X' * X;
        % Cost value of MUSIC
        cost(nnode) = max(eig(A));
    end
end

% Display the running time
t = etime(clock,t0);
fprintf('Total running time: %f seconds\n---------------------------\n',t);
