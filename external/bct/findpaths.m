function [Pq,tpath,plq,qstop,allpths,util] = findpaths(CIJ,sources,qmax,savepths)
%FINDPATHS      Network paths
%
%   [Pq,tpath,plq,qstop,allpths,util] = findpaths(CIJ,sources,qmax,savepths);
%
%   Paths are sequences of linked nodes, that never visit a single node
%   more than once. This function finds all paths that start at a set of 
%   source nodes, up to a specified length. Warning: very memory-intensive.
%
%   Inputs:     CIJ,        binary (directed/undirected) connection matrix
%               qmax,       maximal path length
%               sources,    source units from which paths are grown
%               savepths,   set to 1 if all paths are to be collected in
%                           'allpths'
%
%   Outputs:    Pq,         3D matrix, with P(i,j,q) = number of paths from
%                           'i' to 'j' of length 'q'.
%               tpath,      total number of paths found (lengths 1 to 'qmax')
%               plq,        path length distribution as a function of 'q'
%               qstop,      path length at which 'findpaths' is stopped
%               allpths,    a matrix containing all paths up to 'qmax'
%               util,       node use index
%
%   Note that Pq(:,:,N) can only carry entries on the diagonal, as all
%   "legal" paths of length N-1 must terminate.  Cycles of length N are
%   possible, with all vertices visited exactly once (except for source and
%   target). 'qmax = N' can wreak havoc (due to memory problems).
%
%   Note: Weights are discarded.
%   Note: I am certain that this algorithm is rather inefficient -
%   suggestions for improvements are welcome.
%
%   Olaf Sporns, Indiana University, 2002/2007/2008/2010

%   2010 version:
%   -- a bug affecting the calculation of 'util' was fixed -- thanks to
%      Steve Williams
%   -- better pre-allocation for 'npths'
%   -- note that this code assumes a directed graph as input - calculation
%      of paths and 'util' indices can be easily adapted to undirected
%      graphs. 

% ensure CIJ is binary...
CIJ = double(CIJ~=0);

% initialize some variables
N = size(CIJ,1); K = sum(sum(CIJ));
pths = [];
Pq = zeros(N,N,qmax);
util = zeros(N,qmax);

% this code is for pathlength = 1
% paths are seeded from 'sources'
q = 1;
for j=1:N
    for i=1:length(sources)
        is = sources(i);
        if (CIJ(is,j) == 1)
            pths = [pths [is j]'];
        end;
    end;
end;

% calculate the use index per vertex (for paths of length 1)
util(1:N,q) = util(1:N,q)+hist(reshape(pths,1,size(pths,1)*size(pths,2)),1:N)';
% now enter the found paths of length 1 into the pathmatrix Pq
for np=1:size(pths,2)
    Pq(pths(1,np),pths(q+1,np),q) = Pq(pths(1,np),pths(q+1,np),q) + 1;
end;

% begin saving all paths
if (savepths==1)
    allpths = pths;
end;
if (savepths~=1)
    allpths = [];
end;

% initialize
npthscnt = K;

% "big loop" for all other pathlengths 'q'
% ----------------------------------------------------------------------
for q=2:qmax

    % to keep track of time...
    disp(['current pathlength (q) = ',num2str(q),'   number of paths so far (up to q-1)= ',num2str(sum(sum(sum(Pq))))])

    % old paths are now in 'pths'
    % new paths are about to be collected in 'npths'
    % estimate needed allocation for new paths
    len_npths = min(ceil(1.1*npthscnt*K/N),100000000);
    npths = zeros(q+1,len_npths);

    % find the unique set of endpoints of 'pths'
    endp = unique(pths(q,:));
    npthscnt = 0;

    for ii=1:length(endp)  % set of endpoints of previous paths
        i = endp(ii);
        % in 'pb' collect all previous paths with 'i' as their endpoint
        [pa,pb] = find(pths(q,:) == i);
        % find the outgoing connections from 'i' ("breadth-first")
        nendp = find(CIJ(i,:)==1);
        % if 'i' is not a dead end
        if (~isempty(nendp))
            for jj=1:length(nendp)   % endpoints of next edge
                j = nendp(jj);
                % find new paths - only "legal" ones, i.e. no vertex is visited twice
                pb_temp = pb(sum(j==pths(2:q,pb),1)==0);
                % add new paths to 'npths'
                npths(:,npthscnt+1:npthscnt+length(pb_temp)) = [pths(:,pb_temp)' ones(length(pb_temp),1)*j]';
                npthscnt = npthscnt+length(pb_temp);
                % count new paths and add the number to 'P'
                Pq(1:N,j,q) = Pq(1:N,j,q)+(hist(pths(1,pb_temp),1:N))';
            end;
        end;
    end;

    % note: 'npths' now contains a list of all the paths of length 'q'
    if (len_npths>npthscnt)
        npths = npths(:,1:npthscnt);
    end;

    % append the matrix of all paths
    if (savepths==1)
        allpths = [allpths; zeros(1,size(allpths,2))];
        allpths = [allpths npths(:,1:npthscnt)];
    end;

    % calculate the use index per vertex (correct for cycles, count
    % source/target only once)
    util(1:N,q) = util(1:N,q) + hist(reshape(npths(:,1:npthscnt),1,size(npths,1)*npthscnt),1:N)' - diag(Pq(:,:,q));
    % eliminate cycles from "making it" to the next level, so that
    % 'pths' contains all the paths that have a chance of being continued
    if (~isempty(npths))
        pths = npths(:,npths(1,:)~=npths(q+1,:));
    else
        pths = [];
    end;

    % if there are no 'pths' paths left, end the search
    if (isempty(pths))
        qstop = q;
        tpath = sum(sum(sum(Pq)));
        plq = reshape(sum(sum(Pq)),1,qmax);
        return;
    end;

end;  % q
% ----------------------------------------------------------------------
qstop = q;

% total number of paths
tpath = sum(sum(sum(Pq)));

% path length distribution
plq = reshape(sum(sum(Pq)),1,qmax);
