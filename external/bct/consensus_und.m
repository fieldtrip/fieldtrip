function ciu = consensus_und(d,tau,reps)
%CONSENSUS      consensus clustering
%
%   CIU = CONSENSUS(D,TAU,REPS) seeks a consensus partition of the 
%   agreement matrix D. The algorithm used here is almost identical to the
%   one introduced in Lancichinetti & Fortunato (2012): The agreement
%   matrix D is thresholded at a level TAU to remove an weak elements. The
%   resulting matrix is then partitions REPS number of times using the
%   Louvain algorithm (in principle, any clustering algorithm that can
%   handle weighted matrixes is a suitable alternative to the Louvain
%   algorithm and can be substituted in its place). This clustering
%   produces a set of partitions from which a new agreement is built. If
%   the partitions have not converged to a single representative partition,
%   the above process repeats itself, starting with the newly built
%   agreement matrix.
%
%   NOTE: In this implementation, the elements of the agreement matrix must
%   be converted into probabilities.
%
%   NOTE: This implementation is slightly different from the original
%   algorithm proposed by Lanchichinetti & Fortunato. In its original
%   version, if the thresholding produces singleton communities, those
%   nodes are reconnected to the network. Here, we leave any singleton
%   communities disconnected.
%
%   Inputs:     D,      agreement matrix with entries between 0 and 1
%                       denoting the probability of finding node i in the
%                       same cluster as node j
%               TAU,    threshold which controls the resolution of the
%                       reclustering
%               REPS,   number of times that the clustering algorithm is
%                       reapplied
%
%   Outputs:    CIU,    consensus partition
%
%   References: Lancichinetti & Fortunato (2012). Consensus clustering in
%   complex networks. Scientific Reports 2, Article number: 336.
%
%   Richard Betzel, Indiana University, 2012
%
%   modified on 3/2014 to include "unique_partitions"

n = length(d); flg = 1;
while flg == 1
    
    flg = 0;
    dt = d.*(d >= tau).*~eye(n);
    if nnz(dt) == 0
        ciu = (1:n)';
    else
        ci = zeros(n,reps);
        for iter = 1:reps
            ci(:,iter) = community_louvain(dt);
        end
        ci = relabel_partitions(ci);
        ciu = unique_partitions(ci);
        nu = size(ciu,2);
        if nu > 1
            flg = 1;
            d = agreement(ci)./reps;
        end
    end
    
end

function cinew = relabel_partitions(ci)
[n,m] = size(ci);
cinew = zeros(n,m);
for i = 1:m
    c = ci(:,i);
    d = zeros(size(c));
    count = 0;
    while sum(d ~= 0) < n
        count = count + 1;
        ind = find(c,1,'first');
        tgt = c(ind);
        rep = c == tgt;
        d(rep) = count;
        c(rep) = 0;
    end
    cinew(:,i) = d;
end

function ciu = unique_partitions(ci)
ci = relabel_partitions(ci);
ciu = [];
count = 0;
c = 1:size(ci,2);
while ~isempty(ci)
    count = count + 1;
    tgt = ci(:,1);
    ciu = [ciu,tgt];                %#ok<AGROW>
    dff = sum(abs(bsxfun(@minus,ci,tgt))) == 0;
    ci(:,dff) = [];
    c(dff) = [];
end