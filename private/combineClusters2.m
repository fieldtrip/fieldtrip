function cluster = combineClusters2(labelmat, C)

% COMBINECLUSTERS2 is a helper function for FINDCLUSTER. It searches for
% adjacent clusters in neighbouring channels and combines them. This
% version is an alternative to the original m-file, and an alternative to
% the mex-file. In some cases it may be faster, specifically if the spatial
% adjacency matrix is large (size(C,1)>10000), and the number of columns in
% the input data labelmat is low. Note also, that this version
% supports a sparse adjacency matrix.
%
% A mex-file is available for this functionality, as it can take quite long.
% This mex-file is preferred if the spatial dimension is relatively small
% (<1000). 
%
% Note that there is a minute difference in functional behaviour between
% COMBINECLUSTERS and COMBINECLUSTERS2, where the former occasionally
% labels isolated pixels (i.e. singleton 'clusters' in a spatial slice of
% the labelmat) as a separate clusters, whereas COMBINECLUSTERS2 merges these
% isolated pixels with the cluster in neighbouring slices (provided the
% said isolated pixels is also suprathrehold in the neighbouring slice).

labelmat = permute(labelmat, [2 1]);

% compute an upper bound on the number of edges
nnb = sum(C,1);
max_edges = nnz(labelmat) * max(nnb);

% preallocate some variables
edges1 = zeros(max_edges, 1);
edges2 = zeros(max_edges, 1);
edge_count = 0;

% build a [edges1 edges2] matrix, reflecting the spatially connected non-zero elements of labelmat
for i = 1:size(C,2)
  neighbours = [i find(C(:,i))'];

  mask = sum(labelmat(:,i)~=0&labelmat(:,neighbours)~=0,2)>0;
  a = labelmat(mask, i)*ones(1,nnb(i)+1);
  b = labelmat(mask, neighbours);
  n = nnz(b);

  edges1(edge_count+(1:n)) = a(b>0);
  edges2(edge_count+(1:n)) = b(b>0);
  edge_count = edge_count+n;
end

if isempty(edges1)
  cluster = labelmat;  % nothing to merge
  return
end

% build graph + get the connected components
G = graph(edges1(1:edge_count), edges2(1:edge_count));
comp = conncomp(G);  % 1..N components

% relabel
cluster = zeros(size(labelmat));
mask = labelmat > 0;
cluster(mask) = comp(labelmat(mask));

cluster = permute(cluster, [2 1]);
