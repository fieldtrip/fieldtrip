function C_rand = randomise_graph_bin_dir(C)

% randomises the edges of a binary directed graph. Diagonal elements are
% ignored.
%
% Mark Drakesmith (Cardiff, University)

n=size(C);
C=C~=0;
sort_idx=find(ones(n)-eye(n));
C_rand=zeros(n);
[dum rand_idx]=sort(rand(1,length(sort_idx)));

C_rand(sort_idx(rand_idx));
