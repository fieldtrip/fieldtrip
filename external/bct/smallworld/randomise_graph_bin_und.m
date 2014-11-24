function C_rand = randomise_graph_bin_und(C)

% randomises the edges of a binary undirected graph. Diagonal elements are
% ignored.
%
% Mark Drakesmith (Cardiff, University)

n=size(C);
C=C~=0;
sort_idx=find(triu(ones(n),1));
C_rand=zeros(n);
[dum rand_idx]=sort(rand(1,length(sort_idx)));

C_rand(sort_idx(rand_idx));

C_rand=C_rand+C_rand';


