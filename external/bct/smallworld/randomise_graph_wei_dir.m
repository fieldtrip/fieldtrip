function C_rand = randomise_graph_wei_dir(C)

% randomises the edges of a weighted directed graph. Diagonal elements are
% ignored. The random graph has the same strength distribution as the
% original graph.
%
% Mark Drakesmith (Cardiff, University)

n=size(C,1);

sort_idx=find(ones(n)-eye(n));

C_rand=zeros(n);
[dum rand_idx]=sort(rand(1,length(sort_idx)));


C_rand(sort_idx(rand_idx));

% compute a probability distribtuion
[pdf,wei]=hist(C(sort_idx),n);

% compute cdf
cdf=cumsum(pdf);
cdf=cdf./max(cdf);

% remove duplicades in cdf (interp doesnt like that!)
remove_idx=find(diff(cdf)==0);
cdf(remove_idx)=[];
wei(remove_idx)=[];


% generate random values from uniform distribtuion
rand_uni=(rand(1,length(sort_idx)));

% sample from cdf using interpoation
rand_wei = interp1(cdf,wei,rand_uni);

% put the new edge weights inot the random graph
C_rand(sort_idx(rand_idx))=rand_wei;


