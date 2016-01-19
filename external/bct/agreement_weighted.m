function D = agreement_weighted(CI,Wts)
%WEIGHTED_AGREEMENT     weights agreement matrix
%
%   D = AGREEMENT_WEIGHTED(CI,WTS) is identical to AGREEMENT, with the 
%   exception that each partitions contribution is weighted according to 
%   the corresponding scalar value stored in the vector WTS. As an example,
%   suppose CI contained partitions obtained using some heuristic for 
%   maximizing modularity. A possible choice for WTS might be the Q metric
%   (Newman's modularity score). Such a choice would add more weight to 
%   higher modularity partitions.
%
%   NOTE: Unlike AGREEMENT, this script does not have the input argument
%   BUFFSZ.
%
%   Inputs:     CI,     set of partitions
%               WTS,    relative weight of importance of each paritition
%
%   Outputs:    D,      weighted agreement matrix
%
%   Richard Betzel, Indiana University, 2013

Wts = Wts./sum(Wts);
[N,M] = size(CI);
D = zeros(N);
for i = 1:M
    d = dummyvar(CI(:,i));
    D = D + (d*d')*Wts(i);
end