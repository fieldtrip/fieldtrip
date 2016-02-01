function  [Erange,eta,Eshort,fs] = erange(CIJ)
%ERANGE     Shortcuts
%
%   [Erange,eta,Eshort,fs] = erange(CIJ);
%
%   Shorcuts are central edges which significantly reduce the
%   characteristic path length in the network.
%
%   Input:      CIJ,        binary directed connection matrix
%
%   Outputs:    Erange,     range for each edge, i.e. the length of the 
%                           shortest path from i to j for edge c(i,j) AFTER
%                           the edge has been removed from the graph.
%               eta         average range for entire graph.
%               Eshort      entries are ones for shortcut edges.
%               fs          fraction of shortcuts in the graph.
%
%   Follows the treatment of 'shortcuts' by Duncan Watts
%
%
%   Olaf Sporns, Indiana University, 2002/2007/2008


N = size(CIJ,1);
K = length(nonzeros(CIJ));
Erange = zeros(N,N);
[i,j] = find(CIJ==1);

for c=1:length(i)
   CIJcut = CIJ;
   CIJcut(i(c),j(c)) = 0;
   [R,D] = reachdist(CIJcut);                   %#ok<ASGLU>
   Erange(i(c),j(c)) = D(i(c),j(c));
end;

% average range (ignore Inf)
eta = sum(Erange((Erange>0)&(Erange<Inf)))/length(Erange((Erange>0)&(Erange<Inf)));

% Original entries of D are ones, thus entries of Erange 
% must be two or greater.
% If Erange(i,j) > 2, then the edge is a shortcut.
% 'fshort' is the fraction of shortcuts over the entire graph.

Eshort = Erange>2;
fs = length(nonzeros(Eshort))/K;