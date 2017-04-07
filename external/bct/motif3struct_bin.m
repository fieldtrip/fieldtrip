function [f,F]=motif3struct_bin(A)
%MOTIF3STRUCT_BIN       Frequency of structural class-3 motifs
%
%   [f,F] = motif3struct_bin(A);
%
%   Structural motifs are patterns of local connectivity in complex
%   networks. Such patterns are particularly diverse in directed networks.
%   The motif frequency of occurrence around an individual node is known as
%   the motif fingerprint of that node. The total motif frequency of
%   occurrence in the whole network is correspondingly known as the
%   motif fingerprint of the network. 
%
%   Input:      A,      binary directed connection matrix
%
%   Output:     F,      node motif frequency fingerprint
%               f,      network motif frequency fingerprint
%
%   Note: The function find_motif34.m outputs the motif legend.
%
%   References: Milo et al. (2002) Science 298:824-827
%               Sporns O, Kötter R (2004) PLoS Biol 2: e369
%
%
%   Mika Rubinov, UNSW/U Cambridge, 2007-2015

%   Modification History:
%   2007: Original
%   2015: Improved documentation

persistent M3n ID3
if isempty(ID3)
    load motif34lib M3n ID3              	%load motif data
end

n=length(A);                                %number of vertices in A
F=zeros(13,n);                              %motif count of each vertex
f=zeros(13,1);                              %motif count for whole graph
As=A|A.';                                   %symmetrized adjacency matrix


for u=1:n-2                               	%loop u 1:n-2
    V1=[false(1,u) As(u,u+1:n)];         	%v1: neibs of u (>u)
    for v1=find(V1)
        V2=[false(1,u) As(v1,u+1:n)];       %v2: all neibs of v1 (>u)
        V2(V1)=0;                           %not already in V1
        V2=([false(1,v1) As(u,v1+1:n)])|V2; %and all neibs of u (>v1)
        for v2=find(V2)

            s=uint32(sum(10.^(5:-1:0).*[A(v1,u) A(v2,u) A(u,v1)...
                A(v2,v1) A(u,v2) A(v1,v2)]));
            ind=ID3(s==M3n);
            if nargout==2; F(ind,[u v1 v2])=F(ind,[u v1 v2])+1; end
            f(ind)=f(ind)+1;
        end
    end
end