function [f,F]=motif4funct_bin(A)
%MOTIF4FUNCT_BIN       Frequency of functional class-4 motifs
%
%   [f,F] = motif4funct_bin(A);
%
%   *Structural motifs* are patterns of local connectivity in complex
%   networks. In contrast, *functional motifs* are all possible subsets of
%   patterns of local connectivity embedded within structural motifs. Such
%   patterns are particularly diverse in directed networks. The motif
%   frequency of occurrence around an individual node is known as the motif
%   fingerprint of that node. The total motif frequency of occurrence in
%   the whole network is correspondingly known as the motif fingerprint of
%   the network.
%
%   Input:      A,      binary directed connection matrix
%
%   Output:     F,      node motif frequency fingerprint
%               f,      network motif frequency fingerprint
%
%   Notes: 
%       1. The function find_motif34.m outputs the motif legend.
%       2. There is a source of possible confusion in motif terminology.
%          Motifs ("structural" and "functional") are most frequently
%          considered only in the context of anatomical brain networks
%          (Sporns and Kötter, 2004). On the other hand, motifs are not
%          commonly studied in undirected networks, due to the paucity of
%          local undirected connectivity patterns.
%
%   References: Milo et al. (2002) Science 298:824-827
%               Sporns O, Kötter R (2004) PLoS Biol 2: e369
%
%
%   Mika Rubinov, UNSW/U Cambridge, 2007-2015

%   Modification History:
%   2007: Original
%   2015: Improved documentation

persistent M4 ID4 N4
if isempty(N4)
    load motif34lib M4 ID4 N4                 	%load motif data
end

n=length(A);                                    %number of vertices in A
f=zeros(199,1);
F=zeros(199,n);                                 %frequency

A=1*(A~=0);                                     %adjacency matrix
As=A|A.';                                       %symmetrized adjacency

for u=1:n-3                                     %loop u 1:n-2
    V1=[false(1,u) As(u,u+1:n)];                %v1: neibs of u (>u)
    for v1=find(V1)
        V2=[false(1,u) As(v1,u+1:n)];           %v2: all neibs of v1 (>u)
        V2(V1)=0;                               %not already in V1
        V2=V2|([false(1,v1) As(u,v1+1:n)]);     %and all neibs of u (>v1)
        for v2=find(V2)
            vz=max(v1,v2);                      %vz: largest rank node
            V3=([false(1,u) As(v2,u+1:n)]);     %v3: all neibs of v2 (>u)
            V3(V2)=0;                           %not already in V1&V2
            V3=V3|([false(1,v2) As(v1,v2+1:n)]);%and all neibs of v1 (>v2)
            V3(V1)=0;                           %not already in V1
            V3=V3|([false(1,vz) As(u,vz+1:n)]); %and all neibs of u (>vz)
            for v3=find(V3)

                a=[A(v1,u);A(v2,u);A(v3,u);A(u,v1);A(v2,v1);A(v3,v1);...
                    A(u,v2);A(v1,v2);A(v3,v2);A(u,v3);A(v1,v3);A(v2,v3)];
                ind=(M4*a)==N4;                 %find all contained isomorphs
                id=ID4(ind);

                [idu,j]=unique(id);             %unique motif occurences
                j=[0;j];                        %#ok<AGROW>
                mu=length(idu);                 %number of unique motifs
                f2=zeros(mu,1);

                for h=1:mu                      %for each unique motif
                    f2(h)=j(h+1)-j(h);              %and frequencies
                end

                %then add to cumulative count
                f(idu)=f(idu)+f2;
                if nargout==2
                    F(idu,[u v1 v2 v3])=F(idu,[u v1 v2 v3])+[f2 f2 f2 f2];
                end
            end
        end
    end
end