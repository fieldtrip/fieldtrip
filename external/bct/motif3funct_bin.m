function [f F]=motif3funct_bin(W)
%MOTIF3FUNCT_BIN       Frequency of functional class-3 motifs
%
%   [f F] = motif3funct_bin(A);
%
%   Functional motifs are subsets of connection patterns embedded within 
%   anatomical motifs. Motif frequency is the frequency of occurrence of 
%   motifs around a node.
%
%   Input:      A,      binary directed connection matrix
%
%   Output:     F,      motif frequency matrix
%               f,      motif frequency vector (averaged over all nodes)
%
%
%   Reference: Onnela et al. (2005) Phys Rev E 71:065103
%
%
%   Mika Rubinov, UNSW, 2007-2010


persistent M3 ID3 N3
if isempty(N3)
    load motif34lib M3 ID3 N3            	%load motif data
end

n=length(W);                                %number of vertices in W
f=zeros(13,1);                              %motif count for whole graph
F=zeros(13,n);                          	%frequency

A=1*(W~=0);                                 %adjacency matrix
As=A|A.';                                   %symmetrized adjacency

for u=1:n-2                               	%loop u 1:n-2
    V1=[false(1,u) As(u,u+1:n)];         	%v1: neibs of u (>u)
    for v1=find(V1)
        V2=[false(1,u) As(v1,u+1:n)];       %v2: all neibs of v1 (>u)
        V2(V1)=0;                           %not already in V1
        V2=([false(1,v1) As(u,v1+1:n)])|V2; %and all neibs of u (>v1)
        for v2=find(V2)
            a=[A(v1,u);A(v2,u);A(u,v1);A(v2,v1);A(u,v2);A(v1,v2)];
            ind=(M3*a)==N3;                 %find all contained isomorphs
            id=ID3(ind);

            [idu j]=unique(id);             %unique motif occurences
            j=[0;j];
            mu=length(idu);                 %number of unique motifs
            f2=zeros(mu,1);

            for h=1:mu                      %for each unique motif
                f2(h)=j(h+1)-j(h);              %and frequencies
            end

            %then add to cumulative count
            f(idu)=f(idu)+f2;
            if nargout==2
                F(idu,[u v1 v2])=F(idu,[u v1 v2])+[f2 f2 f2];
            end
        end
    end
end