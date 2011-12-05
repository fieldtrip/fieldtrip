function [I Q F]=motif3struct_wei(W)
%MOTIF3STRUCT_WEI       Intensity and coherence of structural class-3 motifs
%
%   [I Q F] = motif3struct_wei(W);
%
%   Structural motifs are patterns of local connectivity. Motif frequency
%   is the frequency of occurrence of motifs around a node. Motif intensity
%   and coherence are weighted generalizations of motif frequency. 
%
%   Input:      W,      weighted directed connection matrix
%                       (all weights must be between 0 and 1)
%
%   Output:     I,      motif intensity matrix
%               Q,      motif coherence matrix
%               F,      morif frequency matrix
%
%   Note: Average intensity and coherence are given by I./F and Q./F.
%
%   Reference: Onnela et al. (2005), Phys Rev E 71:065103
%
%
%   Mika Rubinov, UNSW, 2007-2010


persistent M3 M3n ID3 N3
if isempty(N3)
    load motif34lib M3 M3n ID3 N3         	%load motif data
end

n=length(W);                                %number of vertices in W
I=zeros(13,n);                              %intensity
Q=zeros(13,n);                              %coherence
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
            w=[W(v1,u) W(v2,u) W(u,v1) W(v2,v1) W(u,v2) W(v1,v2)];
            s=uint32(sum(10.^(5:-1:0).*[A(v1,u) A(v2,u) A(u,v1)...
                A(v2,v1) A(u,v2) A(v1,v2)]));
            ind=(s==M3n);

            M=w.*M3(ind,:);
            id=ID3(ind);
            l=N3(ind);
            x=sum(M,2)/l;                	%arithmetic mean
            M(M==0)=1;                      %enable geometric mean
            i=prod(M,2)^(1/l);              %intensity
            q=i/x;                          %coherence

            %then add to cumulative count
            I(id,[u v1 v2])=I(id,[u v1 v2])+[i i i];
            Q(id,[u v1 v2])=Q(id,[u v1 v2])+[q q q];
            F(id,[u v1 v2])=F(id,[u v1 v2])+[1 1 1];
        end
    end
end