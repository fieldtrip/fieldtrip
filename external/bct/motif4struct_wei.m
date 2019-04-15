function [I,Q,F]=motif4struct_wei(W)
%MOTIF4STRUCT_WEI       Intensity and coherence of structural class-4 motifs
%
%   [I,Q,F] = motif4struct_wei(W);
%
%   Structural motifs are patterns of local connectivity in complex
%   networks. Such patterns are particularly diverse in directed networks.
%   The motif frequency of occurrence around an individual node is known as
%   the motif fingerprint of that node. The motif intensity and coherence
%   are weighted generalizations of the motif frequency. The motif
%   intensity is equivalent to the geometric mean of weights of links
%   comprising each motif. The motif coherence is equivalent to the ratio
%   of geometric and arithmetic means of weights of links comprising each
%   motif.
%
%   Input:      W,      weighted directed connection matrix
%                       (all weights must be between 0 and 1)
%
%   Output:     I,      node motif intensity fingerprint
%               Q,      node motif coherence fingerprint
%               F,      node motif frequency fingerprint
%
%   Notes: 
%       1. The function find_motif34.m outputs the motif legend.
%       2. Average intensity and coherence are given by I./F and Q./F
%       3. All weights must be between 0 and 1. This may be achieved using
%          the weight_conversion.m function, as follows: 
%          W_nrm = weight_conversion(W, 'normalize');
%
%   References: Onnela et al. (2005), Phys Rev E 71:065103
%               Milo et al. (2002) Science 298:824-827
%               Sporns O, Kötter R (2004) PLoS Biol 2: e369%
%
%
%   Mika Rubinov, UNSW/U Cambridge, 2007-2015

%   Modification History:
%   2007: Original
%   2015: Improved documentation

persistent M4 M4n ID4 N4
if isempty(N4)
    load motif34lib M4 M4n ID4 N4           	%load motif data
end

n=length(W);                                    %number of vertices in W
I=zeros(199,n);                                 %intensity
Q=zeros(199,n);                                 %coherence
F=zeros(199,n);                                 %frequency

A=1*(W~=0);                                     %adjacency matrix
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
            V3=V3|([false(1,v2) As(v1,v2+1:n)]); %and all neibs of v1 (>v2)
            V3(V1)=0;                           %not already in V1
            V3=V3|([false(1,vz) As(u,vz+1:n)]); %and all neibs of u (>vz)
            for v3=find(V3)

                w=[W(v1,u) W(v2,u) W(v3,u) W(u,v1) W(v2,v1) W(v3,v1)...
                    W(u,v2) W(v1,v2) W(v3,v2) W(u,v3) W(v1,v3) W(v2,v3)];
                s=uint64(sum(10.^(11:-1:0).*[A(v1,u) A(v2,u) A(v3,u)...
                    A(u,v1) A(v2,v1) A(v3,v1) A(u,v2) A(v1,v2)...
                    A(v3,v2) A(u,v3) A(v1,v3) A(v2,v3)]));
                ind=(s==M4n);

                M=w.*M4(ind,:);
                id=ID4(ind);
                l=N4(ind);
                x=sum(M,2)/l;                   %arithmetic mean
                M(M==0)=1;                      %enable geometric mean
                i=prod(M,2)^(1/l);              %intensity
                q=i/x;                          %coherence

                %then add to cumulative count
                I(id,[u v1 v2 v3])=I(id,[u v1 v2 v3])+[i i i i];
                Q(id,[u v1 v2 v3])=Q(id,[u v1 v2 v3])+[q q q q];
                F(id,[u v1 v2 v3])=F(id,[u v1 v2 v3])+[1 1 1 1];
            end
        end
    end
end