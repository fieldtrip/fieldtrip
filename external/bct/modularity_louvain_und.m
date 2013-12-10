function [Ci Q]=modularity_louvain_und(W,hierarchy)
%MODULARITY_LOUVAIN_UND     Optimal community structure and modularity
%
%   Ci = modularity_louvain_und(W);
%   [Ci Q] = modularity_louvain_und(W);
%   [Ci_h Q_h] = modularity_louvain_und(W,1);
%
%   The optimal community structure is a subdivision of the network into
%   nonoverlapping groups of nodes in a way that maximizes the number of
%   within-group edges, and minimizes the number of between-group edges. 
%   The modularity is a statistic that quantifies the degree to which the
%   network may be subdivided into such clearly delineated groups. 
%
%   The Louvain algorithm is a fast and accurate community detection 
%   algorithm (as of writing). The algorithm may also be used to detect
%   hierarchical community structure.
%
%   Input:      W       undirected (weighted or binary) connection matrix.
%               h,      optional argument
%                       h=1 enables hierarchical output (see below).
%
%   Outputs:    1. Classic
%                       Ci,     community structure
%                       Q,      modularity
%               2. Hierarchical (if h=1)
%                       Ci_h,   community structure at each hierarchy
%                               (access as Ci_h{1}, Ci_h{2}, ...)
%                       Q_h,    modularity at each hierarhcy
%                               (access as Q_h{1}, Q_h{2}, ...)
%
%   Note: Ci and Q may vary from run to run, due to heuristics in the
%   algorithm. Consequently, it may be worth to compare multiple runs.
%
%   Reference: Blondel et al. (2008)  J. Stat. Mech. P10008.
%
%   Mika Rubinov, UNSW, 2010

%   Modification History:
%   Feb 2010: Original
%   Jun 2010: Fix infinite loops: replace >/< 0 with >/< 1e-10


W=+W;                                       %convert from logical
n=length(W);                                %number of nodes
s=sum(W(:));                                %weight of edges
h=1;                                        %hierarchy index
Ci{h}=1:n;                       	        %hierarchical module assignments
Q{h}=-1;                           	        %hierarchical modularity values
n0=n;                                       %number of nodes

while 1
    K=sum(W);                               %node degree
    Km=K;                                   %module degree
    Knm=W;                                  %node-to-module degree

    M=1:n;                                  %initial module assignments
    Nm=ones(1,n);                           %number of nodes in modules

    flag=true;                              %flag for within-hierarchy search
    while flag
        flag=false;

        for i=randperm(n)                   %loop over all nodes in random order
            dQ=(Knm(i,:)-Knm(i,M(i))+W(i,i)) - K(i).*(Km-Km(M(i))+K(i))/s;
            dQ(M(i))=0;                     %(line above) algorithm condition

            max_dQ=max(dQ);                 %find maximal increase in modularity
            if max_dQ>1e-10;                %if maximal increase is positive
                j=find(dQ==max_dQ,1);

                Knm(:,j)=Knm(:,j)+W(:,i);   %change node-to-module degrees
                Knm(:,M(i))=Knm(:,M(i))-W(:,i);

                Km(j)=Km(j)+K(i);           %change module degrees
                Km(M(i))=Km(M(i))-K(i);

                Nm(j)=Nm(j)+1;              %change number of nodes in modules
                Nm(M(i))=Nm(M(i))-1;

                M(i)=j;                     %reassign module
                flag=true;
            end
        end
    end

    [x x M1]=unique(M);                     %new module assignments (NB: size(M1)=size(M))

    h=h+1;
    Ci{h}=zeros(1,n0);
    for i=1:n                               %loop through initial module assignments
        Ci{h}(Ci{h-1}==i)=M1(i);            %assign new modules
    end

    n=max(M1);                              %new number of modules
    W1=zeros(n);                            %new weighted matrix
    for i=1:n
        for j=i:n
            w=sum(sum(W(M1==i,M1==j)));     %pool weights of nodes in same module
            W1(i,j)=w;
            W1(j,i)=w;
        end
    end
    W=W1;

    Q{h}=sum(diag(W))/s-sum(sum((W/s)^2));  %compute modularity
    if Q{h}-Q{h-1}<1e-10                    %if modularity does not increase
        break
    end
end

Ci([1 end])=[];
Q([1 end])=[];

if nargin==1 || hierarchy==0
    Ci=Ci{end};
    Q=Q{end};
end