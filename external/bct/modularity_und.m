function [Ci Q]=modularity_und(A)
%MODULARITY_UND     Optimal community structure and modularity
%
%   Ci = modularity_und(W);
%   [Ci Q] = modularity_und(W);
%
%   The optimal community structure is a subdivision of the network into
%   nonoverlapping groups of nodes in a way that maximizes the number of
%   within-group edges, and minimizes the number of between-group edges.
%   The modularity is a statistic that quantifies the degree to which the
%   network may be subdivided into such clearly delineated groups.
%
%   Input:      W,      undirected (weighted or binary) connection matrix.
%
%   Outputs:    Ci,     optimal community structure
%               Q,      maximized modularity
%
%   Note: Ci and Q may vary from run to run, due to heuristics in the
%   algorithm. Consequently, it may be worth to compare multiple runs.
%   Also see Good et al. (2010) Phys. Rev. E 81:046106.
%
%   Reference: Newman (2006) -- Phys Rev E 74:036104, PNAS 23:8577-8582.
%
%
%   2008-2010
%   Mika Rubinov, UNSW
%   Jonathan Power, WUSTL
%   Alexandros Goulas, Maastricht University
%   Dani Bassett, UCSB


%   Modification History:
%   Jul 2008: Original (Mika Rubinov)
%   Oct 2008: Positive eigenvalues made insufficient for division (Jonathan Power)
%   Dec 2008: Fine-tuning made consistent with Newman's description (Jonathan Power)
%   Dec 2008: Fine-tuning vectorized (Mika Rubinov)
%   Sep 2010: Node identities permuted (Dani Bassett)

N=length(A);                            %number of vertices
n_perm = randperm(N);                   %DB: randomly permute order of nodes
A = A(n_perm,n_perm);                   %DB: use permuted matrix for subsequent analysis
K=sum(A);                               %degree
m=sum(K);                               %number of edges
B=A-(K.'*K)/m;                          %modularity matrix
Ci=ones(N,1);                           %community indices
cn=1;                                   %number of communities
U=[1 0];                                %array of unexamined communites

ind=1:N;
Bg=B;
Ng=N;

while U(1)                              %examine community U(1)
    [V D]=eig(Bg);
    [d1 i1]=max(diag(D));               %most positive eigenvalue of Bg
    v1=V(:,i1);                         %corresponding eigenvector

    S=ones(Ng,1);
    S(v1<0)=-1;
    q=S.'*Bg*S;                         %contribution to modularity

    if q>1e-10                       	%contribution positive: U(1) is divisible
        qmax=q;                         %maximal contribution to modularity
        Bg(logical(eye(Ng)))=0;      	%Bg is modified, to enable fine-tuning
        indg=ones(Ng,1);                %array of unmoved indices
        Sit=S;
        while any(indg);                %iterative fine-tuning
            Qit=qmax-4*Sit.*(Bg*Sit); 	%this line is equivalent to:
            qmax=max(Qit.*indg);        %for i=1:Ng
            imax=(Qit==qmax);           %	Sit(i)=-Sit(i);
            Sit(imax)=-Sit(imax);       %	Qit(i)=Sit.'*Bg*Sit;
            indg(imax)=nan;             %	Sit(i)=-Sit(i);
            if qmax>q;                  %end
                q=qmax;
                S=Sit;
            end
        end

        if abs(sum(S))==Ng              %unsuccessful splitting of U(1)
            U(1)=[];
        else
            cn=cn+1;
            Ci(ind(S==1))=U(1);         %split old U(1) into new U(1) and into cn
            Ci(ind(S==-1))=cn;
            U=[cn U];
        end
    else                                %contribution nonpositive: U(1) is indivisible
        U(1)=[];
    end

    ind=find(Ci==U(1));                 %indices of unexamined community U(1)
    bg=B(ind,ind);
    Bg=bg-diag(sum(bg));                %modularity matrix for U(1)
    Ng=length(ind);                     %number of vertices in U(1)
end

s=Ci(:,ones(1,N));                      %compute modularity
Q=~(s-s.').*B/m;
Q=sum(Q(:));
Ci_corrected = zeros(N,1);              % DB: initialize Ci_corrected
Ci_corrected(n_perm) = Ci;              % DB: return order of nodes to the order used at the input stage.
Ci = Ci_corrected;                      % DB: output corrected community assignments