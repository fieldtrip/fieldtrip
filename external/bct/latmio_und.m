function [Rlatt,Rrp,ind_rp,eff] = latmio_und(R,ITER,D)
%LATMIO_UND     Lattice with preserved degree distribution
%
%   [Rlatt,Rrp,ind_rp,eff] = latmio_und(R,ITER,D);
%
%   This function "latticizes" an undirected network, while preserving the 
%   degree distribution. The function does not preserve the strength 
%   distribution in weighted networks.
%
%   Input:      R,      undirected (binary/weighted) connection matrix
%               ITER,   rewiring parameter
%                       (each edge is rewired approximately ITER times)
%               D,      distance-to-diagonal matrix
%
%   Output:     Rlatt,  latticized network in original node ordering
%               Rrp,    latticized network in node ordering used for
%                       latticization
%               ind_rp, node ordering used for latticization
%               eff,    number of actual rewirings carried out
%
%   References: Maslov and Sneppen (2002) Science 296:910
%               Sporns and Zwi (2004) Neuroinformatics 2:145
%
%   2007-2012
%   Mika Rubinov, UNSW
%   Jonathan Power, WUSTL
%   Olaf Sporns, IU

%   Modification History:
%   Jun 2007: Original (Mika Rubinov)
%   Apr 2008: Edge c-d is flipped with 50% probability, allowing to explore
%             all potential rewirings (Jonathan Power)
%   Feb 2012: limit on number of attempts, distance-to-diagonal as input,
%             count number of successful rewirings (Olaf Sporns)
%   Feb 2012: permute node ordering on each run, to ensure lattices are
%             shuffled across mutliple runs (Olaf Sporns)

n=size(R,1);

% randomly reorder matrix
ind_rp = randperm(n);
R = R(ind_rp,ind_rp);
    
% create 'distance to diagonal' matrix
if nargin<3 %if D is not specified by user
    D=zeros(n);
    u=[0 min([mod(1:n-1,n);mod(n-1:-1:1,n)])];
    for v=1:ceil(n/2)
        D(n-v+1,:)=u([v+1:n 1:v]);
        D(v,:)=D(n-v+1,n:-1:1);
    end
end
%end create

[i,j]=find(tril(R));
K=length(i);
ITER=K*ITER;

% maximal number of rewiring attempts per 'iter'
maxAttempts= round(n*K/(n*(n-1)/2));
% actual number of successful rewirings
eff = 0;

for iter=1:ITER
    att=0;
    while (att<=maxAttempts)                                     %while not rewired
        while 1
            e1=ceil(K*rand);
            e2=ceil(K*rand);
            while (e2==e1),
                e2=ceil(K*rand);
            end
            a=i(e1); b=j(e1);
            c=i(e2); d=j(e2);

            if all(a~=[c d]) && all(b~=[c d]);
                break           %all four vertices must be different
            end
        end

        if rand>0.5
            i(e2)=d; j(e2)=c; 	%flip edge c-d with 50% probability
            c=i(e2); d=j(e2); 	%to explore all potential rewirings
        end
        
        %rewiring condition
        if ~(R(a,d) || R(c,b))
            %lattice condition
            if (D(a,b)*R(a,b)+D(c,d)*R(c,d))>=(D(a,d)*R(a,b)+D(c,b)*R(c,d))
                R(a,d)=R(a,b); R(a,b)=0;
                R(d,a)=R(b,a); R(b,a)=0;
                R(c,b)=R(c,d); R(c,d)=0;
                R(b,c)=R(d,c); R(d,c)=0;

                j(e1) = d;          %reassign edge indices
                j(e2) = b;
                eff = eff+1;
                break;
            end %lattice condition
        end %rewiring condition
        att=att+1;
    end %while not rewired
end %iterations

% lattice in node order used for latticization
Rrp = R;
% reverse random permutation of nodes
[~,ind_rp_reverse] = sort(ind_rp);
Rlatt = Rrp(ind_rp_reverse,ind_rp_reverse);
