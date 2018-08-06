function [Rlatt,Rrp,ind_rp,eff] = latmio_dir_connected(R,ITER,D)
%LATMIO_DIR_CONNECTED     Lattice with preserved in/out degree distribution
%
%   [Rlatt,Rrp,ind_rp,eff] = latmio_dir_connected(R,ITER,D);
%
%   This function "latticizes" a directed network, while preserving the in-
%   and out-degree distributions. In weighted networks, the function
%   preserves the out-strength but not the in-strength distributions. The 
%   function also ensures that the randomized network maintains
%   connectedness, the ability for every node to reach every other node in
%   the network. The input network for this function must be connected.
%
%   Input:      R,      directed (binary/weighted) connection matrix
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
%   Mika Rubinov, UNSW, 2007-2010
%   Olaf Sporns, Indiana University, 2012

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

[i,j]=find(R);
K=length(i);
ITER=K*ITER;

% maximal number of rewiring attempts per 'iter'
maxAttempts= round(n*K/(n*(n-1)));
% actual number of successful rewirings
eff = 0;

for iter=1:ITER
    att=0;
    while (att<=maxAttempts)                                     %while not rewired
        rewire=1;
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

        %rewiring condition
        if ~(R(a,d) || R(c,b))
            %lattice condition
            if (D(a,b)*R(a,b)+D(c,d)*R(c,d))>=(D(a,d)*R(a,b)+D(c,b)*R(c,d))
                %connectedness condition
                if ~(any([R(a,c) R(d,b) R(d,c)]) && any([R(c,a) R(b,d) R(b,a)]))
                    P=R([a c],:);
                    P(1,b)=0; P(1,d)=1;
                    P(2,d)=0; P(2,b)=1;
                    PN=P;
                    PN(1,a)=1; PN(2,c)=1;

                    while 1
                        P(1,:)=any(R(P(1,:)~=0,:),1);
                        P(2,:)=any(R(P(2,:)~=0,:),1);
                        P=P.*(~PN);
                        PN=PN+P;
                        if ~all(any(P,2))
                            rewire=0;
                            break
                        elseif any(PN(1,[b c])) && any(PN(2,[d a]))
                            break
                        end
                    end
                end %connectedness testing

                if rewire               %reassign edges
                    R(a,d)=R(a,b); R(a,b)=0;
                    R(c,b)=R(c,d); R(c,d)=0;

                    j(e1) = d;          %reassign edge indices
                    j(e2) = b;
                    eff = eff+1;
                    break;
                end %edge reassignment
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