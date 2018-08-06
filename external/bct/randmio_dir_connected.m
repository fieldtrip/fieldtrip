function [R,eff] = randmio_dir_connected(R, ITER)
%RANDMIO_DIR_CONNECTED    Random graph with preserved in/out degree distribution
%
%   R = randmio_dir_connected(W, ITER);
%   [R eff] = randmio_dir_connected(W, ITER);
%
%   This function randomizes a directed network, while preserving the in-
%   and out-degree distributions. In weighted networks, the function
%   preserves the out-strength but not the in-strength distributions. The
%   function also ensures that the randomized network maintains
%   connectedness, the ability for every node to reach every other node in
%   the network. The input network for this function must be connected.
%
%   Input:      W,      directed (binary/weighted) connection matrix
%               ITER,   rewiring parameter
%                       (each edge is rewired approximately ITER times)
%
%   Output:     R,      randomized network
%               eff,    number of actual rewirings carried out
%
%   References: Maslov and Sneppen (2002) Science 296:910
%
%
%   2007-2012
%   Mika Rubinov, UNSW
%   Olaf Sporns, IU

%   Modification History:
%   Jun 2007: Original (Mika Rubinov)
%   Mar 2012: Limit number of rewiring attempts, count number of successful
%             rewirings (Olaf Sporns)


n=size(R,1);
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
        end %rewiring condition
        att=att+1;
    end %while not rewired
end %iterations