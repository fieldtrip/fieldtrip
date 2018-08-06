function [R,eff] = randmio_und_signed(W, ITER)
%RANDMIO_UND     Random graph with preserved signed degree distribution
%
%   R       = randmio_und_signed(W,ITER);
%   [R,eff] = randmio_und_signed(W,ITER);
%
%   This function randomizes an undirected network with positively and
%   negatively signed connections, while preserving the positively and
%   negatively signed degree distribution. The function does not preserve
%   the strength distribution in weighted networks.
%
%   Input:      W,      undirected (binary/weighted) connection matrix
%               ITER,   rewiring parameter
%                       (each edge is rewired approximately ITER times)
%
%   Output:     R,      randomized network
%               eff,    number of actual rewirings carried out
%
%   Reference:  Maslov and Sneppen (2002) Science 296:910
%
%
%   2011-2015
%   Dani Bassett, UCSB
%   Olaf Sporns,  Indiana U
%   Mika Rubinov, U Cambridge

%   Modification History:
%   Mar 2011: Original (Dani Bassett, based on randmio_und.m)
%   Mar 2012: Limit number of rewiring attempts,
%             count number of successful rewirings (Olaf Sporns)
%   Dec 2015: Rewritten the core of the rewiring algorithm to allow
%             unbiased exploration of all network configurations. The new
%             algorithm allows positive-positive/negative-negative
%             rewirings, in addition to the previous positive-positive/0-0
%             and negative-negative/0-0 rewirings (Mika Rubinov). 

if nargin('randperm')==1
    warning('This function requires a recent (>2011) version of MATLAB.')
end

R     = double(W);              % sign function requires double input
n     = size(R,1);
ITER  = ITER*n*(n-1)/2;

% maximal number of rewiring attempts per 'iter'
maxAttempts = round(n/2);
% actual number of successful rewirings
eff = 0;

for iter=1:ITER
    att=0;
    while (att<=maxAttempts)    %while not rewired
        %select four distinct vertices
        nodes = randperm(n,4);
        a = nodes(1);
        b = nodes(2);
        c = nodes(3);
        d = nodes(4);
        
        r0_ab = R(a,b);
        r0_cd = R(c,d);
        r0_ad = R(a,d);
        r0_cb = R(c,b);
        
        %rewiring condition
        if      (sign(r0_ab)==sign(r0_cd)) && ...
                (sign(r0_ad)==sign(r0_cb)) && ...
                (sign(r0_ab)~=sign(r0_ad))
            
            R(a,d)=r0_ab; R(a,b)=r0_ad;
            R(d,a)=r0_ab; R(b,a)=r0_ad;
            R(c,b)=r0_cd; R(c,d)=r0_cb;
            R(b,c)=r0_cd; R(d,c)=r0_cb;
            
            eff = eff+1;
            break;
        end %rewiring condition
        att=att+1;
    end %while not rewired
end %iterations