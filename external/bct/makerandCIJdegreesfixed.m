function [cij,flag] = makerandCIJdegreesfixed(in,out)
%MAKERANDCIJDEGREESFIXED        Synthetic directed random network
%
%   CIJ = makerandCIJdegreesfixed(N,K);
%
%   This function generates a directed random network with a specified 
%   in-degree and out-degree sequence. The function returns a flag, 
%   denoting whether the algorithm succeeded or failed.
%
%   Inputs:     in,     indegree vector
%               out,    outdegree vector
%
%   Output:     CIJ,    binary directed connectivity matrix
%               flag,   flag=1 if the algorithm succeeded; flag=0 otherwise
%
%
%   Notes:  Necessary conditions include:
%               length(in) = length(out) = n
%               sum(in) = sum(out) = k
%               in(i), out(i) < n-1
%               in(i) + out(j) < n+2
%               in(i) + out(i) < n
%
%           No connections are placed on the main diagonal
%
%
% Aviad Rubinstein, Indiana University 2005/2007

% intialize
n = length(in);
k = sum(in);
inInv = zeros(k,1);
outInv = inInv;
iIn = 1; iOut = 1;

for i = 1:n
    inInv(iIn:iIn+in(i) - 1) = i;
    outInv(iOut:iOut+out(i) - 1) = i;
    iIn = iIn+in(i);
    iOut = iOut+out(i);
end

cij = eye(n);
edges = [outInv(1:k)'; inInv(randperm(k))'];

% create cij, and check for double edges and self-connections
for i = 1:k
    if cij(edges(1,i),edges(2,i)),
        warningCounter = 1;
        while (1)
            switchTo = ceil(k*rand);
            if ~(cij(edges(1,i),edges(2,switchTo)) || cij(edges(1,switchTo),edges(2,i))),
                cij(edges(1,i),edges(2,switchTo)) = 1;
                if switchTo < i,
                    cij(edges(1,switchTo),edges(2,switchTo)) = 0;
                    cij(edges(1,switchTo),edges(2,i)) = 1;
                end
                temp = edges(2,i);
                edges(2,i) = edges(2,switchTo);
                edges(2,switchTo) = temp;
                break
            end
            warningCounter = warningCounter+1;
            % If there is a legitimate subtitution, it has a probability of 1/k of being done.
            % Thus it is highly unlikely that it will not be done after 2*k^2 attempts.
            % This is an indication that the given indegree / outdegree
            % vectors may not be possible.
            if warningCounter == 2*k^2
                flag = 0;  % no valid solution found
                return;
            end
        end
    else
        cij(edges(1,i),edges(2,i)) = 1;
    end
end

cij = cij - eye(n);

% a valid solution was found
flag = 1;