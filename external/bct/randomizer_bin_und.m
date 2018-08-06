function [R] = randomizer_bin_und(R,alpha)
%RANDOMIZER_BIN_UND     Random graph with preserved in/out degree distribution
%
%   R = randomizer_bin_und(A,alpha);
%
%   This function randomizes a binary undirected network, while preserving 
%   the degree distribution. The function directly searches for rewirable 
%   edge pairs (rather than trying to rewire edge pairs at random), and 
%   hence avoids long loops and works especially well in dense matrices.
%
%   Inputs:     A,          binary undirected connection matrix
%               alpha,      fraction of edges to rewire
%
%   Outputs:    R,          randomized network
%
%   References: Maslov and Sneppen (2002) Science 296:910
%
%
%   Jonathan Power, WUSTL. 3/1/10.

%#ok<*ASGLU>

% make binary
R=ceil(R);

% ensure that matrix is binary
if (max(R(:))~=1) || (min(R(:))~=0)
    error('Matrix should be binary');
end

% ensure that matrix is undirected
if ~isequal(R,R.')
    error('Matrix should be undirected');
end

% find how many edges are possible in the network
[a,b]=size(R);
numpossibleedges=((a*a)-a)/2;

% excise the diagonal and replace it with 9999
savediag=R.*(eye(size(R,1)));
R=R.*(~eye(size(R,1)));
R=R+(eye(size(R,1)).*9999);

% if there are more edges than non-edges we invert the matrix to reduce
% computation time, then revert at the end of the script
inverted=0;
[i,j]=find(triu(R,1)); 
K=size(i,1);
if K>(numpossibleedges/2)
    inverted=1;
    R=double(~R);
    R=R.*(~eye(size(R,1)));
    R=R+(eye(size(R,1)).*9999);
end

% find edges
[i,j]=find(triu(R,1));
K=size(i,1);

% exclude fully connected nodes. will replace later
fullnode=find((sum(triu(R,1),1)+(sum(triu(R,1),2))')==(a-1));
if ~isempty(fullnode)
    R(fullnode,:)=0; R(:,fullnode)=0;
    R=R.*(~eye(size(R,1)));
    R=R+(eye(size(R,1)).*9999);
end

% find the edges
[i,j]=find(triu(R,1));
K=size(i,1);

if (isempty(K) || K==(numpossibleedges) || (K==numpossibleedges-1))
    fprintf('No possible randomization.\n')
else
    for iter=1:K % for every edge
        if rand<=alpha % rewire ~alpha% of edges

            % this is the chosen edge
            a=i(iter);
            b=j(iter);
            
            % for selected edge, see where each end can connect to
            alliholes=find(R(:,i(iter))==0);
            alljholes=find(R(:,j(iter))==0);
            
            % we can only use edges with connection to neither node
            iintersect=intersect(alliholes,alljholes);            
            
            % find which of these nodes are connected
            [ii,jj]=find(R(iintersect,iintersect)==1);         

            % if there an edge to switch
            if ~isempty(ii)

                % choose one randomly
                nummates=size(ii,1);
                mate=ceil(rand*nummates);
                
                % randomly orient the second edge
                if rand<0.5
                    c=iintersect(ii(mate));
                    d=iintersect(jj(mate));
                else
                    d=iintersect(ii(mate));
                    c=iintersect(jj(mate));
                end

                % make the changes in the matrix
                R(a,b)=0; R(c,d)=0;
                R(b,a)=0; R(d,c)=0;
                R(a,c)=1; R(b,d)=1;
                R(c,a)=1; R(d,b)=1;
                
                % update the edge index
                for m=1:K
                    if ((i(m)==d) && (j(m)==c))
                        j(iter)=c; j(m)=b;
                    elseif ((i(m)==c) && (j(m)==d))
                        j(iter)=c; i(m)=b;
                    end
                end
            end % rewiring
        end % if rand<alpha
    end % for every edge
end % if enough edges to flip

% restore full columns
if ~isempty(fullnode)
    R(fullnode,:)=1; R(:,fullnode)=1;
end

% if we did non-edges switch it back to edges
if inverted==1
    R=double(~R);
end

% clear and restore the diagonal
R=R.*(~eye(size(R,1)));
R=R+savediag;
