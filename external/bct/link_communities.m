function M=link_communities(W,type_clustering)
%LINK_COMMUNITIES     Optimal overlapping community structure
%
%   M = link_communities(W)
%   M = link_communities(W,'complete');
%
%   The optimal community structure is a subdivision of the network into
%   groups of nodes which have a high number of within-group connections
%   and a low number of between group connections. 
%
%   This algorithm uncovers overlapping community structure via
%   hierarchical clustering of network links. This algorith is generalized
%   for weighted/directed/fully-connected networks.
%
%   Input:      W,                  directed (weighted or binary) connection matrix.
%               type_clustering,    type of hierarchical clustering (optional)
%                                       'single'        single-linkage (default)
%                                       'complete'      complete-linkage
%                                   
%   Output:     M,                  nodal community-affiliation matrix
%                                   binary matrix of size CxN [communities x nodes]
%
%   NB: The algorithm can be slow and memory intensive.
%
%   Reference: Ahn, Bagrow and Lehmann (2010) Nature 466, 761–764.
%
%   Mika Rubinov, U Cambridge, 2014-2015

%% initialize

n=size(W,1);                                                        % number of nodes
W(1:n+1:end)=0;
W=W./max(W(:));                                                     % normalize weights

if ~exist('type_clustering','var')
    type_clustering='single';
end

%% get node similarity

W(1:n+1:end) = ( sum(W)/sum(W~=0) + sum(W.')/sum(W.'~=0) )/2;       % mean weight on diagonal
No=sum(W.^2,2);                                                     % out-norm squared
Ni=sum(W.^2,1);                                                     % in-norm squared

Jo=zeros(n);                                                        % weighted in-Jaccard
Ji=zeros(n);                                                        % weighted ou-Jaccard
for b=1:n
    for c=1:n;       
        Do=W(b,:)*W(c,:).';
        Jo(b,c)=Do./(No(b)+No(c)-Do);
        
        Di=W(:,b).'*W(:,c);
        Ji(b,c)=Di./(Ni(b)+Ni(c)-Di);
    end
end

%% get link similarity

[A,B]=find( (W|W.') & triu(ones(n),1));
m=length(A);
Ln=zeros(m,2);                                                      % link nodes
Lw=zeros(m,1);                                                      % link weights
for i=1:m;
    Ln(i,:) = [A(i) B(i)];                                          % link nodes
    Lw(i) = (W(A(i),B(i))+W(B(i),A(i)))/2;                          % link weight
end

ES=zeros(m,m,'single');                                             % link similarity
for i=1:m;
    for j=1:m;
        if      Ln(i,1)==Ln(j,1); a=Ln(i,1); b=Ln(i,2); c=Ln(j,2);
        elseif  Ln(i,1)==Ln(j,2); a=Ln(i,1); b=Ln(i,2); c=Ln(j,1);
        elseif  Ln(i,2)==Ln(j,1); a=Ln(i,2); b=Ln(i,1); c=Ln(j,2);
        elseif  Ln(i,2)==Ln(j,2); a=Ln(i,2); b=Ln(i,1); c=Ln(j,1);
        else    continue
        end
        
        ES(i,j) = (W(a,b)*W(a,c)*Ji(b,c) + W(b,a)*W(c,a)*Jo(b,c))/2;
    end
end
ES(1:m+1:end)=0;

%% perform hierarchical clustering

C=zeros(m,m,'single');                                              % community affiliation matrix
Nc=C; Mc=C; Dc=C;                                                   % communities nodes, links and density
U=1:m;                                                              % initial community assignments
C(1,:)=U;                                                           % as above, in the matrix

for i=1:m-1; fprintf('hierarchy%8d\n',i)                            % hierarchy level
    
    % compute densities
    for j=1:length(U);                                              % loop over communities
        idx = C(i,:)==U(j);                                         % get link indices
        links = sort(Lw(idx));                                      % sort link weights
        nodes = sort(reshape(Ln(idx,:),2*nnz(idx),1));
        nodes = nodes([true;nodes(2:end)~=nodes(1:end-1)]);         % get unique nodes
        
        nc = numel(nodes);                                          % community nodes
        mc = sum(links);                                            % community weights
        min_mc = sum(links(1:nc-1));                                % minimal weight
        dc = (mc - min_mc) / (nc.*(nc-1)/2 - min_mc);               % community density
        
        Nc(i,j)=nc;
        Mc(i,j)=mc;
        Dc(i,j)=dc;
    end
    
    % cluster
    C(i+1,:)=C(i,:);                                                % copy current partition
    [u1,u2]=find(ES(U,U)==max(max(ES(U,U))));                       % on this line MAXs MUST BE MAXs
    
    V=U(unique(sortrows(sort([u1 u2],2)),'rows'));                  % get unique links
    for j=1:size(V,1)
        switch type_clustering
            case 'single';      x = max(ES(V(j,:),:),[],1);         % max -> single linkage
            case 'complete';    x = min(ES(V(j,:),:),[],1);         % min -> complete linkage
            otherwise; error('Unknown clustering type.');
        end
        ES(V(j,:),:) = [x;x];                                       % assign distances to whole clusters
        ES(:,V(j,:)) = [x;x].';
        ES(V(j,1),V(j,1)) = 0;                                      % clear diagonal
        ES(V(j,2),V(j,2)) = 0;                                      % clear diagonal
        
        C(i+1,C(i+1,:)==V(j,2)) = V(j,1);                           % merge communities
        V(V==V(j,2)) = V(j,1);                                      % merge indices
    end
    
    U=unique(C(i+1,:));                                             % get unique communities
    if numel(U)==1
        break;
    end
end

%%

Dc(isnan(Dc))=0;
[~,i]=max(sum(Dc.*Mc,2));                                           % get maximal density

U=unique(C(i,:));                                                   % unique communities
M=zeros(1,n);                                                       % nodal affiliations
for j=1:length(U);
    M(j,unique( Ln(C(i,:)==U(j),:)) )=1;
end
M=M(sum(M,2)>2,:);

% M2=zeros(n);                                                      % two dimensional nodal affiliation
% for i=1:size(M,1);
%     M2=M2+(M(i,:).'*ones(1,n) & ones(n,1)*M(i,:));
% end
