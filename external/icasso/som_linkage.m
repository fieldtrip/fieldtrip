function [Z,order,Md] = som_linkage(sM,varargin)

%SOM_LINKAGE Make a hierarchical linkage of the SOM map units.
%
% [Z,order,Dist] = som_linkage(sM, [[argID,] value, ...])
%  
%  Z = som_linkage(sM);
%  Z = som_linkage(D,'complete');
%  Z = som_linkage(sM,'single','ignore',find(~som_hits(sM,D)));
%  Z = som_linkage(sM,pdist(sM.codebook,'mahal'));
%  som_dendrogram(Z); 
%
%  Input and output arguments ([]'s are optional):
%   sM       (struct) map or data struct to be clustered
%            (matrix) size dlen x dim, a data set: the matrix must not
%                     contain any NaN's!
%   [argID,  (string) See below. The values which are unambiguous can 
%    value]  (varies) be given without the preceeding argID.
%
%   Z        (matrix) size dlen-1 x 3, the linkage info
%                     Z(i,1) and Z(i,2) hold the indeces of clusters 
%                     combined on level i (starting from bottom). The new
%                     cluster has index dlen+i. The initial cluster 
%                     index of each unit is its linear index in the 
%                     original data matrix. Z(i,3) is the distance
%                     between the combined clusters. See LINKAGE
%                     function in the Statistics Toolbox.
%                     The ignored samples are listed at the 
%                     end of the Z-matrix and have Z(*,3) == Inf
%   Dist     (matrix) size dlen x dlen, pairwise distance matrix
%
% Here are the valid argument IDs and corresponding values. The values 
% which are unambiguous (marked with '*') can be given without the
% preceeding argID.
%   'linkage' *(string) the linkage criteria to use: 'single' (the
%                       default), 'average' or 'complete' 
%   'topol'   *(struct) topology struct
%   'connect' *(string) 'neighbors' or 'any' (default), whether the
%                       connections should be allowed only between 
%                       neighbors or between any vectors
%              (matrix) size dlen x dlen indicating the connections
%                       between vectors
%              (scalar) the N-neighborhood upto which the connections
%                       should be formed (implies 'neighbors')
%   'ignore'   (vector) the units/vectors which should be ignored 
%   'dist'     (matrix) size dlen x dlen, pairwise distance matrix to 
%                       be used instead of euclidian distances
%              (vector) as the output of PDIST function
%              (scalar) distance norm to use (euclidian = 2)
%   'mask'     (vector) size dim x 1, the search mask used to 
%                       weight distance calculation. By default 
%                       sM.mask or a vector of ones is used.
%
% Note that together 'connect'='neighbors' and 'ignore' may form
% areas on the map which will never be connected: connections
% across the ignored map units simply do not exist.
%
% See also KMEANS_CLUSTERS, LINKAGE, PDIST, DENDROGRAM. 

% Copyright (c) 2000 by Juha Vesanto
% Contributed to SOM Toolbox on June 16th, 2000 by Juha Vesanto
% http://www.cis.hut.fi/projects/somtoolbox/
 
% Version 2.0beta juuso 160600

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input arguments

% the data
if isstruct(sM), 
  if isfield(sM,'data'), D = sM.data; sT = []; mask = []; 
  else D = sM.codebook; sT = sM.topol; mask = sM.mask;
  end
else
  D = sM; sT = []; mask = []; 
end
[dlen dim] = size(D);
if isempty(mask), mask = ones(dim,1); end
if any(isnan(D(:))), error('Data matrix must not have any NaNs.'); end

% varargin
Md = 2; 
linkage = 'single';
ignore_units = []; 
constrained = 0;
i=1; 
while i<=length(varargin), 
  argok = 1; 
  if ischar(varargin{i}), 
    switch varargin{i}, 
      % argument IDs
     case {'topol','som_topol','sTopol'}, i=i+1; sT = varargin{i};
     case 'connect', i=i+1; 
       if ischar(varargin{i}), constrained = ~strcmp(varargin{i},'any');
       else constrained = varargin{i}; end
     case 'ignore',  i=i+1; ignore_units = varargin{i}; 
     case 'dist',    i=i+1; Md = varargin{i};
     case 'linkage', i=i+1; linkage = varargin{i};
     case 'mask',    i=i+1; mask = varargin{i};
     case 'tracking',i=i+1; tracking = varargin{i}; 
      % unambiguous values
     case 'neighbors', constrained = 1; 
     case 'any',       constrained = 0; 
     case {'single','average','complete'}, linkage = varargin{i};
     otherwise argok=0; 
    end
  elseif isstruct(varargin{i}) & isfield(varargin{i},'type'), 
    switch varargin{i}(1).type, 
     case 'som_topol', sTopol = varargin{i}; 
     otherwise argok=0; 
    end
  else
    argok = 0; 
  end
  if ~argok, 
    disp(['(som_linkage) Ignoring invalid argument #' num2str(i+1)]); 
  end
  i = i+1; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% distance matrix

% given distance matrix % jh corrected this place totally 27.3. 03
if (prod(size(Md))==1), % no explicit distance matrix, set flag
  q=2; % 17.2.03 kr added some brackets
else
  if (prod(size(Md))<dlen^2), % check pdist form
    Md = squareform(Md);   % transform to ordinary square diastance matrix
  end
  % jh: 27.3. 03 "calculate pairwise dist. matrix" see approx. 20 lines below
  % sets self-distance to Inf! This must be set here also,
  % otherwise clustering fails for user defined dist. matrix! 
  Md(eye(dlen)==1)=Inf; 
end

% neighborhood constraint
if length(constrained)==1 & constrained>0, 
  Ne1 = som_unit_neighs(sT); 
  Conn = som_neighborhood(Ne1,constrained); 
  Conn(~isfinite(Conn(:))) = 0; 
else Conn = constrained; end
if ~isempty(Conn), for i=1:dlen, C(i,i) = 1; end, end

% pairwise distance matrix across connected units
n = size(D,1);
if prod(size(Md))>1,   
  % remove distances between non-neighbors
  if constrained, for i = 1:n, Md(i,find(Conn(i,:)==0)) = Inf; end, end
else    
  % calculate pairwise distance matrix
  q = Md; 
  Md = zeros(n,n)+Inf;
  if ~constrained & q==2, % fast for the default case 
    for i = 1:n-1,
      x = D(i,:);
      inds = [(i+1):n]; 
      Diff = D(inds,:) - x(ones(n-i,1),:);
      Md(inds,i) = sqrt((Diff.^2)*mask);
      Md(i,inds) = Md(inds,i)';
    end  
  else
    for i = 1:n-1, 
      inds = find(Conn(i,:)==1); 
      inds = inds(find(inds>i)); 
      Diff = abs(D(inds,:) - D(i*ones(length(inds),1),:));
      switch q, 
       case 1,    dist = Diff*mask;
       case 2,    dist = sqrt((Diff.^2)*mask);
       case Inf,  dist = max(Diff,[],2);
       otherwise, dist = ((Diff.^q)*mask).^(1/q);
      end
      Md(inds,i) = dist; 
      Md(i,inds) = dist'; 
    end
  end
end

% set distances to ignored units to Inf
if ~isempty(ignore_units), 
  Md(ignore_units,:) = Inf;
  Md(:,ignore_units) = Inf;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% construct dendrogram

Z = zeros(n-1,3)+NaN;  % merged clusters and distance for each step
clusters = 1:dlen;     % each vector is at first in its own cluster
Cd = Md;               % distances between clusters

h = waitbar(0,'Constructing hierarchical clustering'); 
  
for i=1:n-1,   

  % tracking
  waitbar(i/(n-1),h); 

  %% combine two closest clusters  
  % find the clusters which are closest to each other (c1 and c2)
  [d,ind] = min(min(Cd));  
  if ~isfinite(d), break; end   % no more connected clusters
  [d,c1] = min(Cd(:,ind));      % cluster1
  c2 = clusters(ind);           % cluster2    
  % combine the two clusters
  c1_inds = find(clusters==c1); % vectors belonging to c1
  c2_inds = find(clusters==c2); % vectors belonging to c2
  c_inds = [c1_inds, c2_inds];  % members of the new cluster  
  % new cluster index = bigger cluster
  if length(c2_inds)>length(c1_inds), c=c2; k=c1; else c=c1; k=c2; end
  clusters(c_inds) = c;         % update cluster info
  Z(i,:) = [c, k, d];           % save info into Z
  
  %% update cluster distances  
  % remove the subclusters from the Cd table  
  Cd(c_inds,c_inds) = Inf;      % distance of the cluster to its members = Inf
  k_inds = c_inds(c_inds ~= c); % vectors of the smaller cluster
  Cd(k_inds,:) = Inf;           % distance of the subclusters to 
  Cd(:,k_inds) = Inf;           % other clusters = Inf
  % update the distance of this cluster to the other clusters
  cl = unique(clusters(clusters ~= c)); % indeces of all other clusters
  if ~isempty(cl), % added since v6.5 works differently than 6.1
    for l=cl,   
      o_inds = find(clusters==l); % indeces belonging to cluster k
      vd = Md(c_inds,o_inds);     % distances between vectors in c and k
      vd = vd(isfinite(vd(:)));   % remove infinite distances (no connection)
      len = length(vd);
      if ~len, % if the two clusters are not connected, their distance in Inf
	sd = Inf;
      else   % otherwise, calculate the distance between them
	switch linkage,
	 case 'single',   sd = min(vd);
	 case 'average',  sd = sum(vd)/len; 
	 case 'complete', sd = max(vd);
	 otherwise, error(['Unknown set distance: ' linkage]);
	end
      end
      Cd(c,l) = sd; Cd(l,c) = sd;
    end  
  end
end
close(h); 

last = Z(i,1); 
if isnan(last), 
  last = Z(i-1,1); 
  rest = setdiff(unique(clusters),last); 
  Z(i:n-1,1) = rest'; 
  Z(i:n-1,2) = last; 
  Z(i:n-1,3) = Inf; 
  i = i - 1; 
else 
  rest = []; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% return values

% calculate the order of samples
order = last; 
% go through the combination matrix from top to down
for k=i:-1:1, 
  c = Z(k,1); k = Z(k,2);                     % what (k) change to what (c)
  j = find(order==c);                         % the occurance of c in order    
  if j == length(order), order = [order k];   % put k behind c
  else order = [order(1:j) k order(j+1:end)]; 
  end
end  
order = [rest, order]; 

% to maintain compatibility with Statistics Toolbox, the values in 
% Z must be yet transformed so that they are similar to the output
% of LINKAGE function

Zs = Z;
current_cluster = 1:dlen;
for i=1:size(Z,1),
  Zs(i,1) = current_cluster(Z(i,1));
  Zs(i,2) = current_cluster(Z(i,2));
  current_cluster(Z(i,[1 2])) = dlen+i;  
end

Z = Zs;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




