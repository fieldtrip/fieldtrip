function seg=bbxflatsegment(node,loop)
%
% seg=bbxflatsegment(node,loop)
%
% decompose edge loops into flat segments along the x/y/z 
% planes of the bounding box
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2008/04/08
%
% input:   
%    node:  x,y,z coordinates of each node of the mesh
%    loop:  input, a single vector separated by NaN, each segment
%             is a close-polygon consisted by node IDs 
% output:
%    seg:   output, a single vector separated by NaN, each segment
%             is a close-polygon on x/y/z plane 
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

pos=node(loop,:);

% get the bounding box
mi=min(pos);
ma=max(pos);

% extract nodes on the bounding box
idx0=find(abs(pos(:,1)-mi(1))<1e-6)';
idx1=find(abs(pos(:,1)-ma(1))<1e-6)';

idy0=find(abs(pos(:,2)-mi(2))<1e-6)';
idy1=find(abs(pos(:,2)-ma(2))<1e-6)';

idz0=find(abs(pos(:,3)-mi(3))<1e-6)';
idz1=find(abs(pos(:,3)-ma(3))<1e-6)';

% need to be more than 3 points to make a flat polygon

if(length(idx0)<=3) idx0=[]; end
if(length(idx1)<=3) idx1=[]; end
if(length(idy0)<=3) idy0=[]; end
if(length(idy1)<=3) idy1=[]; end
if(length(idz0)<=3) idz0=[]; end
if(length(idz1)<=3) idz1=[]; end

nn=length(loop);

% if the original is a flat polygon, return

if(unique(length(idx0))==nn | unique(length(idx1))==nn ...
  |unique(length(idy0))==nn | unique(length(idy1))==nn ...
  |unique(length(idz0))==nn | unique(length(idz1))==nn) 
    seg=loop(:)';
    return;
end

% otherwise, find the combination that split the loop

if(length(unique([idx0 idy0 idz0]))==nn)
   seg= [loop(idx0),nan,loop(idy0),nan,loop(idz0)];
elseif(length(unique([idx0 idy1 idz0]))==nn)
   seg= [loop(idx0),nan,loop(idy1),nan,loop(idz0)];
elseif(length(unique([idx0 idy0 idz1]))==nn)
   seg= [loop(idx0),nan,loop(idy0),nan,loop(idz1)];
elseif(length(unique([idx0 idy1 idz1]))==nn)
   seg= [loop(idx0),nan,loop(idy1),nan,loop(idz1)];
elseif(length(unique([idx1 idy0 idz0]))==nn)
   seg= [loop(idx1),nan,loop(idy0),nan,loop(idz0)];
elseif(length(unique([idx1 idy1 idz0]))==nn)
   seg= [loop(idx1),nan,loop(idy1),nan,loop(idz0)];
elseif(length(unique([idx1 idy0 idz1]))==nn)
   seg= [loop(idx1),nan,loop(idy0),nan,loop(idz1)];
elseif(length(unique([idx1 idy1 idz1]))==nn)
   seg= [loop(idx1),nan,loop(idy1),nan,loop(idz1)];
else
    seg=[];
end

% remove pattern [ ... nan nan ...] in the result

if(length(seg) & any(isnan(seg)))
    id=regexp(sprintf('%d',isnan(seg)),'11');
    if(length(id))
        seg(id+1)=[];
    end
end
